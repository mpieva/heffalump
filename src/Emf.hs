{-# LANGUAGE ScopedTypeVariables #-}
module Emf ( main_emf
           , Genome
           , encodeGenome
           , emptyGenome
           , parseMaf
           ) where

-- ^ Getting variant calls from EMF file (Compara EPO whole genome
-- alignments).  EMF contains a list of blocks, and each block aligns
-- contiguous stretches of genomes; sometimes multiple parts of one
-- genome appear in the same block.
--
-- To scan for variants, we compute a difference between the first
-- column and every other column.  Each column gets mapped to some
-- species:  another ape or one of the human-ape common ancestors.
-- (Note that the human-gorilla ancestor is usually, but not necessarily
-- also the ancestor of the chimp!)  Multiple columns can correspond to
-- the same species, and we only call variants for a given species if at
-- least one column maps to it and all such columns agree.


import Bio.Prelude               hiding ( compl )
import Streaming
import System.Console.GetOpt
import System.Directory                 ( getDirectoryContents )
import System.FilePath
import System.IO

import qualified Data.Attoparsec.ByteString.Char8   as A
import qualified Data.ByteString.Char8              as B
import qualified Data.ByteString.Lazy               as L
import qualified Data.ByteString.Streaming.Char8    as S
import qualified Data.IntMap.Strict                 as I
import qualified Streaming.Prelude                  as Q

import Genome
import Lump
import Util

-- | We get disjointed blocks, and they are not guaranteed to be in order.
-- (There is an important use case where they are guaranteed to be
-- jumbled.)  So we need an intermediate, compact representation of the
-- fragmented genome.  We map chromosome indices to maps from starting
-- position to a pair of a 'PackedLump' and its length.

newtype Genome = Genome (I.IntMap (I.IntMap (Int, PackedLump)))

emptyGenome :: Genome
emptyGenome = Genome I.empty

-- | Inserts an encoded block into a genome.  Blocks are checked for
-- unexpected overlaps and warnings are printed.  In MAF files generated
-- in an unspecified manner from EMF files, I get lots of strange
-- overlaps.  They are "fixed" by ignoring the offending blocks.  If
-- this happens, live with the damage or fix your input files.  In
-- contrast, reading the same EMF files directly, I only get warnings
-- about duplicated blocks.

insertGenome :: Genome -> Int -> Int -> Int -> PackedLump -> (String, Genome)
insertGenome (Genome m1) !c !p !l !s = ( warning, Genome m1' )
  where
    !m2 = I.findWithDefault I.empty c m1
    !m2' = I.insert p (l,s) m2
    !m1' = I.insert c m2' m1

    warning = case I.splitLookup p m2 of
        (_, Just (l',s'), _)
            | l == l' && unpackLump s == unpackLump s'
                -> "\27[32;1mduplication at " ++ shows c ":" ++ shows p "\27[0m"
            | otherwise
                -> "\27[33;1mcollision at " ++ shows c ":" ++ shows p
                    " with different " ++
                    (if l == l' then "sequences" else "blocks") ++ "\27[0m"

        (lm, _, rm) -> case (I.maxViewWithKey lm, I.minViewWithKey rm) of

            (Just ((p1, (l1, _)), _), _)
                | p1 + l1 > p
                    -> overlap p1 (p1+l1) p (p+l)

            (_, Just ((p2, (l2, _)), _))
                | p + l > p2
                    -> overlap p (p+l) p2 (p2+l2)

            _ -> ""

    overlap u1 v1 u2 v2
        = "\27[31;1moverlap between " ++ shows c ":" ++ shows u1 "-" ++ shows v1
          " and " ++ shows c ":" ++ shows u2 "-" ++ shows v2 "\27[0m"


-- | To parse MAF, we shouldn't need a reference, and in fact we don't
-- look at the reference sequence at all.  We use the reference only to
-- map chromsome names to indices.  (We allow UCSC-style names to match
-- NCBI-style names.  This happens all the time for the human genome.)
--
-- | Reads MAF.  We discard the header lines starting with '#'.  The
-- rest must be blocks.  Each block must begin with an "a" line and
-- extends to the next empty line.  The "s" lines in the blocks are
-- scanned for the two expected species, the coordinates are extracted
-- from the line for the reference sequence.  Encoded blocks are
-- accumulated in memory and finally emitted in the correct order.
-- Chromosome names are mapped to target indices according to the
-- reference genome.

parseMaf :: (B.ByteString, B.ByteString) -> [B.ByteString] -> Genome -> S.ByteString IO r -> IO Genome
parseMaf (from_spc, to_spc) chrs gnm0 = doParse gnm0 . mapsM S.toStrict . S.lines
  where
    doParse :: Genome -> Stream (Of B.ByteString) IO r -> IO Genome
    doParse gnm = inspect . Q.dropWhile (\l -> B.all isSpace l || B.head l == '#') >=> \case
        Left _                                      -> return gnm
        Right (ahdr :> lns) | "a":_ <- B.words ahdr -> do
            slns :> lns' <- Q.toList $ Q.span (not . B.all isSpace) lns
            gnm' <- foldl (>>=) (return gnm) $ do
                        ( [from_chr], from_pos, from_len, from_rev, from_seq ) <- extract from_spc slns
                        ( _to_chr,   _to_pos,  _to_len,  _to_rev,     to_seq ) <- extract to_spc   slns
                        let (ref', smp') = B.unzip $ filter ((/=) '-' . fst) $ B.zip from_seq to_seq

                        return $ \g -> do
                            hPrintf stderr "%d:%d-%d\n" from_chr from_pos (from_pos+from_len)
                            insert1 g from_rev from_chr from_pos from_len ref' smp'

            doParse gnm' lns'

        Right (l :> ls) -> do ls' :> _ <- Q.toList $ Q.take 2 ls
                              unexpected (B.unpack (B.unlines $ l : ls'))


    insert1 :: Genome -> Bool -> Int -> Int -> Int -> B.ByteString -> B.ByteString -> IO Genome
    insert1 g False chrm pos len ref smp = do
            plump :> _ <- encodeLumpToMem $ diff ref smp
            let (w,g') = insertGenome g chrm pos len plump
            unless (null w) $ hPutStrLn stderr w
            return g'

    insert1 g True chrm pos len ref smp = do
            plump :> _ <- encodeLumpToMem $ diff (revcompl ref) (revcompl smp)
            let (w,g') = insertGenome g chrm (pos-len+1) len plump
            unless (null w) $ hPutStrLn stderr w
            return g'


    extract :: B.ByteString -> [B.ByteString] -> [( [Int], Int, Int, Bool, B.ByteString )]
    extract spc lns =
        [ ( chrom, pos, len, rev, seq_ )
        | sln <- lns
        , "s":name:pos_:len_:str_:_olen:seq_:_ <- [ B.words sln ]
        , B.snoc spc '.' `B.isPrefixOf` name || spc == name
        , let chr_ = B.drop (fromIntegral $ B.length spc + 1) name
        , let chrom = take 1 $                                       elemIndices chr_ chrs
                            ++ do ("chr",cc) <- [B.splitAt 3 chr_] ; elemIndices  cc  chrs
                            ++ do "chrM"     <- [chr_]             ; elemIndices "MT" chrs
        , (pos, "") <- maybeToList $ B.readInt pos_
        , (len, "") <- maybeToList $ B.readInt len_
        , rev <- [ False | str_ == "+" ] ++ [ True | str_ == "-" ] ]


encodeGenome :: Genome -> L.ByteString
encodeGenome (Genome m1)
    | I.null m1 = error "effectively empty genome"
    | otherwise = L.concat [ maybe (unpackLump noLump) encodeChr $ I.lookup i m1
                           | i <- [ 0 .. fst (I.findMax m1) ] ]
  where
    encodeChr :: I.IntMap (Int, PackedLump) -> L.ByteString
    encodeChr m2 = I.foldrWithKey step (const L.empty) m2 0

    step :: Int -> (Int, PackedLump) -> (Int -> L.ByteString) -> Int -> L.ByteString
    step start (len, PackLump lump) k pos
        -- Could we get duplicated or overlapping blocks?  We won't deal
        -- with them properly, but we skip them just in case, so we get
        -- valid, if not correct output.
        | start >= pos = lots_of_Ns (start - pos) <> lump <> k (start+len)
        | otherwise    =                                     k start

    -- Argh, this wasn't supposed to be repeated here, but it's the most
    -- straight forward way to do it.  If the encoding of 'Lump's (see
    -- 'encodeLump') ever changes, this needs to be adapted.
    lots_of_Ns :: Int -> L.ByteString
    lots_of_Ns n
        | n == 0        = L.empty
        | n < 0x3C      = L.singleton $ 0x40 .|. fromIntegral n
        | n < 0x100     = L.pack [ 0x40 .|. 0x3F
                                 , fromIntegral  n ]
        | n < 0x10000   = L.pack [ 0x40 .|. 0x3E
                                 , fromIntegral (n             .&. 0xff)
                                 , fromIntegral (n `shiftR`  8 .&. 0xff) ]
        | n < 0x1000000 = L.pack [ 0x40 .|. 0x3D
                                 , fromIntegral (n             .&. 0xff)
                                 , fromIntegral (n `shiftR`  8 .&. 0xff)
                                 , fromIntegral (n `shiftR` 16 .&. 0xff) ]
        | otherwise     = L.pack [ 0x40 .|. 0x3C
                                 , fromIntegral (n             .&. 0xff)
                                 , fromIntegral (n `shiftR`  8 .&. 0xff)
                                 , fromIntegral (n `shiftR` 16 .&. 0xff)
                                 , fromIntegral (n `shiftR` 24 .&. 0xff) ]


data EmfBlock = EmfBlock { emf_tree :: Tree Label
                         , emf_seqs :: [B.ByteString] }
    deriving Show

scanBlocks :: Monad m => Stream (Of B.ByteString) m r -> Stream (Of EmfBlock) m r
scanBlocks inp = do
        rngs :> inp1 <- lift . Q.toList .
                        Q.zipWith parse_seq_name (Q.enumFrom 0) .
                        Q.span (B.isPrefixOf "SEQ ") .              -- names and indices of sequences
                        Q.dropWhile (B.all A.isSpace) .             -- and empty lines
                        Q.dropWhile ("#" `B.isPrefixOf`) $ inp      -- get rid of comments

        lift (Q.next inp1) >>= \case
            -- We detect EOF here.  ('rngs' will be empty, too, which is not a problem.)
            Left r -> pure r
            Right (treeln, inp1a) -> do
                Right ("DATA", inp2) <- lift $ Q.next inp1a
                let tree = relabel rngs $ parse_tree treeln         -- topology of tree w/ indices
                datalns :> remainder <- lift $ Q.toList $ Q.map (B.map toUpper) $ Q.span ("//" /=) inp2
                EmfBlock tree (B.transpose datalns) `Q.cons` scanBlocks (Q.drop 1 remainder)
  where
    parse_seq_name i s = Label i (spc, B.drop 1 race) (str $ Range (Pos sqn (beg-1)) (end-beg+1))
      where
        ("SEQ":sr:sqn:beg':end':str':_) = B.words s
        (spc,race) = B.break ('_' ==) sr
        Just (beg,"") = B.readInt beg'
        Just (end,"") = B.readInt end'
        str = if str' == "1" then id else reverseRange

    parse_tree s = case A.parseOnly (A.string "TREE " *> treep <* A.char ';') s of
        Left  e -> error $ "parse of TREE line failed with " ++ e ++ " on " ++ show s
        Right r -> r

    treep =   Branch <$> (A.char '(' *> treep) <*> (A.char ',' *> treep) <*> (A.char ')' *> label)
            <|> Leaf <$> label

    label = (\x y z -> (x `B.append` B.pack ['[',y,']'], z))
            <$> A.takeWhile (\c -> A.isAlpha_ascii c || A.isDigit c || c == '_' || c == '.')
            <*> (A.char '[' *> (A.char '+' <|> A.char '-'))
            <*> (A.char ']' *> A.char ':' *> A.double)


data Label = Label Int Species Range deriving Show
data Tree label = Branch (Tree label) (Tree label) label | Leaf label deriving Show

relabel :: [Label] -> Tree (B.ByteString, Double) -> Tree Label
relabel ls (Leaf       (lbl,_)) = Leaf                                 (find_label lbl ls)
relabel ls (Branch u v (lbl,_)) = Branch (relabel ls u) (relabel ls v) (find_label lbl ls)

find_label :: B.ByteString -> [Label] -> Label
find_label lbl ls = case filter ((==) lbl . get_short) ls of
    [l] -> l
    []  -> error $ "not found: " ++ show lbl ++ " in " ++ show ls
    _   -> unexpected (show lbl ++ " in " ++ show ls)

get_short :: Label -> B.ByteString
get_short (Label _ (spc,race) (Range (Pos chrom pos) len)) = B.concat $
    [ B.map toUpper (B.take 1 spc),
      B.take 3 race, B.singleton '_',
      chrom, B.singleton '_' ] ++
    if pos >= 0 then [ B.pack (show $ 1+pos), B.singleton '_',
                       B.pack (show $ pos+len), "[+]" ]
                else [ B.pack (show $ 1-pos-len), B.singleton '_',
                       B.pack (show $ -pos), "[-]" ]

type Species = (B.ByteString,B.ByteString)

-- | Find closely related sequence in another species.  We look for a
-- tree that has homo on one side, but no chimp, and at least one chimp
-- on the other.  For each find of home, we generate one result, but all
-- chimps in the other branch are collected.
species_to_species :: Species -> Species -> Tree Label -> [(Label, [Label])]
species_to_species homo pan = subtrees
  where
    subtrees (Leaf       _) = []
    subtrees (Branch u v _) = case ( tree_has homo u, tree_has pan u, tree_has homo v, tree_has pan v ) of
        (   hs,    [   ],    _, cs@(_:_) ) -> [ (h,cs) | h <- hs ] ++ subtrees v
        (    _, cs@(_:_),   hs,    [   ] ) -> [ (h,cs) | h <- hs ] ++ subtrees u
        (    _,        _,    _,        _ ) -> subtrees u ++ subtrees v

-- Find common ancestor.  The common ancestor (with, say, chimp) is the
-- root of a tree that has homo on one side but no chimp, and chimp on
-- the other.  The ancestral sequence may well apply to more than one
-- coordinate on homo, and we may find more than one suitable tree.
--
-- We return two labels, the first is the human sequence, the second the
-- ancestor.  Multiple hits to different subtrees are possible.

ancestor_to_species :: Species -> Species -> Tree Label -> [(Label, [Label])]
ancestor_to_species homo pan = subtrees
  where
    subtrees (Leaf          _) = []
    subtrees (Branch u v albl) = case ( tree_has homo u, tree_has pan u, tree_has homo v, tree_has pan v ) of
        (   hs,  [],    _, _:_ ) -> [ (h,[albl]) | h <- hs ] ++ subtrees v
        (    _, _:_,   hs,  [] ) -> [ (h,[albl]) | h <- hs ] ++ subtrees u
        (    _,   _,    _,   _ ) -> subtrees u ++ subtrees v

tree_has :: Species -> Tree Label -> [Label]
tree_has s (Leaf lbl@(Label _ spc _)) = [ lbl | s == spc ]
tree_has s (Branch           u' v' _) = tree_has s u' ++ tree_has s v'


-- | Returns a 'Fragment' for a subset of the species (say, the Chimps)
-- in a tree relative to a specific one (usually the Human).  If we find
-- multiple Human sequences, we return mutiple results.  If we find
-- multiple Chimps, we turn them into a consensus sequence and return
-- the result for the consensus sequence.
collectTrees :: [B.ByteString] -> ( Tree Label -> [(Label, [Label])] ) -> Genome -> EmfBlock -> IO Genome
collectTrees chrs select genome block
    = foldlM (\g f -> f g) genome
            [ \g -> do plump :> _ <- encodeLumpToMem $ diff ref_seq smp_seq
                       let (w,g') = insertGenome g ref_idx (p_start (r_pos rref)) (r_length rref) plump
                       hPutStrLn stderr $ "Inserted " ++ shows ref_idx ":" ++
                            shows (p_start (r_pos rref)) "-" ++
                            shows (p_start (r_pos rref) + r_length rref) ", " ++
                            shows (L.length (unpackLump plump)) " bytes."
                       unless (null w) $ hPutStrLn stderr $ "Warning: " ++ w
                       return $! g'
            | (Label iref _ rref_, tgt_lbls) <- select $ emf_tree block
            , ref_idx <- elemIndices (p_seq (r_pos rref_)) chrs
            , let ref_seq_ = emf_seqs block !! iref
            , let smp_seq_ = consensus_seq [ emf_seqs block !! itgt | Label itgt _ _ <- tgt_lbls ]
            , let (rref, ref_seq, smp_seq) =
                    if p_start (r_pos rref_) >= 0
                    then (             rref_,          ref_seq_,          smp_seq_)
                    else (reverseRange rref_, revcompl ref_seq_, revcompl smp_seq_) ]

revcompl :: B.ByteString -> B.ByteString
revcompl = B.reverse . B.map compl
  where
    compl 'A' = 'T' ; compl 'T' = 'A' ; compl 'a' = 't' ; compl 't' = 'a'
    compl 'C' = 'G' ; compl 'G' = 'C' ; compl 'c' = 'g' ; compl 'g' = 'c'
    compl  x  =  x


consensus_seq :: [B.ByteString] -> B.ByteString
consensus_seq [    ] = B.empty
consensus_seq [  x ] = x
consensus_seq (x:xs)
    | all (== B.length x) (map B.length xs) =
        B.pack [ foldl' (\r c -> if r == c then r else 'N')
                        (B.index x i) (map (`B.index` i) xs)
               | i <- [0 .. B.length x -1] ]
    | otherwise =
        error "consensus_seq: sequences have uneuqal lengths."


data OptsEmf = OptsEmf {
    emf_output :: FilePath,
    emf_reference :: FilePath,
    emf_ref_species :: Species,
    emf_select :: Species -> Tree Label -> [(Label, [Label])] }

opts_emf_default :: OptsEmf
opts_emf_default = OptsEmf
    (error "no output specified")
    (error "no reference specified")
    ("homo","sapiens")
    (`ancestor_to_species` ("pan","troglodytes"))

opts_emf :: [ OptDescr ( OptsEmf -> IO OptsEmf ) ]
opts_emf =
    [ Option "o" ["output"]         (ReqArg set_output  "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"]      (ReqArg set_ref     "FILE") "Read reference from FILE (.2bit)"
    , Option "R" ["ref-species"]    (ReqArg set_ref_sp  "NAME") "Set reference species to NAME"
    , Option "S" ["sample-species"] (ReqArg set_species "NAME") "Import species NAME"
    , Option "A" ["ancestor-of"]    (ReqArg set_anc     "NAME") "Import ancestor with species NAME" ]
  where
    set_output  a c = return $ c { emf_output = a }
    set_ref     a c = return $ c { emf_reference  = a }
    set_ref_sp  a c = return $ c { emf_ref_species = w2 a }
    set_species a c = return $ c { emf_select = flip species_to_species  (w2 a) }
    set_anc     a c = return $ c { emf_select = flip ancestor_to_species (w2 a) }

    w2 s = case B.split ',' (B.pack (map toLower s)) of a:b:_ -> (a,b) ; a:_ -> (a,B.empty) ; _ -> (B.empty,B.empty)


main_emf :: [String] -> IO ()
main_emf args = do
    ( emfs, OptsEmf{..} ) <- parseFileOpts opts_emf_default (mk_opts "emf" "[emf-file...]" opts_emf) args
    ref <- readTwoBit emf_reference

    let cons = collectTrees (rss_chroms ref) (emf_select emf_ref_species)
    !genome <- foldM (\g fp -> withFile fp ReadMode $
                                    Q.foldM_ cons (return g) return . scanBlocks .
                                    mapsM S.toStrict . S.lines . decomp . S.fromHandle
                     ) emptyGenome
               -- this is just so the one path at MPI that makes sense
               -- doesn't need to be typed over and over again
               =<< case emfs of
                    [] -> let dir = "/mnt/expressions/martin/sequence_db/epo/epo_6_primate_v66"
                          in map (dir </>) . filter (".emf.gz" `isSuffixOf`) <$> getDirectoryContents dir
                    _  -> return emfs

    withFile emf_output WriteMode $ \hdl ->
            L.hPut hdl $ encodeGenome genome

