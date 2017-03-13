{-# LANGUAGE OverloadedStrings, BangPatterns, ScopedTypeVariables #-}
module Emf where

-- ^ Getting variant calls from EMF file (Compara EPO whole genome
-- alignments).  EMF contains a list of blocks, and each block align
-- contiguous stretches of genomes, sometime multiple parts of one
-- genome appear in the same block.
--
-- To scan for variants, we compute a difference between the first
-- column and every other column in such a format that it can be plugged
-- into VCF.  Each column gets mapped to some species:  another ape or
-- one of the human-ape common ancestors.  (Note that the human-gorilla
-- ancestor is usually, but not necessarily also the ancestor of the
-- chimp!)  Multiple columns can correspond to the same species, and we
-- only call variants for some species if at least column maps to it and
-- all these columns agree.


-- XXX  If our 'Stretch' data type could encode strings of odd length,
-- we could simplify the whole strategy in here:  We could scan the
-- input in its given order, build small streatches, the sort them.  As
-- is, this is awkward.  *sigh*  Maybe in version 2.

import BasePrelude
import Codec.Compression.GZip
import Data.ByteString.Builder          ( hPutBuilder )
import Data.ByteString.Short            ( ShortByteString )
import System.Console.GetOpt
import System.Directory                 ( getDirectoryContents )
import System.FilePath
import System.IO

import qualified Codec.Compression.Snappy           as Snappy
import qualified Data.Attoparsec.ByteString.Char8   as A
import qualified Data.ByteString                    as SB
import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Lazy               as LB
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.ByteString.Short              as BS
import qualified Data.HashMap.Strict                as H
import qualified Data.IntMap.Strict                 as I

import Stretch
import Util

-- Read MAF.  We discard the header lines starting with '#'.  The rest
-- must be alignment blocks.  Each a-block must contain its own header
-- line and exactly two s-lines.  The first line must have human
-- coordinates, and the file must be sorted by human coordinates.   We
-- ignore the chromosome name.
--
-- We wish to skip about between blocks.  But sometimes blocks join
-- immediately, and then we get into trouble because 'diff' only works
-- correctly for stretches of even length.  So instead we return a new
-- data type to make that explicit.

data MafStretch = Skip !Int MafStretch | Aln !L.ByteString !L.ByteString MafStretch | MafDone

parseMaf :: L.ByteString -> Stretch -> Stretch
parseMaf inp done = normMaf done . parse1 0 $ L.lines inp
  where
    parse1 oldpos lns = case dropWhile (\l -> L.null l || L.all isSpace l || L.head l == '#') lns of
        [] -> MafDone
        (ahdr:shum:sother:lns')
            | "a " `L.isPrefixOf` ahdr
            , ["s", _, hpos_, _hlen, "+",  _clen, hum_seq] <- L.words shum
            , ["s", _, _opos, _olen,   _, _oclen, oth_seq] <- L.words sother
            , Just (hpos, "") <- L.readInt hpos_
            , hpos >= oldpos ->
                -- We have to deal with gaps.  Gaps in the reference can
                -- simply be removed.  If the sample is gapped, we want
                -- to encode an N.  Thankfully, that's precisely what
                -- diff already does when it sees a gap.
                let (ref', smp') = LB.unzip $ map (up *** up) $ filter ((/=) 45 . fst) $ LB.zip hum_seq oth_seq
                in Skip (hpos - oldpos) $ Aln ref' smp' $
                   parse1 (fromIntegral (L.length ref') + hpos) lns'

        ls -> error $ "WTF?! (" ++ show oldpos ++ ") near \n" ++ L.unpack (L.unlines $ take 3 ls)

-- remove degenerate stretches
normMaf :: Stretch -> MafStretch -> Stretch
normMaf done = norm
  where
    norm (Skip  0 ms)            = norm ms
    norm (Aln r _ ms) | L.null r = norm ms

    -- join adjacent similar stretches
    norm (Skip n1   (Skip n2   ms)) = norm $ Skip (n1+n2) ms
    norm (Aln r1 s1 (Aln r2 s2 ms)) = norm $ Aln (L.append r1 r2) (L.append s1 s2) ms

    -- emit 'Skip' or 'Aln', but transfer an N to make the length even
    norm (Skip n1 (Aln ref smp ms))
        | L.null ref = norm $ Skip n1 ms
        | even n1    = Ns (fromIntegral $ n1 `div` 2) $ norm $ Aln             ref              smp  ms
        | n1 == 1    =                                  norm $ Aln (L.cons 'X' ref) (L.cons 'N' smp) ms
        | otherwise  = Ns (fromIntegral $ n1 `div` 2) $ norm $ Aln (L.cons 'X' ref) (L.cons 'N' smp) ms

    norm (Aln r1 s1 (Skip n2 ms))
        | n2 == 0            = norm $ Aln r1 s1 ms
        | even (L.length r1) = diff         r1              s1      (norm $ Skip       n2  ms)
        | otherwise          = diff (L.snoc r1 'X') (L.snoc s1 'N') (norm $ Skip (pred n2) ms)

    -- anything else runs into the end
    norm (Aln r1 s1 MafDone) = diff r1 s1 done
    norm (Skip n1   MafDone) = Ns (fromIntegral $ succ n1 `div` 2) done
    norm            MafDone  = done

data EmfBlock = EmfBlock { emf_tree :: Tree Label
                         , emf_seqs :: [S.ByteString] }
    deriving Show

scanOneBlock :: [S.ByteString] -> ( EmfBlock, [S.ByteString] )
scanOneBlock inp = ( EmfBlock tree seqs, drop 1 remainder )
  where
    (seqlns, treeln : "DATA" : inp2) = span (S.isPrefixOf "SEQ ") .           -- names and indices of sequences
                                       dropWhile (S.all A.isSpace) $          -- and empty lines
                                       dropWhile ("#" `S.isPrefixOf`) inp     -- get rid of comments

    rngs = zipWith parse_seq_name [0..] seqlns
    tree = relabel rngs $ parse_tree treeln                                   -- topology of tree w/ indices

    (datalns, remainder) = span ("//" /=) inp2
    seqs = S.transpose . map (S.map toUpper) $ datalns

    parse_seq_name i s = Label i (spc, S.drop 1 race) (str $ Range (BS.toShort sqn) (beg-1) (end-beg+1))
      where
        ("SEQ":sr:sqn:beg':end':str':_) = S.words s
        (spc,race) = S.break ('_' ==) sr
        Just (beg,"") = S.readInt beg'
        Just (end,"") = S.readInt end'
        str = if str' == "1" then id else reverseRange

    parse_tree s = case A.parseOnly (A.string "TREE " *> treep <* A.char ';') s of
        Left  e -> error $ "parse of TREE line failed with " ++ e ++ " on " ++ show s
        Right r -> r

    treep = Branch <$> (A.char '(' *> treep) <*> (A.char ',' *> treep) <*> (A.char ')' *> label) <|> Leaf <$> label
    label = (\x y z -> (x `S.append` S.pack ['[',y,']'], z))
            <$> A.takeWhile (\c -> A.isAlpha_ascii c || A.isDigit c || c == '_' || c == '.')
            <*> (A.char '[' *> (A.char '+' <|> A.char '-'))
            <*> (A.char ']' *> A.char ':' *> A.double)


data Range    = Range { r_seq    :: {-# UNPACK #-} !ShortByteString
                      , r_pos    :: {-# UNPACK #-} !Int
                      , r_length :: {-# UNPACK #-} !Int } deriving (Show, Eq, Ord)

reverseRange :: Range -> Range
reverseRange (Range sq pos len) = Range sq (-pos-len) len

data Label = Label Int Species Range deriving Show
data Tree label = Branch (Tree label) (Tree label) label | Leaf label deriving Show

relabel :: [Label] -> Tree (S.ByteString, Double) -> Tree Label
relabel ls (Leaf (lbl,_)) = Leaf (find_label lbl ls)
relabel ls (Branch u v (lbl,_)) = Branch (relabel ls u) (relabel ls v) (find_label lbl ls)

find_label :: S.ByteString -> [Label] -> Label
find_label lbl ls = case filter ((==) lbl . get_short) ls of
    [l] -> l
    []  -> error $ "not found: " ++ show lbl ++ " in " ++ show ls
    _   -> error $ "WTF?  " ++ show lbl ++ " in " ++ show ls

get_short :: Label -> S.ByteString
get_short (Label _ (spc,race) (Range chrom pos len)) = S.concat $
    [ S.map toUpper (S.take 1 spc),
      S.take 3 race, S.singleton '_',
      BS.fromShort chrom, S.singleton '_' ] ++
    if pos >= 0 then [ S.pack (show $ 1+pos), S.singleton '_',
                       S.pack (show $ pos+len), "[+]" ]
                else [ S.pack (show $ 1-pos-len), S.singleton '_',
                       S.pack (show $ -pos), "[-]" ]

data BlockDef = BlockDef { human_label :: !Label
                         , anc_label   :: !Label }

type Species = (S.ByteString,S.ByteString)

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
tree_has s (Leaf lbl@(Label _ spc _)) = if s == spc then [lbl] else []
tree_has s (Branch u' v' _) = tree_has s u' ++ tree_has s v'


-- | Returns a 'Fragment' for a subset of the species (say, the Chimps)
-- in a tree relative to a specific one (usually the Human).  If we find
-- multiple Human sequences, we return mutiple results.  If we find
-- multiple Chimps, we turn them into a consensus sequence and return
-- the result for the consensus sequence.
get_trees :: ( Tree Label -> [(Label, [Label])] ) -> EmfBlock -> [Fragment]
get_trees select  block = [ mkFragment rref (emf_seqs block !! iref) $ consensus_seq
                                [ emf_seqs block !! itgt | Label itgt _ _ <- tgt_lbls ]
                          | (Label iref _ rref, tgt_lbls) <- select $ emf_tree block ]


-- Fragmentary alignment reconstructed from EMF.  A list of these could
-- add up to a chromosome; we keep chromosomes neatly apart.
data Fragment = Fragment {-# UNPACK #-} !Range {-# UNPACK #-} !ShortByteString deriving Eq

frag_range :: Fragment -> Range
frag_range (Fragment r _) = r

frag_seqs :: Fragment -> (L.ByteString,L.ByteString)
frag_seqs (Fragment _ s) = frag_seqs1 s

frag_seqs1 :: ShortByteString -> (L.ByteString,L.ByteString)
frag_seqs1 s = L.splitAt (fromIntegral $ S.length s' `div` 2) (L.fromStrict s')
    where s' = Snappy.decompress $ BS.fromShort s

frag_raw :: Fragment -> ShortByteString
frag_raw (Fragment _ s) = s

instance Show Fragment where
    showsPrec _ frag = (++) "[ " . shows (r_seq rng) .  (:)  ':' . shows (r_pos rng) .
                       (++) ".." . shows (r_pos rng + r_length rng) .  (++) " ] " .
                       shows (L.take 25 rs) . (++) ".../" . shows (L.take 25 ss) . (++) "..."
      where
        rng     = frag_range frag
        (rs,ss) = frag_seqs frag


mkFragment :: Range -> S.ByteString -> S.ByteString -> Fragment
mkFragment rng rs ss
    | S.length rs == S.length ss &&
      S.length (S.filter (/='-') rs) == r_length rng
        = if r_pos rng >= 0
            then Fragment rng (compr2 rs ss)
            else Fragment (reverseRange rng) (compr2 (revcom rs) (revcom ss))
  where
    compr2 a b = case LB.unzip $ filter ((/=) (45,45)) $ SB.zip a b of
                    (u,v) -> BS.toShort $ Snappy.compress (L.toStrict $ u <> v)

    revcom = SB.reverse . SB.map complc

    complc 65 = 84
    complc 67 = 71
    complc 71 = 67
    complc 84 = 65

    complc  97 = 116
    complc  99 = 103
    complc 103 =  99
    complc 116 =  97

    complc c = c


consensus_seq :: [S.ByteString] -> S.ByteString
consensus_seq [   x  ] = x
consensus_seq (x:xs) | all (== S.length x) (map S.length xs) =
    S.pack [ foldl' (\r c -> if r == c then r else 'N') (S.index x i) (map (`S.index` i) xs)
           | i <- [0 .. S.length x -1] ]


-- | Scans a directory with gzipped EMF files in arbitrary order.  (The
-- files are not sorted internally, so there is no point in trying any
-- particular order.)
scan_all_emf :: FilePath -> IO [EmfBlock]
scan_all_emf dir = do
    fs <- map (dir </>) . filter (".emf.gz" `isSuffixOf`) <$> getDirectoryContents dir
    foldr (\fp k -> do bs  <- unsafeInterleaveIO k
                       let go [] = bs
                           go ls = case scanOneBlock ls of (b,ls') -> b : go ls'
                       go . map L.toStrict . L.lines . decompress <$> L.readFile fp)
          (return []) fs


apes :: [(S.ByteString, S.ByteString)]
apes = [("pan","troglodytes"),("gorilla","gorilla"),("pongo","abelii")
       ,("macaca","mulatta"),("callithrix","jacchus") ]


smoke_test :: IO ()
smoke_test = mapM_ print .
             concatMap (get_trees $ species_to_species ("homo","sapiens") ("pan","troglodytes")) =<<
             scan_all_emf "/mnt/expressions/martin/sequence_db/epo/epo_6_primate_v66"


data OptsEmf = OptsEmf {
    emf_output :: FilePath,
    emf_input  :: FilePath,
    emf_select :: Tree Label -> [(Label, [Label])] }

opts_emf_default :: OptsEmf
opts_emf_default = OptsEmf
    (error "no output specified")
    "/mnt/expressions/martin/sequence_db/epo/epo_6_primate_v66"
    (ancestor_to_species ("homo","sapiens") ("pan","troglodytes"))

opts_emf :: [ OptDescr ( OptsEmf -> IO OptsEmf ) ]
opts_emf =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "e" ["emf"] (ReqArg set_emf "PATH") "Compara alignments are at PATH"
    , Option "s" ["species"] (ReqArg set_species "NAME") "Import species NAME"
    , Option "a" ["ancestor"] (ReqArg set_anc "NAME") "Import ancestor with species NAME" ]
  where
    set_output  a c = return $ c { emf_output = a }
    set_emf     a c = return $ c { emf_input  = a }
    set_species a c = return $ c { emf_select = species_to_species  ("homo","sapiens") (w2 a) }
    set_anc     a c = return $ c { emf_select = ancestor_to_species ("homo","sapiens") (w2 a) }

    w2 s = case S.words (S.pack s) of (a:b:_) -> (a,b)


main_emf :: [String] -> IO ()
main_emf args = do
    ( _, OptsEmf{..} ) <- parseOpts False opts_emf_default (mk_opts "emf" "[emf-file...]" opts_emf) args
    !frags <- foldl' collect_frags empty_frag_map .
              concatMap (get_trees emf_select) <$>
              scan_all_emf emf_input

    hPutStrLn stderr "I survived!"
    hPrint stderr $ H.keys frags
    hPrint stderr $ map I.size $ H.elems frags

    withFile emf_output WriteMode $ \hdl ->
        hPutBuilder hdl . encode_hap $ catStretches
            [ fragmap1_to_stretch $ H.lookupDefault I.empty c frags
            | c <- map (BS.toShort . L.toStrict) chroms ]

    hPutStrLn stderr "All done!"



type FragMap  = H.HashMap ShortByteString FragMap1
type FragMap1 = I.IntMap ShortByteString

empty_frag_map :: FragMap
empty_frag_map = H.empty

collect_frags :: FragMap -> Fragment -> FragMap
collect_frags hm !frag =
    let m1 = H.lookupDefault I.empty (r_seq (frag_range frag)) hm
    in case I.lookup (r_pos (frag_range frag)) m1 of
            Just f' | f' == frag_raw frag -> hm
                    | otherwise  -> error $ "not-quite-duplicate: " ++ show (f',frag)
            Nothing -> H.insert (r_seq (frag_range frag))
                                (I.insert (r_pos (frag_range frag))
                                          (frag_raw frag) m1) hm

-- | Takes 'Fragment's and turns them into a 'Stretch'.  For good reasons,
-- this only works for one chromosome at a time.  (Recycles the logic
-- already developed for MAF input.)
fragmap1_to_stretch :: FragMap1 -> Stretch -> Stretch
fragmap1_to_stretch fm done = normMaf done $ go 0 $ I.toList fm
  where
    go  _  [        ] = MafDone
    go pos ((p,f):fs) = Skip (p - pos) $ Aln rs ss $ go (p+len) fs
      where
        (rs,ss) = frag_seqs1 f
        len     = fromIntegral $ L.length $ L.filter (/='-') rs


