{-# LANGUAGE BangPatterns, OverloadedStrings, CPP #-}
module NewRef where

-- ^ We'll use a 2bit file as reference.  Chromosomes in the file will
-- be used in exactly this order, everything else is ignored.
--
-- When reading input, we always supply the reference file.  Sometimes we
-- don't actually need the reference sequence, but we still use it to
-- decide which chromomes to include, and sometimes how to reorder them.
-- (We should probably leave the list of terget chromosomes in the hef
-- files themselves; that way, we can ensure that two hefs cn be mnerged
-- safely.)
--
-- Applications for import:
--
--   VCF, BCF files:  A bunch of files, usually in the wrong order.  We
--                    read each, split them into chromosomes, then use
--                    the 2bit file to put them in the correct order.
--
--   BAM files:  must be ordered somehow, we require the reference, we
--               can reorder chromosomes after we're done
--
--   EMF, MAF:  will not be sorted.  We require the reference so we can
--              put everything into a proper order.
--
-- Applications for export:
--
--   D-stats, Treemix:  reference is no longer needed!
--
--   Eigenstrat:  can be done without the reference (but will look
--                ridiculous).  Output of proper SNP files
--                requires the reference.
--
--   VCF:  clearly requires the reference.


import Bio.Prelude                   hiding ( Ns )
import Data.ByteString.Short         hiding ( length, null )
import Foreign.Storable                     ( peekByteOff )
import System.Console.GetOpt
import System.Directory                     ( makeAbsolute )
import System.IO                            ( hPutStrLn, stdout, stderr )
import System.IO.Unsafe                     ( unsafeDupablePerformIO )

#if MIN_VERSION_biohazard(0,6,16)
import Bio.Util.MMap                        ( unsafeMMapFile )
#else
import System.IO.Posix.MMap                 ( unsafeMMapFile )
#endif

import qualified Data.ByteString                as B
import qualified Data.ByteString.Builder        as B
import qualified Data.ByteString.Char8          as C
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.ByteString.Short          as S
import qualified Data.ByteString.Unsafe         as B
import qualified Data.Vector.Unboxed            as U

import Util ( decomp, mk_opts, parseOpts )

-- | This is a reference sequence.  It consists of stretches of Ns and
-- sequence.  Invariant:  the lengths for 'ManyNs' and 'SomeSeq' are
-- always strictly greater than zero.

data NewRefSeq = ManyNs !Int NewRefSeq
                 -- Primitive bases in 2bit encoding:  [0..3] = TCAG
               | SomeSeq {-# UNPACK #-} !B.ByteString  -- ^ underlying data
                         {-# UNPACK #-} !Int           -- ^ offset in bases(!)
                         {-# UNPACK #-} !Int           -- ^ length in bases
                         NewRefSeq
               | NewRefEnd

instance Show NewRefSeq where
    showsPrec _ (ManyNs      n r) = (++) "ManyNs " . shows n . (++) " $ " . shows r
    showsPrec _ (SomeSeq _ _ l r) = (++) "SomeSeq " . shows l . (++) " $ " . shows r
    showsPrec _ NewRefEnd         = (++) "NewRefEnd"

data NewRefSeqs = NewRefSeqs { nrss_chroms  :: [ B.ByteString ]
                             , nrss_lengths :: [ Int ]
                             , nrss_seqs    :: [ () -> NewRefSeq ]
                             , nrss_path    :: !B.ByteString }

readTwoBit :: FilePath -> IO NewRefSeqs
readTwoBit fp = parseTwoBit <$> makeAbsolute fp <*> unsafeMMapFile fp

parseTwoBit :: FilePath -> B.ByteString -> NewRefSeqs
parseTwoBit fp0 raw = case (getW32 0, getW32 4) of (0x1A412743, 0) -> parseEachSeq 16 0
                                                   _ -> error "This does not look like a 2bit file."
  where
    getW32 :: Int -> Word32
    getW32 o | o + 4 >= B.length raw = error $ "out of bounds access: " ++ show (o, B.length raw)
             | otherwise = unsafeDupablePerformIO $ B.unsafeUseAsCString raw $ \p -> peekByteOff p o

    num_seq    = getW32 8

    parseEachSeq  _  n | num_seq == n = NewRefSeqs [] [] [] (fromString fp0)
    parseEachSeq off n = case parseEachSeq (off+5+nmsize) (n+1) of
                            NewRefSeqs cs ls ss fp ->
                                NewRefSeqs (name:cs) (dnasize:ls) ((\() -> the_seq):ss) fp
      where
        !nmsize  = fromIntegral $ B.index raw off
        !name    = B.take nmsize $ B.drop (1+off) raw
        !offset  = fromIntegral . getW32 $ off+1+nmsize

        !dnasize = fromIntegral . getW32 $ offset
        !nBlockCount = fromIntegral . getW32 $ offset + 4
        !maskBlockCount = fromIntegral . getW32 $ offset + 8 + 8*nBlockCount
        !packedDnaOff = offset + 16 + 8 * (nBlockCount+maskBlockCount)

        -- We will need to decode the N blocks, but we ignore the M blocks.
        n_blocks = [ (u,v) | i <- [ 0 .. nBlockCount-1 ]
                           , let u = fromIntegral . getW32 $ offset+8 + 4*i
                                 v = fromIntegral . getW32 $ offset+8 + 4*(i+nBlockCount) ]

        the_seq = unfoldSeq dnasize n_blocks 0 (packedDnaOff*4)


    unfoldSeq dnasize  _ chroff _fileoff | chroff >= dnasize = NewRefEnd
    unfoldSeq dnasize [] chroff  fileoff = SomeSeq raw fileoff (dnasize-chroff) NewRefEnd

    unfoldSeq dnasize ((nstart,nlen):ns) chroff fileoff

        | chroff  < nstart = SomeSeq raw fileoff (nstart-chroff)
                           $ unfoldSeq dnasize ((nstart,nlen):ns) nstart (fileoff+nstart-chroff)

        | chroff == nstart = ManyNs nlen
                           $ unfoldSeq dnasize ns (chroff+nlen) (fileoff+nlen)

        | otherwise = error "WTF?!"



{-# INLINE lengthNRS #-}
lengthNRS :: NewRefSeq -> Int64
lengthNRS = go 0
  where
    go !acc (ManyNs      l s) = go (acc + fromIntegral l) s
    go !acc (SomeSeq _ _ l s) = go (acc + fromIntegral l) s
    go !acc NewRefEnd         = acc


{-# INLINE nullNRS #-}
nullNRS :: NewRefSeq -> Bool
nullNRS = go
  where
    go (ManyNs      0 s) = go s
    go (SomeSeq _ _ 0 s) = go s
    go  NewRefEnd        = True
    go  _                = False


{-# INLINE dropNRS #-}
dropNRS :: Int -> NewRefSeq -> NewRefSeq
dropNRS 0         s = s
dropNRS _ NewRefEnd = NewRefEnd
dropNRS n (ManyNs l s) | n > l     = dropNRS (n-l) s
                       | n < l     = ManyNs (l-n) s
                       | otherwise = s
dropNRS n (SomeSeq r o l s) | n > l     = dropNRS (n-l) s
                            | n < l     = SomeSeq r (o+n) (l-n) s
                            | otherwise = s

{-# INLINE takeNRS #-}
takeNRS :: Int -> NewRefSeq -> NewRefSeq
takeNRS 0         _ = NewRefEnd
takeNRS _ NewRefEnd = NewRefEnd
takeNRS n (ManyNs l s) | n >= l    = ManyNs l $ takeNRS (n-l) s
                       | otherwise = ManyNs n NewRefEnd
takeNRS n (SomeSeq r o l s) | n <= l    = SomeSeq r o n NewRefEnd
                            | otherwise = SomeSeq r o l $ takeNRS (n-l) s

unpackNRS :: NewRefSeq -> L.ByteString
unpackNRS = L.unfoldr $ fmap (first toCode) . unconsNRS
  where
    toCode (N2b 0) = 'T'
    toCode (N2b 1) = 'C'
    toCode (N2b 2) = 'A'
    toCode (N2b 3) = 'G'
    toCode      _  = 'N'

-- | Nucleotide in 2bit encoding: "TCAG" == [0..3], N == 255.
newtype Nuc2b = N2b Word8 deriving Eq

-- | Variant in 2bit encoding: [0..3] for "IOPX"
newtype Var2b = V2b Word8 deriving Eq

isTransversion :: Var2b -> Bool
isTransversion (V2b v) = testBit v 1

toAltCode :: Var2b -> Nuc2b -> Char
toAltCode (V2b v) (N2b r) = C.index "TCAGXPOI" $ fromIntegral (xor r v .&. 7)

toRefCode :: Nuc2b -> Char
toRefCode = toAltCode (V2b 0)


instance Show Nuc2b where
    show (N2b 0) = "T"
    show (N2b 1) = "C"
    show (N2b 2) = "A"
    show (N2b 3) = "G"
    show      _  = "N"

instance Show Var2b where
    show (V2b 0) = "I"
    show (V2b 1) = "O"
    show (V2b 2) = "P"
    show (V2b 3) = "X"
    show      _  = "?"

{-# INLINE unconsNRS #-}
unconsNRS :: NewRefSeq -> Maybe ( Nuc2b, NewRefSeq )
unconsNRS NewRefEnd = Nothing
unconsNRS (ManyNs 1 s) = Just (N2b 255, s)
unconsNRS (ManyNs l s) = Just (N2b 255, ManyNs (l-1) s)
unconsNRS (SomeSeq raw off len s) = Just (c, s')
  where
    c  = N2b . fromIntegral $ (B.index raw (off `shiftR` 2) `shiftR` (6 - 2 * (off .&. 3))) .&. 3
    s' = if len == 1 then s else SomeSeq raw (off+1) (len-1) s

data RefSeqView a = !Int :== a | !Nuc2b :^ a | NilRef

viewNRS :: NewRefSeq -> RefSeqView NewRefSeq
viewNRS  NewRefEnd              = NilRef
viewNRS (ManyNs            l s) = l :== s
viewNRS (SomeSeq raw off len s) = c :^  s'
  where
    c  = N2b . fromIntegral $ (B.index raw (off `shiftR` 2) `shiftR` (6 - 2 * (off .&. 3))) .&. 3
    s' = if len == 1 then s else SomeSeq raw (off+1) (len-1) s

-- A variant call.  Has a genomic position, ref and alt alleles, and a
-- bunch of calls.
data Variant = Variant { v_chr   :: !Int                -- chromosome number
                       , v_pos   :: !Int                -- 0-based
                       , v_ref   :: !Nuc2b
                       , v_alt   :: !Var2b
                       , v_calls :: !(U.Vector Word8) } -- Variant codes:  #ref + 4 * #alt
  deriving Show

addRef :: NewRefSeqs -> [Variant] -> [Variant]
addRef ref = go 0 (nrss_seqs ref)
  where
    go _ [     ] = const []
    go c (r0:rs) = go1 (r0 ()) 0
      where
        go1 _ _ [    ] = []
        go1 r p (v:vs)
            | c /= v_chr v = go (c+1) rs (v:vs)
            | p == v_pos v = case unconsNRS r of
                Just (c',_) -> v { v_ref = c' } : go1 r p vs
                Nothing     ->                    go1 r p vs
            | p < v_pos v  = go1 (dropNRS (v_pos v - p) r) (v_pos v) (v:vs)
            | otherwise    = error "expected sorted variants!"


has_help :: [String] -> Bool
has_help args = null args || any (`elem` h_opts) (takeWhile (/= "--") args)
  where
    h_opts = ["-h","--help","-?","--usage"]

eff_args :: [String] -> [String]
eff_args args = case break (== "--") args of (l,r) -> l ++ drop 1 r

main_2bitinfo :: [String] -> IO ()
main_2bitinfo fs
    | has_help fs = do
        pn <- getProgName
        hPutStrLn stderr $ "Usage: " ++ pn ++ " twobitinfo [2bit-file...]\n"
        exitFailure

    | otherwise = forM_ (eff_args fs) $ \f -> do
        ref <- readTwoBit f
        sequence_ $ zipWith (printf "%s\t%d\n" . C.unpack)
                            (nrss_chroms ref) (nrss_lengths ref)

main_2bittofa :: [String] -> IO ()
main_2bittofa fs = case eff_args fs of
    _ | has_help fs -> do
            pn <- getProgName
            hPutStrLn stderr $ "Usage: " ++ pn ++
                " twobittofa [2bit-file] [[name|name:start-end]...]\n"
            exitFailure

    [] -> return () -- can't happen

    [fp] -> do ref <- readTwoBit fp
               forM_ (nrss_chroms ref `zip` nrss_seqs ref) $ \(ch,sq) ->
                    putStrLn ('>' : C.unpack ch) >> twoBitToFa (sq ())

    (fp:rns) -> do ref <- readTwoBit fp
                   forM_ rns $ \rn ->
                        case findIndex ((==) rn . C.unpack) (nrss_chroms ref) of
                            Just i  -> do putStrLn $ '>' : rn
                                          twoBitToFa $ (nrss_seqs ref !! i) ()
                            Nothing -> do let (ch,':':s1) = break ((==) ':') rn
                                          [(start,'-':s2)] <- return $ reads s1
                                          [(end,"")] <- return $ reads s2
                                          Just i <- return $ findIndex ((==) ch . C.unpack) (nrss_chroms ref)

                                          printf ">%s:%d-%d\n" ch start end
                                          twoBitToFa $ dropNRS start $ takeNRS end $ (nrss_seqs ref !! i) ()

twoBitToFa :: NewRefSeq -> IO ()
twoBitToFa = splitLns . unpackNRS
  where
    splitLns s
        | L.null s = return ()
        | otherwise = case L.splitAt 50 s of { (l,r) -> do
                            L.hPut stdout l
                            L.hPut stdout "\n"
                            splitLns r }

opts_fato2bit :: [ OptDescr (FilePath -> IO FilePath) ]
opts_fato2bit = [ Option "o" ["output"] (ReqArg (\a _ -> return a) "FILE") "Write output to FILE" ]

main_fato2bit :: [String] -> IO ()
main_fato2bit args = do
    ( fs, fp ) <- parseOpts True (error "no output file given")
                            (mk_opts "fatotwobit" "[fasta-file...]" opts_fato2bit) args

    L.writeFile fp . faToTwoBit . L.concat . map decomp =<<
        mapM l'readFile (if null fs then ["-"] else fs)
  where
    l'readFile "-" = L.getContents
    l'readFile  f  = L.readFile f


-- List of pairs of 'Word32's.  Specialized and unpacked to conserve
-- space.  Probably overkill...
data L2i = L2i {-# UNPACK #-} !Word32 {-# UNPACK #-} !Word32 L2i | L2i_Nil

encodeL2i :: L2i -> B.Builder
encodeL2i = go 0 mempty mempty
  where
    go !n ss ls  L2i_Nil     = B.word32LE n <> ss <> ls
    go !n ss ls (L2i s e rs) = go (succ n) (B.word32LE s <> ss) (B.word32LE (e-s) <> ls) rs


-- Strategy:  We can oly write the packedDNA after we wrote the nBlocks
-- and mBlocks.  So packedDNA needs to be buffered.  We have to do three
-- simultaneous strict folds of the input, all of which result in reasonably
-- compact structures, which get concatenated at the end.
--
-- We also have to buffer everything, since the header with the sequence
-- names must be written first.  Oh joy.

faToTwoBit :: L.ByteString -> L.ByteString
faToTwoBit s0 = L.concat $ B.toLazyByteString header : map snd seqs
  where
    seqs = get_each s0

    offset0 = 16 + 5 * length seqs + sum (map (S.length . fst) seqs)
    offsets = scanl (\a b -> a + fromIntegral (L.length b)) offset0 $ map snd seqs

    header = B.word32LE 0x1A412743 <> B.word32LE 0 <>
             B.word32LE (fromIntegral $ length seqs) <> B.word32LE 0 <>
             fold (zipWith (\nm off -> B.word8 (fromIntegral (S.length nm)) <>
                                       B.shortByteString nm <>
                                       B.word32LE (fromIntegral off))
                           (map fst seqs) offsets)

    get_each s = let s1 = L.dropWhile (/= '>') s
                     nm = L.takeWhile (not . isSpace) $ L.drop 1 s1
                     s2 = L.dropWhile (/= '\n') s1
                 in if L.null s1 then []
                    else get_one (toShort $ L.toStrict nm) 0
                                 (maxBound :!: L2i_Nil) (maxBound :!: L2i_Nil) (0 :!: 0 :!: []) s2


    get_one !nm !pos !_ns !_ms !_bs _s
        | pos .&. 0xFFFFFF == 0 && trace (show (nm,pos)) False = undefined

    get_one !nm !pos !ns !ms !bs s = case L.uncons s of
        Nothing -> fin
        Just (c,s')
            | isSpace c -> get_one nm pos ns ms bs s'
            | c == '>'  -> fin
            | otherwise -> get_one nm (succ pos)
                                   (collect_Ns ns pos c)
                                   (collect_ms ms pos c)
                                   (collect_bases bs c) s'
      where
        fin = let ss' = case bs of (0 :!: _ :!: ss) -> ss
                                   (n :!: w :!: ss) -> B.singleton (w `shiftL` (6-2*n)) : ss
                  raw = B.toLazyByteString $
                            B.word32LE pos <>
                            encodeL2i (case ns of p :!: rs | p == maxBound -> rs ; p :!: rs -> L2i p pos rs) <>
                            encodeL2i (case ms of p :!: rs | p == maxBound -> rs ; p :!: rs -> L2i p pos rs) <>
                            B.word32LE 0 <>
                            foldMap B.byteString (reverse ss')
              in L.length raw `seq` (nm, raw) : get_each s

    collect_Ns (spos :!: rs) pos c
        | spos == maxBound && c `C.elem` "ACGTacgt" = maxBound :!: rs
        | spos == maxBound                          =      pos :!: rs
        |                     c `C.elem` "ACGTacgt" = maxBound :!: L2i spos pos rs
        | otherwise                                 =     spos :!: rs

    collect_ms (spos :!: rs) pos c
        | spos == maxBound && isUpper c = maxBound :!: rs
        | spos == maxBound              =      pos :!: rs
        |                     isUpper c = maxBound :!: L2i spos pos rs
        | otherwise                     =     spos :!: rs

    -- collect 4 bases in w, then collect bytes in a list of byte
    -- strings of increasing length
    -- packedDna - the DNA packed to two bits per base, represented as
    --             so: T - 00, C - 01, A - 10, G - 11. The first base is
    --             in the most significant 2-bit byte; the last base is
    --             in the least significant 2 bits. For example, the
    --             sequence TCAG is represented as 00011011.
    collect_bases (n :!: w :!: ss) c
        = let code = case c of 'C'->1;'c'->1;'A'->2;'a'->2;'G'->3;'g'->3;_->0
              w' = shiftL w 2 .|. code
          in if n == 3 then 0 :!: 0 :!: put w' ss else succ n :!: w' :!: ss

    -- Keep appending bytes to a collection of 'ByteString's in such a
    -- way that the 'ByteString's keep doubling in size.  (XXX  This is
    -- a recurring problem; could use a reusable solution.)
    put w = go 1 [B.singleton w]
      where
        go l acc (s:ss)
            | B.length s <= l = go (l+B.length s) (s:acc) ss
            | otherwise = let !s' = B.concat acc in s' : s : ss
        go _ acc [    ] = let !s' = B.concat acc in [s']

