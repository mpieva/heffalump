{-# LANGUAGE CPP #-}
module Genome where

-- ^ We'll use a 2bit file as reference.  Chromosomes in the file will
-- be used in exactly this order, everything else is ignored.
--
-- When reading input, we always supply the reference file.  Sometimes we
-- don't actually need the reference sequence, but we still use it to
-- decide which chromomes to include, and sometimes how to reorder them.
-- We also compute a checksum over the list of target chromosomes and
-- leave it in the hef files themselves; that way, we can ensure that
-- two hefs can be merged safely.
--
-- Applications for import:
--
--   VCF, BCF files:  A bunch of files, usually in the wrong order.  We
--                    read each, split them into chromosomes, then use
--                    the 2bit file to put them in the correct order.
--
--   BAM files:  Must be sorted, but the order of target sequences is
--               arbitrary.  We require the reference, and reorder valid
--               chromosomes after we're done importing.
--
--   EMF, MAF:  Will not be sorted.  We parse pieces in arbitrary order,
--              and require the reference so we can put everything into
--              a proper order at the end.
--
-- Applications for export:
--
--   D-stats, Treemix:  The reference is no longer needed!
--
--   Eigenstrat:  Can be done without the reference, but will look
--                ridiculous.  Output of proper SNP files requires the
--                reference.
--
--   VCF:  Clearly requires the reference.


import Bio.Prelude                   hiding ( Ns )
import Data.ByteString.Short         hiding ( length, null )
import Foreign.Storable                     ( peekByteOff )
import Streaming
import System.Console.GetOpt
import System.Directory                     ( makeAbsolute )
import System.IO                            ( hPutStrLn, stdout, stderr, withFile, IOMode(..) )
import System.IO.Unsafe                     ( unsafeDupablePerformIO )

#if MIN_VERSION_biohazard(0,6,16)
import Bio.Util.MMap                        ( unsafeMMapFile )
#else
import System.IO.Posix.MMap                 ( unsafeMMapFile )
#endif

import qualified Data.ByteString                    as B
import qualified Data.ByteString.Builder            as B
import qualified Data.ByteString.Char8              as C
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.ByteString.Short              as H
import qualified Data.ByteString.Streaming.Char8    as S
import qualified Data.ByteString.Unsafe             as B
import qualified Data.Vector.Unboxed                as U
import qualified Streaming.Prelude                  as Q

import Util ( decomp, mk_opts, parseFileOpts, unexpected )

-- | This is a reference sequence.  It consists of stretches of Ns and
-- stretches of sequence.  Invariant:  the lengths for 'ManyNs' and
-- 'SomeSeq' are always strictly greater than zero.

data RefSeq = ManyNs !Int RefSeq
              -- Primitive bases in 2bit encoding:  [0..3] = TCAG
            | SomeSeq {-# UNPACK #-} !B.ByteString  -- ^ underlying data
                      {-# UNPACK #-} !Int           -- ^ offset in bases(!)
                      {-# UNPACK #-} !Int           -- ^ length in bases
                      RefSeq
            | RefEnd

instance Show RefSeq where
    showsPrec _ (ManyNs      n r) = (++) "ManyNs " . shows n . (++) " $ " . shows r
    showsPrec _ (SomeSeq _ _ l r) = (++) "SomeSeq " . shows l . (++) " $ " . shows r
    showsPrec _  RefEnd           = (++) "RefEnd"

data RefSeqs = RefSeqs { rss_chroms  :: [ B.ByteString ]
                       , rss_lengths :: [ Int ]
                       , rss_seqs    :: [ () -> RefSeq ]
                       , rss_path    :: !B.ByteString }

readTwoBit :: FilePath -> IO RefSeqs
readTwoBit fp = parseTwoBit <$> makeAbsolute fp <*> unsafeMMapFile fp

-- In theory, there could be 2bit files in big endian format out there.
-- We don't support them, since I've never seen one in the wild.
parseTwoBit :: FilePath -> B.ByteString -> RefSeqs
parseTwoBit fp0 raw = case (getW32 0, getW32 4) of (0x1A412743, 0) -> parseEachSeq 16 0
                                                   _ -> error "This does not look like a 2bit file."
  where
    getW32 :: Int -> Word32
    getW32 o | o + 4 >= B.length raw = error $ "out of bounds access: " ++ show (o, B.length raw)
             | otherwise = unsafeDupablePerformIO $ B.unsafeUseAsCString raw $ \p -> peekByteOff p o

    num_seq  = getW32 8

    parseEachSeq  _  n | num_seq == n = RefSeqs [] [] [] (fromString fp0)
    parseEachSeq off n = case parseEachSeq (off+5+nmsize) (n+1) of
                            RefSeqs cs ls ss fp ->
                                RefSeqs (name:cs) (dnasize:ls) ((\() -> the_seq):ss) fp
      where
        !nmsize  = fromIntegral $ B.index raw off
        !name    = B.take nmsize $ B.drop (1+off) raw
        !offset  = fromIntegral . getW32 $ off+1+nmsize

        !dnasize        = fromIntegral . getW32 $ offset
        !nBlockCount    = fromIntegral . getW32 $ offset + 4
        !maskBlockCount = fromIntegral . getW32 $ offset + 8 + 8*nBlockCount
        !packedDnaOff   = offset + 16 + 8 * (nBlockCount+maskBlockCount)

        -- We will need to decode the N blocks, but we ignore the M blocks.
        n_blocks = [ (u,v) | i <- [ 0 .. nBlockCount-1 ]
                           , let u = fromIntegral . getW32 $ offset+8 + 4*i
                           , let v = fromIntegral . getW32 $ offset+8 + 4*(i+nBlockCount) ]

        the_seq = unfoldSeq dnasize n_blocks 0 (packedDnaOff*4)


    unfoldSeq dnasize  _ chroff _fileoff | chroff >= dnasize = RefEnd
    unfoldSeq dnasize [] chroff  fileoff = SomeSeq raw fileoff (dnasize-chroff) RefEnd

    unfoldSeq dnasize ((nstart,nlen):ns) chroff fileoff

        | chroff  < nstart = SomeSeq raw fileoff (nstart-chroff)
                           $ unfoldSeq dnasize ((nstart,nlen):ns) nstart (fileoff+nstart-chroff)

        | chroff == nstart = ManyNs nlen
                           $ unfoldSeq dnasize ns (chroff+nlen) (fileoff+nlen)

        | otherwise = unexpected ""


{-# INLINE lengthRS #-}
lengthRS :: RefSeq -> Int64
lengthRS = go 0
  where
    go !acc (ManyNs      l s) = go (acc + fromIntegral l) s
    go !acc (SomeSeq _ _ l s) = go (acc + fromIntegral l) s
    go !acc  RefEnd           = acc

{-# INLINE nullRS #-}
nullRS :: RefSeq -> Bool
nullRS = go
  where
    go (ManyNs      0 s) = go s
    go (SomeSeq _ _ 0 s) = go s
    go  RefEnd           = True
    go  _                = False

{-# INLINE dropRS #-}
dropRS :: Int -> RefSeq -> RefSeq
dropRS 0                s              = s
dropRS _  RefEnd                       = RefEnd
dropRS n (ManyNs      l s) | n > l     = dropRS (n-l) s
                           | n < l     = ManyNs (l-n) s
                           | otherwise = s
dropRS n (SomeSeq r o l s) | n > l     = dropRS (n-l) s
                           | n < l     = SomeSeq r (o+n) (l-n) s
                           | otherwise = s

{-# INLINE takeRS #-}
takeRS :: Int -> RefSeq -> RefSeq
takeRS 0                _              = RefEnd
takeRS _  RefEnd                       = RefEnd
takeRS n (ManyNs      l s) | n >= l    = ManyNs l $ takeRS (n-l) s
                           | otherwise = ManyNs n RefEnd
takeRS n (SomeSeq r o l s) | n <= l    = SomeSeq r o n RefEnd
                           | otherwise = SomeSeq r o l $ takeRS (n-l) s

unpackRS :: RefSeq -> L.ByteString
unpackRS = L.unfoldr $ fmap (first toCode) . unconsRS
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
toAltCode (V2b 255)    _  = '.'
toAltCode (V2b v) (N2b r) = C.index "TCAGXPOI" $ fromIntegral (xor r v .&. 7)

toRefCode :: Nuc2b -> Char
toRefCode (N2b r) = C.index "NTCAG" . fromIntegral $ r+1


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

{-# INLINE unconsRS #-}
unconsRS :: RefSeq -> Maybe ( Nuc2b, RefSeq )
unconsRS  RefEnd                 = Nothing
unconsRS (ManyNs            1 s) = Just (N2b 255, s)
unconsRS (ManyNs            l s) = Just (N2b 255, ManyNs (l-1) s)
unconsRS (SomeSeq raw off len s) = Just (c, s')
  where
    c  = N2b . fromIntegral $ (B.index raw (off `shiftR` 2) `shiftR` (6 - 2 * (off .&. 3))) .&. 3
    s' = if len == 1 then s else SomeSeq raw (off+1) (len-1) s

data RefSeqView a = !Int :== a | !Nuc2b :^ a | NilRef

viewRS :: RefSeq -> RefSeqView RefSeq
viewRS  RefEnd                 = NilRef
viewRS (ManyNs            l s) = l :== s
viewRS (SomeSeq raw off len s) = c :^  s'
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

addRef :: Monad m => RefSeqs -> Stream (Of Variant) m r -> Stream (Of Variant) m r
addRef ref = go 0 (rss_seqs ref)
  where
    go _ [     ] = lift . Q.effects
    go c (r0:rs) = go1 (r0 ()) 0
      where
        go1 r p = lift . Q.next >=> \case
            Left x              -> pure x
            Right (v,vs)
                | c /= v_chr v  -> go (c+1) rs (Q.cons v vs)
                | p == v_pos v -> case unconsRS r of
                    Just (c',_) -> v { v_ref = c' } `Q.cons` go1 r p vs
                    Nothing     ->                           go1 r p vs
                | p < v_pos v   -> go1 (dropRS (v_pos v - p) r) (v_pos v) (Q.cons v vs)
                | otherwise     -> error "expected sorted variants!"


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
        zipWithM_ (printf "%s\t%d\n" . C.unpack)
                  (rss_chroms ref) (rss_lengths ref)

main_2bittofa :: [String] -> IO ()
main_2bittofa fs = case eff_args fs of
    _ | has_help fs -> do
            pn <- getProgName
            hPutStrLn stderr $ "Usage: " ++ pn ++
                " twobittofa [2bit-file] [[name|name:start-end]...]\n"
            exitFailure

    [] -> return () -- can't happen

    [fp] -> do ref <- readTwoBit fp
               forM_ (rss_chroms ref `zip` rss_seqs ref) $ \(ch,sq) ->
                    putStrLn ('>' : C.unpack ch) >> twoBitToFa (sq ())

    (fp:rns) -> do ref <- readTwoBit fp
                   forM_ rns $ \rn ->
                        case elemIndex (C.pack rn) (rss_chroms ref) of
                            Just i  -> do putStrLn $ '>' : rn
                                          twoBitToFa $ (rss_seqs ref !! i) ()
                            Nothing -> sequence_
                                          [ do printf ">%s:%d-%d\n" ch start end
                                               twoBitToFa $ dropRS start $ takeRS end $ (rss_seqs ref !! i) ()
                                          | (ch,':':s1)    <- [break (':' ==) rn]
                                          , (start,'-':s2) <- reads s1
                                          , (end,"")       <- reads s2
                                          , i              <- elemIndices (C.pack ch) (rss_chroms ref) ]

twoBitToFa :: RefSeq -> IO ()
twoBitToFa = splitLns . unpackRS
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
    ( fs, fp ) <- parseFileOpts (error "no output file given")
                                (mk_opts "fatotwobit" "[fasta-file...]" opts_fato2bit) args

    withFiles (if null fs then ["-"] else fs) $
        L.writeFile fp <=< faToTwoBit
  where
    withFiles [      ] k = k $ pure ()
    withFiles ("-":fs) k = withFiles fs $ \s ->
                               k $ decomp (S.fromHandle stdin) >> s
    withFiles ( f :fs) k = withFile f ReadMode $ \h ->
                             withFiles fs $ \s ->
                               k $ decomp (S.fromHandle h) >> s


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

faToTwoBit :: Monad m => S.ByteString m r -> m L.ByteString
faToTwoBit s0 = do
    seqs <- get_each [] s0

    let offset0 = 16 + 5 * length seqs + sum (map (H.length . fst) seqs)
        offsets = scanl (\a b -> a + fromIntegral (L.length b)) offset0 $ map snd seqs

        header = B.word32LE 0x1A412743 <> B.word32LE 0 <>
                 B.word32LE (fromIntegral $ length seqs) <> B.word32LE 0 <>
                 fold (zipWith (\nm off -> B.word8 (fromIntegral (H.length nm)) <>
                                           B.shortByteString nm <>
                                           B.word32LE (fromIntegral off))
                               (map fst seqs) offsets)

    return $ L.concat $ B.toLazyByteString header : map snd seqs

  where
    get_each acc s = do s1 <- S.uncons $ S.dropWhile (/= '>') s
                        case s1 of
                            Left    _    -> return $ reverse acc
                            Right (_,s2) -> do
                                nm :> s' <- S.toStrict $ S.break isSpace s2
                                get_one acc (toShort nm) 0 (maxBound :!: L2i_Nil)
                                        (maxBound :!: L2i_Nil) (0 :!: 0 :!: [])
                                        (S.dropWhile (/= '\n') s')

    get_one acc !nm !pos !ns !ms !bs = S.uncons >=> \case
        Left    r       -> fin (pure r)
        Right (c,s')
            | isSpace c -> get_one acc nm pos ns ms bs s'
            | c == '>'  -> fin (S.cons c s')
            | otherwise -> get_one acc nm (succ pos)
                                   (collect_Ns ns pos c)
                                   (collect_ms ms pos c)
                                   (collect_bases bs c) s'
      where
        fin s = let ss' = case bs of (0 :!: _ :!: ss) -> ss
                                     (n :!: w :!: ss) -> B.singleton (w `shiftL` (6-2*n)) : ss
                    raw = B.toLazyByteString $
                              B.word32LE pos <>
                              encodeL2i (case ns of p :!: rs | p == maxBound -> rs ; p :!: rs -> L2i p pos rs) <>
                              encodeL2i (case ms of p :!: rs | p == maxBound -> rs ; p :!: rs -> L2i p pos rs) <>
                              B.word32LE 0 <>
                              foldMap B.byteString (reverse ss')
                in L.length raw `seq` get_each ((nm, raw) : acc) s

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
              w'   = shiftL w 2 .|. code
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

