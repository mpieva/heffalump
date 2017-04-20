{-# LANGUAGE BangPatterns, OverloadedStrings #-}
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
--                ridiculous).  Output of proper SNP files required the
--                reference.
--
--   VCF:  clearly requires the reference


import BasePrelude
import Foreign.Storable                     ( peekByteOff )
import System.IO.Posix.MMap                 ( unsafeMMapFile )
import System.IO.Unsafe                     ( unsafeDupablePerformIO )

import qualified Data.ByteString                as B
import qualified Data.ByteString.Char8          as C
import qualified Data.ByteString.Unsafe         as B
import qualified Data.Vector.Unboxed            as U

-- | This is a reference sequence.  It consists of stretches of Ns and
-- sequence.

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
                             , nrss_seqs    :: [ NewRefSeq ]
                             , nrss_path    :: !B.ByteString }

readTwoBit :: FilePath -> IO NewRefSeqs
readTwoBit fp = parseTwopBit fp <$> unsafeMMapFile fp

parseTwopBit :: FilePath -> B.ByteString -> NewRefSeqs
parseTwopBit fp0 raw = case (getW32 0, getW32 4) of (0x1A412743, 0) -> parseEachSeq 16 0
                                                    _ -> error "This does not look like a 2bit file."
  where
    getW32 :: Int -> Word32
    getW32 o | o + 4 >= B.length raw = error $ "out of bounds access: " ++ show (o, B.length raw)
             | otherwise = unsafeDupablePerformIO $ B.unsafeUseAsCString raw $ \p -> peekByteOff p o

    num_seq    = getW32 8

    parseEachSeq  _  n | num_seq == n = NewRefSeqs [] [] [] (fromString fp0)
    parseEachSeq off n = case parseEachSeq (off+5+nmsize) (n+1) of
                            NewRefSeqs cs ls ss fp ->
                                NewRefSeqs (name:cs) (dnasize:ls) (the_seq:ss) fp
      where
        !nmsize  = fromIntegral $ B.index raw off
        !name    = B.take nmsize $ B.drop (1+off) raw
        !offset  = fromIntegral . getW32 $ off+1+nmsize

        !dnasize = fromIntegral . getW32 $ offset
        !nBlockCount = fromIntegral . getW32 $ offset + 4
        !maskBlockCount = fromIntegral . getW32 $ offset + 8 + 8*nBlockCount
        !packedDnaOff = offset + 16 + 8 * (nBlockCount+maskBlockCount)

        -- We will need to decode the N blocks, we ignore the m blocks.
        n_blocks = [ ( fromIntegral . getW32 $ offset+8 + 4*i
                     , fromIntegral . getW32 $ offset+8 + 4*(i+nBlockCount) )
                   | i <- [ 0 .. nBlockCount-1 ] ]

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
    go c (r0:rs) = go1 r0 0
      where
        go1 _ _ [    ] = []
        go1 r p (v:vs)
            | c /= v_chr v = go (c+1) rs (v:vs)
            | p == v_pos v = case unconsNRS r of
                Just (c',_) -> v { v_ref = c' } : go1 r p vs
                Nothing     ->                    go1 r p vs
            | p < v_pos v  = go1 (dropNRS (v_pos v - p) r) (v_pos v) (v:vs)
            | otherwise    = error "expected sorted variants!"


