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

import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Unsafe as B
import System.IO.Posix.MMap
import System.IO.Unsafe ( unsafeDupablePerformIO )
import BasePrelude
import Foreign.Storable ( peekByteOff )


-- | This is a reference sequence.  It consists of stretches of Ns and
-- sequence.

data NewRefSeq = ManyNs !Int NewRefSeq | SomeSeq !PrimBases NewRefSeq | NewRefEnd
    deriving Show

-- | Primitive bases in 2bit encoding:  [0..3] = TCAG
data PrimBases = PrimBases {-# UNPACK #-} !B.ByteString  -- ^ underlying data
                           {-# UNPACK #-} !Int           -- ^ offset in bases(!)
                           {-# UNPACK #-} !Int           -- ^ length in bases

instance Show PrimBases where
    show (PrimBases _ _ l) = show l

type NewRefSeqs = [( B.ByteString, Int, NewRefSeq )]

readTwoBit :: FilePath -> IO NewRefSeqs
readTwoBit fp = parseTwopBit <$> unsafeMMapFile fp

parseTwopBit :: B.ByteString -> NewRefSeqs
parseTwopBit raw = case (getW32 0, getW32 4) of (0x1A412743, 0) -> parseEachSeq 16 0
  where
    getW32 :: Int -> Word32
    getW32 o | o + 4 >= B.length raw = error "out of bounds access"
             | otherwise = unsafeDupablePerformIO $ B.unsafeUseAsCString raw $ \p -> peekByteOff p o

    num_seq    = getW32 8

    parseEachSeq off n
        | num_seq == n = []
        | otherwise    = ( name, dnasize, the_seq ) : parseEachSeq (off+5+nmsize) (n+1)
      where
        nmsize  = fromIntegral $ B.index raw off
        name    = B.take nmsize $ B.drop (1+off) raw
        offset  = fromIntegral . getW32 $ off+1+nmsize

        dnasize = fromIntegral . getW32 $ offset
        nBlockCount = fromIntegral . getW32 $ offset + 4
        maskBlockCount = fromIntegral . getW32 $ offset + 8 + 8*nBlockCount
        packedDnaOff = offset + 16 + 8 * (nBlockCount+maskBlockCount)

        -- We will need to decode the N blocks, we ignore the m blocks.
        n_blocks = [ ( fromIntegral . getW32 $ offset+8 + 4*i
                     , fromIntegral . getW32 $ offset+8 + 4*(i+nBlockCount) )
                   | i <- [ 0 .. nBlockCount-1 ] ]

        the_seq = unfoldSeq dnasize n_blocks 0 (packedDnaOff*4)


    unfoldSeq dnasize  _ chroff fileoff | chroff >= dnasize = NewRefEnd
    unfoldSeq dnasize [] chroff fileoff = SomeSeq (PrimBases raw fileoff (dnasize-chroff)) NewRefEnd

    unfoldSeq dnasize ((nstart,nlen):ns) chroff fileoff

        | chroff < nstart = SomeSeq (PrimBases raw fileoff (nstart-chroff))
                          $ unfoldSeq dnasize ((nstart,nlen):ns) nstart (fileoff+nstart-chroff)

        | chroff == nstart = ManyNs nlen
                          $ unfoldSeq dnasize ns (chroff+nlen) (fileoff+nlen)

        | otherwise = error "WTF?!"



{-# INLINE lengthNRS #-}
lengthNRS :: NewRefSeq -> Int64
lengthNRS = go 0
  where
    go !acc (ManyNs l s)                  = go (acc + fromIntegral l) s
    go !acc (SomeSeq (PrimBases _ _ l) s) = go (acc + fromIntegral l) s
    go !acc NewRefEnd                     = acc


{-# INLINE nullNRS #-}
nullNRS :: NewRefSeq -> Bool
nullNRS = go
  where
    go (ManyNs 0 s)                  = go s
    go (SomeSeq (PrimBases _ _ 0) s) = go s
    go  NewRefEnd                    = True
    go  _                            = False


{-# INLINE dropNRS #-}
dropNRS :: Int -> NewRefSeq -> NewRefSeq
dropNRS 0         s = s
dropNRS _ NewRefEnd = NewRefEnd
dropNRS n (ManyNs l s) | n > l     = dropNRS (n-l) s
                       | n < l     = ManyNs (l-n) s
                       | otherwise = s
dropNRS n (SomeSeq (PrimBases r o l) s) | n+o > l   = dropNRS (n-l) s
                                        | n+o < l   = SomeSeq (PrimBases r (o+n) l) s
                                        | otherwise = s


{-# INLINE unconsNRS #-}
unconsNRS :: NewRefSeq -> Maybe ( Char, NewRefSeq )
unconsNRS NewRefEnd = Nothing
unconsNRS (ManyNs 1 s) = Just ('N', s)
unconsNRS (ManyNs l s) = Just ('N', ManyNs (l-1) s)
unconsNRS (SomeSeq (PrimBases raw off len) s) = Just (c, s')
  where
    c = "TCAG" `S.index` fromIntegral ((B.index raw (off `shiftR` 2) `shiftR` (6 - 2 * (off .&. 3))) .&. 3)
    s' = if off+1 == len then s else SomeSeq (PrimBases raw (off+1) len) s


