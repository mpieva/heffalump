module VcfScan
    ( withScanner
    , getVariant
    , RawVariant(..)
    ) where

import Bio.Prelude
import Data.ByteString.Short.Internal   ( ShortByteString, createFromPtr, toShort )
import Foreign.C.String
import Foreign.C.Types
import Foreign.Marshal.Alloc            ( allocaBytes )
import Foreign.Marshal.Utils            ( moveBytes, copyBytes )

import qualified Data.ByteString.Char8      as B
import qualified Data.ByteString.Unsafe     as B ( unsafeUseAsCString )
import qualified Data.ByteString.Streaming  as S
import qualified Data.HashMap.Strict        as M

#include "zlib.h"
#include "vcf-scan.h"

-- | Supports one individual only.  (Project your VCF if you have more
-- than one; better yet, use BCF.)
data RawVariant = RawVariant {
    rv_chrom :: !Int,                                -- chromosome index (0-based)
    rv_pos   :: !Int,                                -- position (1-based)
    rv_vars  :: !Bytes,                              -- vars, comma separated, incl. ref
    rv_gt    :: !Word16 }                            -- genotype
        deriving (Show, Eq)

data CScanner

workbuffer_size :: Int
workbuffer_size = 127*1024

-- | Initialize scanner.  Allocates an input and a working buffer, fills
-- the working buffer from the input, scans the header and sets up the
-- c structure.
withScanner :: [Bytes] -> S.ByteString IO r
            -> (Ptr CScanner -> HashMap ShortByteString Int -> S.ByteString IO r -> IO b)
            -> IO b
withScanner chroms str k =
    allocaBytes (#{const sizeof(scanner)})          $ \psc ->
    allocaBytes workbuffer_size                     $ \work_buf -> do

    (#{ poke scanner, working_buffer } psc work_buf)
    (#{ poke scanner, next_work } psc work_buf)
    (#{ poke scanner, last_work } psc work_buf)

    str' <- fix (\lp s -> do
                    s' <- fillWork psc s
                    p  <- scan_hdr psc
                    if p == nullPtr
                      then lp s'
                      else do smps <- B.packCString p
                              unless (length (B.split '\t' smps) == 1) . error $
                                    "Sorry, this VCF scanner only supports one individual."
                              return s'
                ) (str <* S.singleton 10)       -- trailing newline, just in case


    str'' <- fix (\lp s -> do s' <- fillWork psc s
                              r  <- skip_hdr psc
                              if r == 0 then lp s' else return s') str'

    let ctab = M.fromList $ zip (map toShort chroms) [0..]
    k psc ctab str''

foreign import ccall unsafe "vcf_scan.h scan_hdr"
    scan_hdr :: Ptr CScanner -> IO CString

foreign import ccall unsafe "vcf_scan.h skip_hdr"
    skip_hdr :: Ptr CScanner -> IO CInt

foreign import ccall unsafe "vcf_scan.h scan_vcf1"
    scan_vcf1 :: Ptr CScanner -> IO CInt

-- | Tries to fill the work buffer as far as possible.  Since this
-- function moves stuff in the working buffer around, it shouldn't be
-- called if that buffer is mostly full.
fillWork :: Ptr CScanner -> S.ByteString IO r -> IO (S.ByteString IO r)
fillWork psc str = do
    wbuf  <- #{ peek scanner, working_buffer } psc
    wnext <- #{ peek scanner, next_work } psc
    wlast <- #{ peek scanner, last_work } psc
    moveBytes wbuf wnext (wlast `minusPtr` wnext)
    (#{ poke scanner, next_work } psc wbuf)

    let wlast' = wbuf `plusPtr` (wlast `minusPtr` wnext)
    fill wlast' (wbuf `plusPtr` workbuffer_size `minusPtr` wlast') str

  where
    -- Keep copying until the work area is full.
    fill pb l = S.nextChunk >=> \case
        Left r -> do (#{ poke scanner, last_work } psc pb)
                     return $ pure r
        Right (c,cs)
            | B.length c > l
                -> do B.unsafeUseAsCString c $ \pc ->
                            copyBytes pb pc l
                      (#{ poke scanner, last_work } psc (pb `plusPtr` l))
                      return $ S.chunk (B.drop l c) >> cs

            | otherwise
                -> do B.unsafeUseAsCString c $ \pc ->
                            copyBytes pb pc (B.length c)
                      fill (pb `plusPtr` B.length c) (l - B.length c) cs


getVariant :: Ptr CScanner -> HashMap ShortByteString Int
           -> S.ByteString IO r -> IO (Either r (RawVariant, S.ByteString IO r))
getVariant psc ctab str = scan_vcf1 psc >>= go str
  where
    go s 0 = fillWork psc s >>= S.nextChunk >>= \case
                Left    r    -> return $ Left r
                Right (c,cs) -> getVariant psc ctab (S.chunk c >> cs)

    go s _ = variant s <$> getChrom
                       <*> #{ peek scanner, pos } psc
                       <*> getVars
                       <*> getGt

    getChrom = do a <- #{ peek scanner, refseq } psc
                  e <- #{ peek scanner, erefseq } psc
                  c <- createFromPtr a (e `minusPtr` a)
                  return $ M.lookupDefault (-1) c ctab

    getVars = do a <- #{ peek scanner, alleles } psc
                 e <- #{ peek scanner, ealleles } psc
                 B.packCStringLen (a, e `minusPtr` a)

    getGt = peekByteOff psc (#{ offset scanner, gts })

    variant :: S.ByteString IO r -> Int -> Word32 -> Bytes -> Word16
            -> Either r (RawVariant,S.ByteString IO r)
    variant s r p vs gt = Right (v,s)
      where !v = RawVariant r (fromIntegral p) vs gt


