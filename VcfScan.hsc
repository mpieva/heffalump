{-# LANGUAGE OverloadedStrings #-}
module VcfScan where

import BasePrelude
import Foreign
import Foreign.C
import System.IO

import qualified Data.ByteString        as BB
import qualified Data.ByteString.Char8  as B
import qualified Data.Vector            as V

#include "zlib.h"
#include "vcf-scan.h"

-- | Supports one individual only.  (Project your VCF if you have more
-- than one.)
data RawVariant = RawVariant {
    rv_chrom :: !Int,                                -- chromosome number, hashed
    rv_pos   :: !Int,                                -- position (1-based)
    rv_vars  :: !B.ByteString,                       -- vars, comma separated, incl. ref
    rv_gt    :: !Word16 }                            -- genotype
        deriving (Show, Eq)

-- The way we hash chromosome names.
hashChrom :: B.ByteString -> Int
hashChrom = BB.foldl' (\a b -> 49 * a + fromIntegral b) 0

data CScanner

data Scanner = Scanner { fhandle :: !Handle
                       , filename :: FilePath
                       , c_scanner :: !(ForeignPtr CScanner) }

withScanner :: (Ptr CScanner -> Handle -> IO a) -> Scanner -> IO a
withScanner k sc = withForeignPtr (c_scanner sc) $ \p -> k p (fhandle sc)

ibuffer_size, workbuffer_size :: Int
ibuffer_size = 2000000
workbuffer_size = 127*1024

type Samples = V.Vector B.ByteString

-- | Initialize scanner.  Allocates an input and a working buffer, fills
-- the working buffer from the input, scans the header and sets up the
-- c structure.
initVcf :: FilePath -> IO Scanner
initVcf fp = do
    psc <- mallocBytes (#{const sizeof(scanner)})
    hdl <- openFile fp ReadMode

    raw_buf <- mallocBytes ibuffer_size
    (#{ poke scanner, input_buffer } psc raw_buf)
    (#{ poke scanner, next_input } psc raw_buf)
    (#{ poke scanner, last_input } psc raw_buf)

    work_buf <- mallocBytes workbuffer_size
    (#{ poke scanner, working_buffer } psc work_buf)
    (#{ poke scanner, next_work } psc work_buf)
    (#{ poke scanner, last_work } psc work_buf)

    smps <- let loop = do _ <- fillWork psc hdl
                          p <- scan_hdr psc
                          if p == nullPtr then loop
                                          else B.packCString p
            in B.split '\t' <$> loop

    unless (length smps == 1) . error $
        "Sorry, this VCF scanner only supports one individual."

    let loop = do  _ <- fillWork psc hdl
                   r <- skip_hdr psc
                   when (r == 0) loop
        in loop

    (#{ poke scanner, nsmp } psc (length smps))
    Scanner hdl fp <$> newForeignPtr free_scanner psc

foreign import ccall unsafe "vcf_scan.h &free_scanner"
    free_scanner :: FunPtr (Ptr CScanner -> IO ())

foreign import ccall unsafe "vcf_scan.h scan_hdr"
    scan_hdr :: Ptr CScanner -> IO CString

foreign import ccall unsafe "vcf_scan.h skip_hdr"
    skip_hdr :: Ptr CScanner -> IO CInt

foreign import ccall unsafe "vcf_scan.h scan_vcf1"
    scan_vcf1 :: Ptr CScanner -> IO CInt

-- | Try to fill raw buffer.  Reads as much as it can from the handle.
fillRaw :: Ptr CScanner -> Handle -> IO ()
fillRaw psc hdl = do
    ibuf <- #{ peek scanner, input_buffer } psc
    ibeg <- #{ peek scanner, next_input } psc
    iend <- #{ peek scanner, last_input } psc
    moveBytes ibuf ibeg (iend `minusPtr` ibeg)        -- move remainder to begin of buffer
    l <- hGetBuf hdl (ibuf `plusPtr` (iend `minusPtr` ibeg))
                     (ibuffer_size - (iend `minusPtr` ibeg))
    (#{ poke scanner, next_input } psc ibuf)
    (#{ poke scanner, last_input } psc (ibuf `plusPtr` (iend `minusPtr` ibeg + l)))


-- | Tries to fill the work buffer.  We try to decompress stuff first.
-- If that fails, we call 'fillRaw', then try again.
-- At end of input, this naturally results in nothing happening.  Since
-- this function moves stuff in the working buffer around, it shouldn't
-- be called if that buffer is mostly full.  Returns the number of bytes
-- added to the working buffer.
fillWork :: Ptr CScanner -> Handle -> IO Int
fillWork psc hdl = do
    wbuf <- #{ peek scanner, working_buffer } psc
    wnext <- #{ peek scanner, next_work } psc
    wlast <- #{ peek scanner, last_work } psc
    -- hPutStrLn stderr $ "moving " ++ show (wlast `minusPtr` wnext) ++ " bytes in working buffer."
    moveBytes wbuf wnext (wlast `minusPtr` wnext)
    (#{ poke scanner, next_work } psc wbuf)
    (#{ poke scanner, last_work } psc (wbuf `plusPtr` (wlast `minusPtr` wnext)))

    n1 <- tryDecompressBlocks psc
    if n1 == (-1) then do fillRaw psc hdl ; tryDecompressBlocks psc
                  else return n1

-- | Decompresses as long as at least one complete block is available
-- and the working buffer has space to decompress it.  Returns the
-- number of bytes produces, or -1 if more input is needed.
tryDecompressBlocks :: Ptr CScanner -> IO Int
tryDecompressBlocks psc = do
    ibeg <- #{ peek scanner, next_input } psc
    iend <- #{ peek scanner, last_input } psc
    if iend `minusPtr` ibeg < 26                 -- BGZF header size?
      then return (-1)
      else do
        66 <- peekByteOff ibeg 12 :: IO Word8      -- BC tag
        67 <- peekByteOff ibeg 13 :: IO Word8
        bsize <- peekByteOff ibeg 16 :: IO Word16
        if iend `minusPtr` ibeg < fromIntegral bsize + 1
          then return (-1)
          else do
            crc <- peekByteOff ibeg (fromIntegral bsize - 7) :: IO Word32
            ucompsize <- peekByteOff ibeg (fromIntegral bsize - 3) :: IO Word32
            wbuf <- #{ peek scanner, working_buffer } psc
            wlast <- #{ peek scanner, last_work } psc
            if workbuffer_size - (wlast `minusPtr` wbuf) < fromIntegral ucompsize
              then return 0
              else do
                -- call to zlib here
                -- hPutStrLn stderr $ "inflating " ++ show (bsize + 1) ++ " bytes into " ++ show ucompsize ++ " bytes."
                inflate (ibeg `plusPtr` 18) (fromIntegral bsize + 1) wlast
                                            (fromIntegral ucompsize) (fromIntegral crc)
                (#{ poke scanner, last_work }) psc (wlast `plusPtr` fromIntegral ucompsize)
                (#{ poke scanner, next_input }) psc (ibeg `plusPtr` (fromIntegral bsize + 1))
                (+ fromIntegral ucompsize) <$> tryDecompressBlocks psc

inflate :: Ptr Word8 -> CUInt -> Ptr Word8 -> CUInt -> CUInt -> IO ()
inflate next_in num_in next_out num_out crc =
    allocaBytes #{const sizeof(z_stream)} $ \stream -> do
    (#{poke z_stream, msg}       stream nullPtr)
    (#{poke z_stream, zalloc}    stream nullPtr)
    (#{poke z_stream, zfree}     stream nullPtr)
    (#{poke z_stream, opaque}    stream nullPtr)
    (#{poke z_stream, next_in}   stream next_in)
    (#{poke z_stream, next_out}  stream next_out)
    (#{poke z_stream, avail_in}  stream num_in)
    (#{poke z_stream, avail_out} stream num_out)
    z_check "inflateInit2" =<< c_inflateInit2 stream (-15)
    z_check "inflate" =<< c_inflate stream #{const Z_FINISH}
    z_check "inflateEnd" =<< c_inflateEnd stream

    pe <- #{peek z_stream, next_out} stream
    when (pe `minusPtr` next_out /= fromIntegral num_out) $ error "size mismatch after deflate()"

    crc0 <- c_crc32 0 nullPtr 0
    crc' <- c_crc32 crc0 next_out (fromIntegral num_out)
    when (fromIntegral crc /= crc') $ error "CRC error after deflate()"

data ZStream

{-# INLINE z_check #-}
z_check :: String -> CInt -> IO ()
z_check msg c = when (c /= #{const Z_OK} && c /= #{const Z_STREAM_END}) $
                   error $ msg ++ " failed: " ++ show c

c_inflateInit2 :: Ptr ZStream -> CInt -> IO CInt
c_inflateInit2 z a = withCAString #{const_str ZLIB_VERSION} $ \versionStr ->
    c_inflateInit2_ z a versionStr (#{const sizeof(z_stream)} :: CInt)

foreign import ccall unsafe "zlib.h inflateInit2_" c_inflateInit2_ ::
    Ptr ZStream -> CInt -> Ptr CChar -> CInt -> IO CInt

foreign import ccall unsafe "zlib.h inflate" c_inflate ::
    Ptr ZStream -> CInt -> IO CInt

foreign import ccall unsafe "zlib.h inflateEnd" c_inflateEnd ::
    Ptr ZStream -> IO CInt

foreign import ccall unsafe "zlib.h crc32" c_crc32 ::
    CULong -> Ptr Word8 -> CUInt -> IO CULong


getVariant :: Scanner -> IO (Maybe RawVariant)
getVariant sc = withScanner go sc
  where
    go psc hdl = scan_vcf1 psc >>= go' psc hdl

    go' psc hdl 0 = do
        n <- fillWork psc hdl
        wb <- #{ peek scanner, next_work } psc
        we <- #{ peek scanner, last_work } psc
        if n <= 0 && wb == (we :: Ptr Word8)
            then return Nothing
            else do -- maybe newline is missing?
                    when (n <= 0) $ pokeByteOff we (-1) (10::Word8)
                    go psc hdl

    go' psc  _  _ = variant <$> #{ peek scanner, refseq } psc
                            <*> #{ peek scanner, pos } psc
                            <*> getVars psc
                            <*> getGt psc

    getVars p = do a <- #{ peek scanner, alleles } p
                   e <- #{ peek scanner, ealleles } p
                   B.packCStringLen (a, e `minusPtr` a)

    getGt p = peekByteOff p (#{ offset scanner, gts })

    variant :: Word16 -> Word32 -> B.ByteString -> Word16 -> Maybe RawVariant
    variant r p vs gt = Just $! RawVariant (fromIntegral r) (fromIntegral p) vs gt


