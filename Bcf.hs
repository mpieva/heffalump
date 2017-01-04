{-# LANGUAGE OverloadedStrings, RecordWildCards #-}
module Bcf ( readBcf, decodeBcf ) where

-- Minimalistic BCF reader.  We only care about the genotype of the
-- first individual.  Should allow for a handful of shortcuts...

import BasePrelude
import Util ( decomp )
import VcfScan ( RawVariant(..), hashChrom )

import Foreign.Storable ( peekByteOff )

import qualified Data.ByteString                 as B
import qualified Data.ByteString.Unsafe          as B
import qualified Data.ByteString.Lazy            as L
import qualified Data.ByteString.Lazy.Char8      as C
import qualified Data.Vector.Unboxed             as V

readBcf :: FilePath -> IO [ RawVariant ]
readBcf fp = decodeBcf . decomp <$> L.readFile fp

-- Skip over header, for now we won't parse it.  Then scan records.
-- Missing: - decode the GT field (*should* be the first in the indiv part)
--          - figure out the variant alleles
decodeBcf :: L.ByteString -> [ RawVariant ]
decodeBcf str0
    | "BCF\2" == L.take 4 str0 = let l_text     = slow_word32 (L.drop 5 str0)
                                     (hdr,body) = L.splitAt l_text $ L.drop 9 str0
                                 in case L.toChunks body of
                                     [  ] -> []
                                     s:ss -> getvars (mkSymTab hdr) ss s
    | otherwise = error "not a BCFv2 file"

  where
    slow_word32 s = fromIntegral (L.index s 0) .|.
                    fromIntegral (L.index s 1) `shiftL` 8 .|.
                    fromIntegral (L.index s 2) `shiftL` 16 .|.
                    fromIntegral (L.index s 3) `shiftL` 24

    getvars :: V.Vector Int -> [B.ByteString] -> B.ByteString -> [RawVariant]
    getvars tab strs str
        | B.length str < 32                 = case strs of (s:ss) -> getvars tab ss (B.append str s)
                                                           [    ] | B.null str -> []
                                                                  | otherwise  -> error "Short record."
        | B.length str < fromIntegral l_tot = case strs of (s:ss) -> getvars tab ss (B.append str s)
        | otherwise = unsafeDupablePerformIO (B.unsafeUseAsCString str $ \p -> do
                          !l_shared <- fromIntegral        <$> (peekByteOff p  0 :: IO Word32)
                          !refid    <- fromIntegral        <$> (peekByteOff p  8 :: IO Word32)
                          !rv_pos   <- fromIntegral . succ <$> (peekByteOff p 12 :: IO Word32)
                          !n_allls  <- (`shiftR` 16)       <$> (peekByteOff p 24 :: IO Word32)
                          let !rv_vars = get_als (fromIntegral n_allls) (B.drop 32 str)
                              !rv_gt   = get_gts (B.drop (l_shared + 8) str)
                          return $! RawVariant{ rv_chrom = tab V.! refid, .. })
                      : getvars tab strs (B.drop (fromIntegral l_tot) str)
      where
        l_tot = unsafeDupablePerformIO $ B.unsafeUseAsCString str $ \p -> do
                    a <- peekByteOff p 0
                    b <- peekByteOff p 4
                    return (8 + a + b :: Word32)


    -- skip over variant ID, then get alleles
    get_als :: Int -> B.ByteString -> B.ByteString
    get_als n !s = let !sk = case B.index s 0 of
                                0xF7 -> case B.index s 1 of
                                            0x01 -> 3 + fromIntegral (B.index  s 2) -- should be plenty
                                            0x02 -> 4 + fromIntegral (indexW16 s 2) -- but isn't :(
                                            0x03 -> 6 + fromIntegral (indexW32 s 2) -- but isn't :(
                                            x -> error $ "Huh? " ++ show x
                                tp | tp .&. 0xF == 7 -> 1 + fromIntegral (tp `shiftR` 4)
                                _                    -> error "string expected"
                   in get_als' [] n (B.drop sk s)

    get_als' :: [B.ByteString] -> Int -> B.ByteString -> B.ByteString
    get_als' acc 0 !_ = B.intercalate "," $ reverse acc
    get_als' acc n !s = let (!sk,!ln) = case B.index s 0 of
                                            0xF7 -> case B.index s 1 of
                                                        0x01 -> (3, fromIntegral $ B.index  s 2) -- should be plenty
                                                        0x02 -> (4, fromIntegral $ indexW16 s 2) -- but isn't :(
                                                        0x03 -> (6, fromIntegral $ indexW32 s 2) -- but isn't :(
                                                        x -> error $ "Huh? " ++ show x
                                            tp | tp .&. 0xF == 7 -> (1, fromIntegral $ tp `shiftR` 4)
                                            _                    -> error "string expected"
                        in get_als' (B.take ln (B.drop sk s) : acc) (n-1) (B.drop (sk+ln) s)


    get_gts :: B.ByteString -> Word16
    get_gts str = let !ks = case B.index str 0 of 1 -> 2; 2 -> 3; 3 -> 5; _ -> error "WTF?"  -- key size, value ignored
                      !tp_byte = B.index str ks
                    -- we support haploid and diploid, and Word8 and Word16
                  in case tp_byte of
                    0x11 -> 0xFF00 .|. fromIntegral (B.index str (ks+1))
                    0x12 -> 0xFF00 .|. indexW16 str (ks+1)

                    0x21 -> let x = fromIntegral $ B.index str (ks+1)
                                y = fromIntegral $ B.index str (ks+2)
                            in y `shiftL` 8 .|. x

                    0x22 -> let x = indexW16 str (ks+1)
                                y = indexW16 str (ks+3)
                            in y `shiftL` 8 .|. x

                    _    -> error $ "only haploid or diploid calls are supported " ++ showHex tp_byte []

    indexW16 :: B.ByteString -> Int -> Word16
    indexW16 s i = unsafeDupablePerformIO $ B.unsafeUseAsCString s $ \p -> peekByteOff p i

    indexW32 :: B.ByteString -> Int -> Word32
    indexW32 s i = unsafeDupablePerformIO $ B.unsafeUseAsCString s $ \p -> peekByteOff p i


mkSymTab :: L.ByteString -> V.Vector Int
mkSymTab = V.fromList . mapMaybe parse . C.lines
  where
    parse l = case C.splitAt (C.length key) l of
        (u,v) | u == key  -> Just $! hashChrom (l2s $ C.takeWhile (/=',') v)
              | otherwise -> Nothing
    key = "##contig=<ID="
    l2s = B.concat . L.toChunks

