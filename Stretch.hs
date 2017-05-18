module Stretch (
        NucCode(..),
        Stretch(..),
        decode_dip,
        decode_hap,
        main_dumppatch
    ) where

import BasePrelude
import System.IO                        ( stderr )

import qualified Data.ByteString.Char8          as B
import qualified Data.ByteString.Internal       as B
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.ByteString.Lazy           as LB
import qualified Data.ByteString.Unsafe         as BB

import Util ( decomp )

-- ^ A genome is encoded by taking the difference to the reference and
-- run-length coding the result.

newtype NucCode = NucCode Word8 deriving Eq
instance Show NucCode where show = (:[]) . tr


-- All stretch lengths are even, so we pre-divide them by two.
data Stretch = Ns   !Int64            Stretch
             | Eqs  !Int64            Stretch
             | Eqs1 !Int64            Stretch
             | Chrs !NucCode !NucCode Stretch
             | Break                  Stretch
             | Done
    deriving Show

debugStretch :: Stretch -> IO ()
debugStretch = debugStretch' 0 0

debugStretch' :: Int -> Int64 -> Stretch -> IO ()
debugStretch' _ _  Done        = return ()
debugStretch' c i (Break    s) = do putStrLn $ shows (c,i) "\tBreak" ;                    debugStretch' (c+1) 0 s
debugStretch' c i (Ns     n s) = do putStrLn $ shows (c,i) "\tNs   " ++ show n ;          debugStretch' c (i+2*n) s
debugStretch' c i (Eqs    n s) = do putStrLn $ shows (c,i) "\tEqs  " ++ show n ;          debugStretch' c (i+2*n) s
debugStretch' c i (Eqs1   n s) = do putStrLn $ shows (c,i) "\tEqs1 " ++ show n ;          debugStretch' c (i+2*n) s
debugStretch' c i (Chrs x y s) = do putStrLn $ shows (c,i) "\tChrs " ++ [tr x,' ',tr y] ; debugStretch' c (i+2) s

-- | Main decoder.  Switches behavior based on header.
{-# DEPRECATED decode "Switch to Lumps" #-}
decode :: L.ByteString -> Stretch
decode str | "HEF\0"   `L.isPrefixOf` str = decode_dip (L.drop 4 str)
           | "HEF\1"   `L.isPrefixOf` str = decode_hap (L.drop 4 str)
           | "HEF\3"   `L.isPrefixOf` str = error "I can't handle this new format."

           -- legacy files (could start with anything, but
           -- in practice start with (Ns 5000) or (Ns 5001)
           | 176:113:2:_ <- LB.unpack str = decode_dip str
           | 177:113:2:_ <- LB.unpack str = decode_dip str

           | otherwise                    = error "File format not recognixed."

-- We need to store 11 different symbols: "NACGTMRWSYK"  Long
-- matching stretches and long stretches of no call ('N') shall be
-- run-length encoded.  Since 11^2 == 121 < 128, we pack two symbols and
-- a signal bit into one byte.  Also, to avoid messing with nybbles, all
-- counts shall be multiples of two.
--
-- 0xxxxxxxb: two symbols are encoded as x+11*y where x and y are the
--            positions of the symbols in the string above.
-- 10000000b: 'Break'
-- 100xxxxxb: 1-31 Ns
-- 1010xxxxb: together with the next byte, 0-8191 Ns
-- 10110xxxb: together with the next two bytes, 0-1048575 Ns
-- 101110xxb: together with the next three bytes, 0-134217727 Ns
-- 1011110xb: together with the next four bytes, 0-17179869183 Ns
--            (this shouldn't be needed)
-- 110xxxxxb: 0-31 matches
-- 1110xxxxb: together with the next byte, 0-8191 matches
-- 11110xxxb: together with the next two bytes, 0-1048575 matches
--            (this will rarely be needed)
-- 111110xxb: together with the next three bytes, 0-134217727 matches
--            (this shouldn't be needed)
-- 1111110xb: together with the next four bytes, 0-17179869183 matches
--            (this shouldn't be needed)

-- Decoding produces one stretch until it finds the end marker or runs
-- out of input.  We shall get one 'Stretch' for each chromosome.

decode_dip, decode_hap :: L.ByteString -> Stretch
decode_dip = decode_v0 NucCode Eqs
decode_hap = decode_v0 nc Eqs1
  where
    nc x |    x == 0 = NucCode 0
         |    x <= 4 = NucCode $ x+10
         | otherwise = NucCode x


-- | Decode for original format, parameterized so it can be used for
-- haploid and diploid individuals (aka "high coverage" and "low
-- coverage" input).
{-# DEPRECATED decode_v0 "Switch to Lumps" #-}
decode_v0 :: (Word8 -> NucCode) -> (Int64 -> Stretch -> Stretch) -> L.ByteString -> Stretch
decode_v0 nucCode eqs = go
  where
    go s
      | L.null  s = Done
      | otherwise = case LB.head s of
        w | w < 0x80 -> let (y,x) = w `divMod` 11 in Chrs (nucCode x) (nucCode y) $ go (LB.tail s)
          | w ==0x80 -> Break $ go (LB.tail s)

          | w < 0xA0 -> Ns (fromIntegral (w .&. 0x1f)) $ go (LB.tail s)
          | w < 0xB0 -> Ns (fromIntegral (w .&. 0x0f) .|.
                            fromIntegral (LB.index s 1) `shiftL` 4) $ go (LB.drop 2 s)
          | w < 0xB8 -> Ns (fromIntegral (w .&. 0x07) .|.
                            fromIntegral (LB.index s 1) `shiftL` 3 .|.
                            fromIntegral (LB.index s 2) `shiftL` 11) $ go (LB.drop 3 s)
          | w < 0xBC -> Ns (fromIntegral (w .&. 0x03) .|.
                            fromIntegral (LB.index s 1) `shiftL` 2 .|.
                            fromIntegral (LB.index s 2) `shiftL` 10 .|.
                            fromIntegral (LB.index s 3) `shiftL` 18) $ go (LB.drop 4 s)
          | w < 0xBE -> Ns (fromIntegral (w .&. 0x01) .|.
                            fromIntegral (LB.index s 1) `shiftL` 1 .|.
                            fromIntegral (LB.index s 2) `shiftL` 9 .|.
                            fromIntegral (LB.index s 3) `shiftL` 17 .|.
                            fromIntegral (LB.index s 4) `shiftL` 25) $ go (LB.drop 5 s)
          | w < 0xC0 -> error $ "WTF?! (too many Ns) " ++ show w

          | w < 0xE0 -> eqs (fromIntegral (w .&. 0x1f)) $ go (L.tail s)
          | w < 0xF0 -> eqs (fromIntegral (w .&. 0x0f) .|.
                            fromIntegral (LB.index s 1) `shiftL` 4) $ go (LB.drop 2 s)
          | w < 0xF8 -> eqs (fromIntegral (w .&. 0x07) .|.
                            fromIntegral (LB.index s 1) `shiftL` 3 .|.
                            fromIntegral (LB.index s 2) `shiftL` 11) $ go (LB.drop 3 s)
          | otherwise -> error $ "WTF?! (too many matches) " ++ show w


-- | We store diploid calls.  For this we need 11 codes:  no call(1),
-- the bases(4), heterozygotes(6).  If we also want to support haploid
-- calls (for Treemix, this does make a difference!), we need four more.
iupac_chars :: B.ByteString
iupac_chars = "NACGTMRWSYKacgtN"

{-# INLINE tr #-}
tr :: NucCode -> Char
tr (NucCode w) = B.w2c . BB.unsafeIndex iupac_chars . fromIntegral $ w .&. 0xF

main_dumppatch :: [String] -> IO ()
main_dumppatch [inf] = debugStretch . decode . decomp =<< L.readFile inf
main_dumppatch     _ = B.hPut stderr "Usage: dumppatch [foo.hef]\n"

