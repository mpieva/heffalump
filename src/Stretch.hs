{-# LANGUAGE Rank2Types #-}
module Stretch (
        NucCode(..),
        Stretch(..),
        decode_dip,
        decode_hap,
        main_dumppatch
    ) where

import Bio.Prelude               hiding ( Ns )
import Streaming
import System.IO                        ( withFile, IOMode(..), stderr )

import qualified Data.ByteString.Char8          as B
import qualified Data.ByteString.Internal       as B
import qualified Data.ByteString.Streaming      as S
import qualified Data.ByteString.Unsafe         as BB
import qualified Streaming.Prelude              as Q

import Util ( decomp )

-- ^ A genome is encoded by taking the difference to the reference and
-- run-length coding the result.

newtype NucCode = NucCode Word8 deriving Eq
instance Show NucCode where show = (:[]) . tr


-- All stretch lengths are even, so we pre-divide them by two.
data Stretch = Ns   !Int64
             | Eqs  !Int64
             | Eqs1 !Int64
             | Chrs !NucCode !NucCode
             | Break
    deriving Show

debugStretch :: MonadIO m => Stream (Of Stretch) m r -> m r
debugStretch = debugStretch' 0 0

debugStretch' :: MonadIO m => Int -> Int64 -> Stream (Of Stretch) m r -> m r
debugStretch' c i = Q.next >=> \case
    Left r -> return r
    Right (Break   ,s) -> do liftIO . putStrLn $ shows (c,i) "\tBreak" ;                    debugStretch' (c+1) 0 s
    Right (Ns     n,s) -> do liftIO . putStrLn $ shows (c,i) "\tNs   " ++ show n ;          debugStretch' c (i+2*n) s
    Right (Eqs    n,s) -> do liftIO . putStrLn $ shows (c,i) "\tEqs  " ++ show n ;          debugStretch' c (i+2*n) s
    Right (Eqs1   n,s) -> do liftIO . putStrLn $ shows (c,i) "\tEqs1 " ++ show n ;          debugStretch' c (i+2*n) s
    Right (Chrs x y,s) -> do liftIO . putStrLn $ shows (c,i) "\tChrs " ++ [tr x,' ',tr y] ; debugStretch' c (i+2) s


-- | Main decoder.  Switches behavior based on header.
{-# DEPRECATED decode "Switch to Lumps" #-}
decode :: Monad m => S.ByteString m r -> Stream (Of Stretch) m r
decode str = do hd :> tl <- lift . S.toStrict $ S.splitAt 4 str
                case B.take 3 hd :> B.drop 3 hd of
                    "HEF" :> "\0" -> decode_dip tl
                    "HEF" :> "\1" -> decode_hap tl
                    "HEF" :> "\3" -> error "I can't handle this new format."

                    -- legacy files could start with anything, but
                    -- in practice start with (Ns 5000) or (Ns 5001)
                    "\176\113\2" :> _ -> decode_dip $ S.chunk hd >> tl
                    "\177\113\2" :> _ -> decode_dip $ S.chunk hd >> tl
                    _                 -> error "File format not recognixed."


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

decode_dip, decode_hap :: Monad m => S.ByteString m r -> Stream (Of Stretch) m r
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
decode_v0 :: Monad m => (Word8 -> NucCode) -> (Int64 -> Stretch)
                     -> S.ByteString m r -> Stream (Of Stretch) m r
decode_v0 nucCode eqs = go
  where
    go = lift . S.nextByte >=> \case
       Left    r     -> pure r
       Right (w,s')
          | w < 0x80 -> let (y,x) = w `divMod` 11 in Chrs (nucCode x) (nucCode y) `Q.cons` go s'
          | w ==0x80 -> Break `Q.cons` go s'

          | w < 0xA0 -> Ns (fromIntegral (w .&. 0x1f)) `Q.cons` go s'
          | w < 0xB0 -> do Right (w1,s1) <- lift $ S.nextByte s'
                           Ns (fromIntegral (w .&. 0x0f) .|.
                               fromIntegral w1 `shiftL` 4) `Q.cons` go s1
          | w < 0xB8 -> do Right (w1,s1) <- lift $ S.nextByte s'
                           Right (w2,s2) <- lift $ S.nextByte s1
                           Ns (fromIntegral (w .&. 0x07) .|.
                               fromIntegral w1 `shiftL` 3 .|.
                               fromIntegral w2 `shiftL` 11) `Q.cons` go s2
          | w < 0xBC -> do Right (w1,s1) <- lift $ S.nextByte s'
                           Right (w2,s2) <- lift $ S.nextByte s1
                           Right (w3,s3) <- lift $ S.nextByte s2
                           Ns (fromIntegral (w .&. 0x03) .|.
                               fromIntegral w1 `shiftL` 2 .|.
                               fromIntegral w2 `shiftL` 10 .|.
                               fromIntegral w3 `shiftL` 18) `Q.cons` go s3
          | w < 0xBE -> do Right (w1,s1) <- lift $ S.nextByte s'
                           Right (w2,s2) <- lift $ S.nextByte s1
                           Right (w3,s3) <- lift $ S.nextByte s2
                           Right (w4,s4) <- lift $ S.nextByte s3
                           Ns (fromIntegral (w .&. 0x01) .|.
                               fromIntegral w1 `shiftL` 1 .|.
                               fromIntegral w2 `shiftL` 9 .|.
                               fromIntegral w3 `shiftL` 17 .|.
                               fromIntegral w4 `shiftL` 25) `Q.cons` go s4
          | w < 0xC0 -> error $ "WTF?! (too many Ns) " ++ show w

          | w < 0xE0 -> eqs (fromIntegral (w .&. 0x1f)) `Q.cons` go s'
          | w < 0xF0 -> do Right (w1,s1) <- lift $ S.nextByte s'
                           eqs (fromIntegral (w .&. 0x0f) .|.
                                fromIntegral w1 `shiftL` 4) `Q.cons` go s1
          | w < 0xF8 -> do Right (w1,s1) <- lift $ S.nextByte s'
                           Right (w2,s2) <- lift $ S.nextByte s1
                           eqs (fromIntegral (w .&. 0x07) .|.
                                fromIntegral w1 `shiftL` 3 .|.
                                fromIntegral w2 `shiftL` 11) `Q.cons` go s2
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
main_dumppatch [inf] = withFile inf ReadMode $
                            debugStretch . decode . decomp . S.fromHandle
main_dumppatch     _ = B.hPut stderr "Usage: dumppatch [foo.hef]\n"

