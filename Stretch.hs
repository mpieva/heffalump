module Stretch (
    NucCode(..),
    Stretch(..),
    debugStretch,
    catStretches,
    encode_v0,
    encode_dip,
    encode_hap,
    encode_mix,
    decode,
    diff,
    Frag(..),
    patch,
    code,
    tr,
    iupac_chars
    ) where

import Data.Bits
import Data.ByteString.Builder          ( word8, Builder, byteString )
import Data.Int                         ( Int64 )
import Data.Monoid
import Data.Word                        ( Word8 )

import qualified Data.ByteString.Char8          as B
import qualified Data.ByteString.Internal       as B
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.ByteString.Lazy           as LB
import qualified Data.ByteString.Unsafe         as B
import qualified Data.Nullable                  as Null
import qualified Data.NullPoint                 as Null

import Util ( low )

-- ^ A genome is encoded by taking the difference to the reference and
-- run-length coding the result.

newtype NucCode = NucCode Word8 deriving Eq
instance Show NucCode where show = (:[]) . tr

data Patch = NewNs   !Int64 Patch
           | NewEqs  !Int64 Patch
           | NewEqs1 !Int64 Patch
           | Call  !NucCode Patch
           | NewBreak       Patch
           | NewDone
    deriving Show

stretchToPatch :: Stretch -> Patch
stretchToPatch (Ns     n k) = NewNs   n $ stretchToPatch k
stretchToPatch (Eqs    n k) = NewEqs  n $ stretchToPatch k
stretchToPatch (Eqs1   n k) = NewEqs1 n $ stretchToPatch k
stretchToPatch (Chrs a b k) = Call a $ Call b $ stretchToPatch k
stretchToPatch (Break    k) = NewBreak $ stretchToPatch k
stretchToPatch  Done        = NewDone

-- All stretch lengths are even, so we pre-divide them by two.
data Stretch = Ns   !Int64            Stretch
             | Eqs  !Int64            Stretch
             | Eqs1 !Int64            Stretch
             | Chrs !NucCode !NucCode Stretch
             | Break                  Stretch
             | Done
    deriving Show

instance Null.NullPoint (Stretch -> Stretch) where empty = id
instance Null.Nullable  (Stretch -> Stretch) where
    nullC f = case f Done of Done -> True ; _ -> False


debugStretch :: Stretch -> IO ()
debugStretch = debugStretch' 0 0

debugStretch' :: Int -> Int64 -> Stretch -> IO ()
debugStretch' c i  Done        = return () --    putStrLn $ shows (c,i) "\tDone"
debugStretch' c i (Break    s) = do putStrLn $ shows (c,i) "\tBreak" ;                    debugStretch' (c+1) 0 s
debugStretch' c i (Ns     n s) = do putStrLn $ shows (c,i) "\tNs   " ++ show n ;          debugStretch' c (i+2*n) s
debugStretch' c i (Eqs    n s) = do putStrLn $ shows (c,i) "\tEqs  " ++ show n ;          debugStretch' c (i+2*n) s
debugStretch' c i (Eqs1   n s) = do putStrLn $ shows (c,i) "\tEqs1 " ++ show n ;          debugStretch' c (i+2*n) s
debugStretch' c i (Chrs x y s) = do putStrLn $ shows (c,i) "\tChrs " ++ [tr x,' ',tr y] ; debugStretch' c (i+2) s

catStretches :: [Stretch -> Stretch] -> Stretch
catStretches = foldr (\a b -> a $ Break b) Done


-- | Main decoder.  Switches behavior based on header.
decode :: L.ByteString -> Stretch
decode str | "HEF\0" `L.isPrefixOf` str = decode_dip (L.drop 4 str)
           | "HEF\1" `L.isPrefixOf` str = decode_hap (L.drop 4 str)
           | "HEF\2" `L.isPrefixOf` str = decode_mix (L.drop 4 str)
           | otherwise                  = decode_dip (L.drop 4 str) -- error "Format not recognixed."

encode_dip :: Stretch -> Builder
encode_dip s = byteString "HEF\0" <> encode_v0 s

encode_hap :: Stretch -> Builder
encode_hap s = byteString "HEF\1" <> encode_v0 s

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

-- | Version 0 encoding.  Deals with diploid calls, stretches of no call,
-- stretches of matches.
encode_v0 :: Stretch -> Builder
encode_v0 (Chrs (NucCode x) (NucCode y) k)
    | 0 <= x && x < 15 && 0 <= y && y < 15 = let x' = if x >= 11 then x-10 else x
                                                 y' = if y >= 11 then y-10 else y
                                             in word8 (fromIntegral $ x' + 11 * y') <> encode_v0 k
    | otherwise = error $ "shouldn't happen: Chrs " ++ show x ++ " " ++ show y

encode_v0 (Ns x k) | x <= 0 = error $ "WTF?  Backwards stretch?! " ++ show x
                | x < 0x20        = word8 (0x80 .|. fromIntegral x) <> encode_v0 k
                | x < 0x1000      = word8 (0xA0 .|. fromIntegral x .&. 0xf)
                                 <> word8 (fromIntegral (x `shiftR` 4)) <> encode_v0 k
                | x < 0x80000     = word8 (0xB0 .|. fromIntegral x .&. 0x7)
                                 <> word8 (fromIntegral (x `shiftR` 3))
                                 <> word8 (fromIntegral (x `shiftR` 11)) <> encode_v0 k
                | x < 0x4000000   = word8 (0xB8 .|. fromIntegral x .&. 0x3)
                                 <> word8 (fromIntegral (x `shiftR` 2))
                                 <> word8 (fromIntegral (x `shiftR` 10))
                                 <> word8 (fromIntegral (x `shiftR` 18)) <> encode_v0 k
                | x < 0x200000000 = word8 (0xBC .|. fromIntegral x .&. 0x1)
                                 <> word8 (fromIntegral (x `shiftR` 1))
                                 <> word8 (fromIntegral (x `shiftR` 9))
                                 <> word8 (fromIntegral (x `shiftR` 17))
                                 <> word8 (fromIntegral (x `shiftR` 25)) <> encode_v0 k
                | otherwise       = error $ "WTF?! (too many Ns: " ++ show (2*x) ++ ")"

encode_v0 (Eqs x k) | x <= 0 = error "WTF?  Backwards matches?!"
                 | x < 0x20     = word8 (0xC0 .|. fromIntegral x) <> encode_v0 k
                 | x < 0x1000   = word8 (0xE0 .|. fromIntegral x .&. 0xf)
                               <> word8 (fromIntegral (x `shiftR` 4)) <> encode_v0 k
                 | x < 0x80000  = word8 (0xF0 .|. fromIntegral x .&. 0x7)
                               <> word8 (fromIntegral (x `shiftR` 3))
                               <> word8 (fromIntegral (x `shiftR` 11)) <> encode_v0 k
                 | otherwise    = error "WTF?! (too many matches)"

encode_v0 (Eqs1 x k) = encode_v0 $ Eqs x k
encode_v0 (Break k) = word8 0x80 <> encode_v0 k
encode_v0 Done      = mempty


decode_mix :: L.ByteString -> Stretch
decode_mix = undefined -- XXX

-- | Version 2 encoding.  Deals with haploid or diploid calls, stretches of no call,
-- stretches of matches, stretches matching only one base.
--
-- 0xF0     Ns, one length byte follows
-- 0xF1     Ns, two length bytes follow
-- 0xF2     Ns, three length bytes follow
-- 0xF3     Ns, four length bytes follow
-- 0xF4     Eqs, one length byte follows
-- 0xF5     Eqs, two length bytes follow
-- 0xF6     Eqs, three length bytes follow
-- 0xF7     Eqs, four length bytes follow
-- 0xF8     Eqs1, one length byte follows
-- 0xF9     Eqs1, two length bytes follow
-- 0xFA     Eqs1, three length bytes follow
-- 0xFB     Eqs1, four length bytes follow
-- 0xFF     Break
-- 0xXY     Two chars, Y and X

encode_mix :: Stretch -> Builder
encode_mix s = byteString "HEF\2" <> go s
  where
    go (Chrs (NucCode x) (NucCode y) k)
        | 0 <= x && x < 15 && 0 <= y && y < 15 = word8 (fromIntegral $ x + 16 * y) <> go k
        | otherwise = error $ "shouldn't happen: Chrs " ++ show x ++ " " ++ show y

    go (Ns x k) | x <= 0 = error $ "WTF?  Backwards stretch?! " ++ show x
                | x < 0x100       = word8 0xF0
                                 <> word8 (fromIntegral x) <> go k
                | x < 0x10000     = word8 0xF1
                                 <> word8 (fromIntegral (x `shiftR` 0))
                                 <> word8 (fromIntegral (x `shiftR` 8)) <> go k
                | x < 0x1000000   = word8 0xF2
                                 <> word8 (fromIntegral (x `shiftR` 0))
                                 <> word8 (fromIntegral (x `shiftR` 8))
                                 <> word8 (fromIntegral (x `shiftR` 16)) <> go k
                | x < 0x100000000 = word8 0xF3
                                 <> word8 (fromIntegral (x `shiftR` 0))
                                 <> word8 (fromIntegral (x `shiftR` 8))
                                 <> word8 (fromIntegral (x `shiftR` 16))
                                 <> word8 (fromIntegral (x `shiftR` 24)) <> go k
                | otherwise       = error $ "WTF?! (too many Ns: " ++ show (2*x) ++ ")"

    go (Eqs x k) | x <= 0 = error $ "WTF?  Backwards stretch?! " ++ show x
                 | x < 0x100       = word8 0xF4
                                  <> word8 (fromIntegral x) <> go k
                 | x < 0x10000     = word8 0xF5
                                  <> word8 (fromIntegral (x `shiftR` 0))
                                  <> word8 (fromIntegral (x `shiftR` 8)) <> go k
                 | x < 0x1000000   = word8 0xF6
                                  <> word8 (fromIntegral (x `shiftR` 0))
                                  <> word8 (fromIntegral (x `shiftR` 8))
                                  <> word8 (fromIntegral (x `shiftR` 16)) <> go k
                 | x < 0x100000000 = word8 0xF7
                                  <> word8 (fromIntegral (x `shiftR` 0))
                                  <> word8 (fromIntegral (x `shiftR` 8))
                                  <> word8 (fromIntegral (x `shiftR` 16))
                                  <> word8 (fromIntegral (x `shiftR` 24)) <> go k
                 | otherwise       = error $ "WTF?! (too many Ns: " ++ show (2*x) ++ ")"

    go (Eqs1 x k) | x <= 0 = error $ "WTF?  Backwards stretch?! " ++ show x
                  | x < 0x100       = word8 0xF8
                                   <> word8 (fromIntegral x) <> go k
                  | x < 0x10000     = word8 0xF9
                                   <> word8 (fromIntegral (x `shiftR` 0))
                                   <> word8 (fromIntegral (x `shiftR` 8)) <> go k
                  | x < 0x1000000   = word8 0xFA
                                   <> word8 (fromIntegral (x `shiftR` 0))
                                   <> word8 (fromIntegral (x `shiftR` 8))
                                   <> word8 (fromIntegral (x `shiftR` 16)) <> go k
                  | x < 0x100000000 = word8 0xFB
                                   <> word8 (fromIntegral (x `shiftR` 0))
                                   <> word8 (fromIntegral (x `shiftR` 8))
                                   <> word8 (fromIntegral (x `shiftR` 16))
                                   <> word8 (fromIntegral (x `shiftR` 24)) <> go k
                  | otherwise       = error $ "WTF?! (too many Ns: " ++ show (2*x) ++ ")"

    go (Break k) = word8 0xFF <> go k
    go Done      = mempty

-- | We store diploid calls.  For this we need 11 codes:  no call(1),
-- the bases(4), heterozygotes(6).  If we also want to support haploid
-- calls (for Treemix, this does make a difference!), we need four more.
iupac_chars :: B.ByteString
iupac_chars = "NACGTMRWSYKacgtN"

{-# INLINE code #-}
code :: Char -> NucCode
code a = NucCode $ maybe 0 fromIntegral $ B.elemIndex a iupac_chars

{-# INLINE tr #-}
tr :: NucCode -> Char
tr (NucCode w) = B.w2c . B.unsafeIndex iupac_chars . fromIntegral $ w

-- We operate on two characters at a time, to maintain alignment.  The
-- output length is that of the reference rounded up to a multiple of
-- two; if necessary, the sample is paddded with Ns.
--
-- This code treats "N" and "-" the same, since they don't seem to have
-- a meaningful difference in hetfa files.  Note that correct parsing of
-- MAF files depends on this behavior!  Also, "Q" means a match to the
-- reference.
diff :: L.ByteString -> L.ByteString -> Stretch -> Stretch
diff r0 s0 done = generic r0 s0
  where
    isN  c = c == 'N' || c == 'n' || c == '-'
    eq a b = b == 'Q' || b == 'q' || low (B.c2w a) == low (B.c2w b)

    -- Scan generic strings
    generic ref smp
        -- corner cases if the reference ends
        | L.null         ref                = done
        | L.null (L.tail ref) && L.null smp = Chrs         (NucCode 0) (NucCode 0) done
        | L.null (L.tail ref)               = Chrs (code $ L.head smp) (NucCode 0) done

        -- corner cases if the sample ends
        | L.null smp                        =                                        Ns (succ (L.length ref) `shiftR` 1) done
        | L.null (L.tail smp)               = Chrs (code $ L.head smp) (NucCode 0) $ Ns       (L.length ref  `shiftR` 1) done

        -- general case, look at two characters
        | isN  x && isN  y = go_n  1 (L.drop 2 ref) (L.drop 2 smp)
        | eq u x && eq v y = go_eq 1 (L.drop 2 ref) (L.drop 2 smp)
        | otherwise        = Chrs (code x) (code y) $ generic (L.drop 2 ref) (L.drop 2 smp)
      where
        x = L.head smp
        y = L.head (L.tail smp)
        u = L.head ref
        v = L.head (L.tail ref)


    -- Scan a stretch of Ns
    go_n !n ref smp
        -- corner cases if the reference ends
        | L.null         ref                 = Ns  n    done
        | L.null (L.tail ref) && L.null smp  = Ns (succ n) done
        | L.null (L.tail ref) && isN x       = Ns (succ n) done
        | L.null (L.tail ref)                = Ns  n $ Chrs (code x) (NucCode 0) done

        -- corner cases if the sample ends
        | L.null smp                   = Ns (succ (L.length ref) `shiftR` 1 + n) done
        | L.null (L.tail smp) && isN x = Ns (succ (L.length ref) `shiftR` 1 + n) done
        | L.null (L.tail smp)          = Ns n $ Chrs (code x) (NucCode 0) done

        -- general case, look at two characters
        | isN x && isN y = go_n (succ n) (L.drop 2 ref) (L.drop 2 smp)
        | otherwise      = Ns n $ generic ref smp
      where
        x = L.head smp
        y = L.head (L.tail smp)

    -- Scan a stretch of matches
    go_eq !n ref smp
        -- corner cases if the reference ends
        | L.null         ref                 = Eqs n done
        | L.null (L.tail ref) && L.null smp  = Eqs n $ Chrs (NucCode 0) (NucCode 0) done
        | L.null (L.tail ref) && eq u x      = Eqs (succ n) done
        | L.null (L.tail ref)                = Eqs n $ Chrs (code $ L.head smp) (NucCode 0) done

        -- corner cases if the sample ends
        | L.null smp || L.null (L.tail smp)  = Eqs n $ generic ref smp

        -- general case, look at two characters
        | eq u x && eq v y = go_eq (succ n) (L.drop 2 ref) (L.drop 2 smp)
        | otherwise        = Eqs n $ generic ref smp
      where
        x = L.head smp
        y = L.head (L.tail smp)
        u = L.head ref
        v = L.head (L.tail ref)

data Frag = Short !Char Frag | Long !L.ByteString Frag | Term Stretch

-- | Applies a 'Stretch' to a sequence.  Returns fragments of the
-- modified sequence (this should restore the sample) and finally the
-- remaining 'Stretch'.  The 'Break' is consumed if present.
patch :: L.ByteString -> Stretch -> Frag
patch ref (Break s) = Long ref $ Term s
patch ref s@Done    = Long ref $ Term s

patch ref s | L.null ref = clear s

patch ref (Ns     1 s) | L.null (L.tail ref) = Short 'N'          (clear s)
patch ref (Eqs    1 s) | L.null (L.tail ref) = Short (L.head ref) (clear s)
patch ref (Eqs1   1 s) | L.null (L.tail ref) = Short (L.head ref) (clear s)
patch ref (Chrs x _ s) | L.null (L.tail ref) = Short (tr x)       (clear s)

patch ref (Chrs x y s) = Short (tr x) $ Short (tr y) $ patch (L.drop  2    ref) s
patch ref (Eqs    n s) = Long     (L.take (2*n) ref) $ patch (L.drop (2*n) ref) s
patch ref (Eqs1   n s) = Long     (L.take (2*n) ref) $ patch (L.drop (2*n) ref) s
patch ref (Ns     n s) = Long   (L.replicate n' 'N') $ patch (L.drop (2*n) ref) s
  where n' = L.length $ L.take (2*n) ref

clear :: Stretch -> Frag
clear (Ns     _ s) = clear s
clear (Eqs    _ s) = clear s
clear (Eqs1   _ s) = clear s
clear (Chrs _ _ s) = clear s
clear (Break    s) = Term  s
clear Done         = Term Done

