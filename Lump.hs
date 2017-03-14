{-# LANGUAGE OverloadedLists #-}
module Lump where

import BasePrelude
import Data.Fix
import Data.ByteString.Builder

import qualified Data.ByteString.Lazy as L
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

-- | Improved diff representation and encoding.  Goals:
--
-- * support one or two known alleles
-- * do not force SNPs into pairs
-- * allow short indels
-- * allow many operations /without access to the reference/
--
-- Outline:
--
-- We store SNP bases as patch to the reference.  Instead of a base, we
-- store the change needed:  none, transition, complement,
-- trans-complement.  Deletions just have a length, insertions have a
-- length and a sequence of bases.  Complex alleles (such as two different
-- alleles) would need to be broken up into separate variants---haven't
-- thought about the details.   XXX

data Lump a
    = Ns !Int a                         -- ^ uncalled stretch

    | Eqs1 !Int a                       -- ^ stretch of monoallelic matches
    | Trans1 a                          -- ^ monoallelic transversion
    | Compl1 a                          -- ^ monoallelic complement
    | TCompl1 a                         -- ^ monoallelic trans-complement

    | Eqs2 !Int a                       -- ^ stretch of diallelic matches
    | RefTrans a                        -- ^ nine more diallelic SNP codes (*groan*)
    | Trans2 a
    | RefCompl a
    | TransCompl a
    | Compl2 a
    | RefTCompl a
    | TransTCompl a
    | ComplTCompl a
    | TCompl2 a

    | Del1 !Int a                       -- ^ deletion with other allele unknown
    | Del2 !Int a                       -- ^ homozygous deletion
    | DelH !Int a                       -- ^ deletion, the other allele is the reference

    | Ins1 !Seq1 a                      -- ^ insertion with other allele unknown
    | Ins2 !Seq1 a                      -- ^ homozygous insertion
    | InsH !Seq1 a                      -- ^ insertion, the other allele is the reference

    | Break a                           -- ^ break marker (end-of-chromosome)
    | Done                              -- ^ end of data stream
  deriving Functor

type Seq1 = U.Vector Word8

encodeLump :: Fix Lump -> Builder
encodeLump = cata encode1
  where
    encode1 (Ns n b)                           = stretch_of 0x40 n <> b
    encode1 (Eqs1 n b)                         = stretch_of 0x80 n <> b
    encode1 (Trans1 b)                         = word8 0x01 <> b
    encode1 (Compl1 b)                         = word8 0x02 <> b
    encode1 (TCompl1 b)                        = word8 0x03 <> b

    encode1 (Eqs2 n b)                         = stretch_of 0xC0 n <> b
    encode1 (RefTrans b)                       = word8 0x05 <> b
    encode1 (Trans2 b)                         = word8 0x06 <> b
    encode1 (RefCompl b)                       = word8 0x07 <> b
    encode1 (TransCompl b)                     = word8 0x08 <> b
    encode1 (Compl2 b)                         = word8 0x09 <> b
    encode1 (RefTCompl b)                      = word8 0x0A <> b
    encode1 (TransTCompl b)                    = word8 0x0B <> b
    encode1 (ComplTCompl b)                    = word8 0x0C <> b
    encode1 (TCompl2 b)                        = word8 0x0D <> b

    encode1 (Del1 n b)                         = indel_of 0x10 n <> b
    encode1 (DelH n b)                         = indel_of 0x20 n <> b
    encode1 (Del2 n b)                         = indel_of 0x30 n <> b

    encode1 (Ins1 s b)                         = indel_of 0x18 (U.length s) <> seq_of s <> b
    encode1 (InsH s b)                         = indel_of 0x28 (U.length s) <> seq_of s <> b
    encode1 (Ins2 s b)                         = indel_of 0x38 (U.length s) <> seq_of s <> b

    encode1 (Break b)                          = word8 0x00 <> b
    encode1  Done                              = mempty

    stretch_of k n
        | n < 0x3C                             = word8 (k .|. fromIntegral n)
        | n < 0x100                            = word8 (k .|. 0x3F) <> word8 (fromIntegral  n)
        | n < 0x10000                          = word8 (k .|. 0x3E) <> word8 (fromIntegral (n             .&. 0xff))
                                                                    <> word8 (fromIntegral (n `shiftR`  8 .&. 0xff))
        | n < 0x1000000                        = word8 (k .|. 0x3D) <> word8 (fromIntegral (n             .&. 0xff))
                                                                    <> word8 (fromIntegral (n `shiftR`  8 .&. 0xff))
                                                                    <> word8 (fromIntegral (n `shiftR` 16 .&. 0xff))
        | otherwise                            = word8 (k .|. 0x3C) <> word8 (fromIntegral (n             .&. 0xff))
                                                                    <> word8 (fromIntegral (n `shiftR`  8 .&. 0xff))
                                                                    <> word8 (fromIntegral (n `shiftR` 16 .&. 0xff))
                                                                    <> word8 (fromIntegral (n `shiftR` 24 .&. 0xff))
    indel_of k n
        | n == 0        = error "empty indel"
        | n < 8         = word8 (k .|. fromIntegral n)
        | n < 0x100     = word8 k <> word8 (fromIntegral n)
        | otherwise     = error $ "long indel: " ++ show n

    seq_of s
        | U.length s == 0 = mempty
        | U.length s == 1 = word8 (U.unsafeIndex s 0)
        | U.length s == 2 = word8 (U.unsafeIndex s 0 .|.
                                   U.unsafeIndex s 1 `shiftL` 2)
        | U.length s == 3 = word8 (U.unsafeIndex s 0 .|.
                                   U.unsafeIndex s 1 `shiftL` 2 .|.
                                   U.unsafeIndex s 2 `shiftL` 4)
        | otherwise       = word8 (U.unsafeIndex s 0 .|.
                                   U.unsafeIndex s 1 `shiftL` 2 .|.
                                   U.unsafeIndex s 2 `shiftL` 4 .|.
                                   U.unsafeIndex s 3 `shiftL` 6)
                            <> seq_of (U.unsafeDrop 4 s)


decodeLump :: L.ByteString -> Fix Lump
decodeLump = ana decode1
  where
    decode1 s = case L.uncons s of
        Nothing                                -> Done
        Just (c,s1)
            | c == 0x00                        -> Break s1
            | c == 0x01                        -> Trans1 s1
            | c == 0x02                        -> Compl1 s1
            | c == 0x03                        -> TCompl1 s1
            | c == 0x04                        -> error $ "unexpected " ++ show c
            | c == 0x05                        -> RefTrans s1
            | c == 0x06                        -> Trans2 s1
            | c == 0x07                        -> RefCompl s1
            | c == 0x08                        -> TransCompl s1
            | c == 0x09                        -> Compl2 s1
            | c == 0x0A                        -> RefTCompl s1
            | c == 0x0B                        -> TransTCompl s1
            | c == 0x0C                        -> ComplTCompl s1
            | c == 0x0D                        -> TCompl2 s1
            | c == 0x0E                        -> error $ "unexpected " ++ show c
            | c == 0x0F                        -> error $ "unexpected " ++ show c

            | c .&. 0xF8 == 0x10               -> del_of Del1 c s1
            | c .&. 0xF8 == 0x20               -> del_of DelH c s1
            | c .&. 0xF8 == 0x30               -> del_of Del2 c s1

            | c .&. 0xF8 == 0x18               -> ins_of Ins1 c s1
            | c .&. 0xF8 == 0x28               -> ins_of InsH c s1
            | c .&. 0xF8 == 0x38               -> ins_of Ins2 c s1

            | c .&. 0xC0 == 0x40               -> stretch_of Ns c s1
            | c .&. 0xC0 == 0x80               -> stretch_of Eqs1 c s1
            | c .&. 0xC0 == 0xC0               -> stretch_of Eqs2 c s1

            | otherwise                        -> error $ "Impossibru! " ++ show c

    stretch_of cons c s1
        | c .&. 0x40 == 0x3F  = cons (fromIntegral (L.index s1 0)) (L.drop 1 s1)
        | c .&. 0x40 == 0x3E  = cons (fromIntegral (L.index s1 0) .|.
                                      fromIntegral (L.index s1 1) `shiftL`  8) (L.drop 2 s1)
        | c .&. 0x40 == 0x3D  = cons (fromIntegral (L.index s1 0) .|.
                                      fromIntegral (L.index s1 1) `shiftL`  8 .|.
                                      fromIntegral (L.index s1 2) `shiftL` 16) (L.drop 3 s1)
        | c .&. 0x40 == 0x3C  = cons (fromIntegral (L.index s1 0) .|.
                                      fromIntegral (L.index s1 1) `shiftL`  8 .|.
                                      fromIntegral (L.index s1 1) `shiftL` 16 .|.
                                      fromIntegral (L.index s1 1) `shiftL` 24) (L.drop 4 s1)
        | otherwise           = cons (fromIntegral (c .&. 0x40)) s1

    del_of cons c s1
        | c .&. 0x07 == 0  = cons (fromIntegral (L.head s1)) (L.tail s1)
        | otherwise        = cons (fromIntegral c .&. 0x07) s1

    ins_of cons c s1
        | c .&. 0x07 == 0  = seq_of cons U.empty (L.head s1) (L.tail s1)
        | otherwise        = seq_of cons U.empty (c .&. 0x07) s1



    seq_of cons acc 0 s1 = cons acc s1

    seq_of cons acc 1 s1 = cons acc' (L.tail s1)
      where
        acc' = acc `U.snoc` (L.head s1 `shiftR` 0 .&. 0x3)

    seq_of cons acc 2 s1 = cons acc' (L.tail s1)
      where
        acc' = acc U.++ [ L.head s1 `shiftR` 0 .&. 0x3
                        , L.head s1 `shiftR` 2 .&. 0x3 ]

    seq_of cons acc 3 s1 = cons acc' (L.tail s1)
      where
        acc' = acc U.++ [ L.head s1 `shiftR` 0 .&. 0x3
                        , L.head s1 `shiftR` 2 .&. 0x3
                        , L.head s1 `shiftR` 4 .&. 0x3 ]

    seq_of cons acc n s1 = seq_of cons acc' (n-4) (L.tail s1)
      where
        acc' = acc U.++ [ L.head s1 `shiftR` 0 .&. 0x3
                        , L.head s1 `shiftR` 2 .&. 0x3
                        , L.head s1 `shiftR` 4 .&. 0x3
                        , L.head s1 `shiftR` 6 .&. 0x3 ]


-- 0b00000000: break
--
-- 0b00000001: transition
-- 0b00000010: complement
-- 0b00000011: trans-complement
--
-- 0b00000101: ref+trans
-- 0b00000110: trans+trans
-- 0b00000111: ref+compl
-- 0b00001000: trans+compl
-- 0b00001001: compl+compl
-- 0b00001010: ref+tcompl
-- 0b00001011: trans+tcompl
-- 0b00001100: compl+tcompl
-- 0b00001101: tcompl+tcompl
--
-- 0b0001yxxx: monoallelic indel
-- 0b0010yxxx: heterozygous indel
-- 0b0011yxxx: diallelic indel
--     if xxx == 0, one length byte follows.
--     if y == 1, it's an insertion and xxx bases follow (2 bits per base,
--                ATGC, padded to full bytes)
--
-- 0b01xxxxxx: short uncalled stretch
-- 0b10xxxxxx: short monoallelic matches
-- 0b11xxxxxx: short diallelic matches
--     if xxxxxx == 111111, one length byte follows
--     if xxxxxx == 111110, two length bytes follow
--     if xxxxxx == 111101, three length bytes follow
--     if xxxxxx == 111100, four length bytes follow
--
--
-- We are left with three "reserved" codes (0x04, 0x0E, 0x0F) and three
-- nonsensical ones (0x40, 0x80, 0xC0).  Should the integers maybe be coded
-- similarly to UTF8?  Loses some codes for short stretches, increases the
-- range for the longer codes.  Complicates the decoder, I think.
--
-- Note the bases: bit 0 codes for transition, bit 1 for complement.
-- Possible codes when outputting SNPs without knowledge of the reference
-- base:  I (Identical), O (transitiOn), P (comPlement), X
-- (trans-complement)


-- A variant call.  Has a genomic position, ref and alt alleles, and a
-- bunch of calls.
data Variant = Variant { v_chr :: !Int                  -- chromosome number
                       , v_pos :: !Int                  -- 0-based
                       , v_alt :: !Word8                -- [0..3] for "IOPX"
                       , v_calls :: !(U.Vector Word8) } -- XXX
  deriving Show

-- | Merging without access to a reference sequence.  This code doesn't
-- believe in Indels and skips over them.
--
-- 'nosplit':  don't split multiallelic variants.  XXX Necessary?
-- 'noutgroups':  number of outgroups

merge_lumps :: Int -> V.Vector (Fix Lump) -> [[Variant]]
merge_lumps !noutgroups = filter (not . null) . go 0 0
    -- Merging stretches.  We define a 'Variant' as anything that is
    -- different from the reference.  Therefore, 'Eqs' ('Eqs1') and 'Ns'
    -- never create a 'Variant' and we can skip forwards.  A 'Done' is
    -- an error.  Else, if at least one sample has some sort of call, we
    -- need to collect the alleles, produce up to four variants, then
    -- output and skip forward by one base.
  where
    go :: Int -> Int -> V.Vector (Fix Lump) -> [[Variant]]
    go !ix !pos !smps
        -- all samples 'Done', no more 'Variant's to produce
        | V.all (isDone . unFix) smps = []

        -- all samples 'Break' or are done, next chromosome
        | V.all (isBreak . unFix) smps = go (succ ix) 0 (V.map (skipBreak . unFix) smps)

        -- We ignore outgroups in computing the longest stretch.  That
        -- way, outgroups don't get to define variants, but participate
        -- in variants found in another way.
        | otherwise = case longestStretch $ V.drop noutgroups smps of

            -- no stretch, have to look for vars
            0 -> mkVar ix    pos                                            smps $
                 go ix (succ pos) (unS $ V.mapM (S . skipStretch 1 . unFix) smps)

            -- a stretch of l non-variants can be skipped over
            l -> go ix (pos + fromIntegral l) (unS $ V.mapM (S . skipStretch l . unFix) smps)


    isDone Done = True
    isDone    _ = False

    isBreak  Done     = True
    isBreak (Break _) = True
    isBreak        _  = False

    skipBreak (Break l) =     l
    skipBreak        l  = Fix l

    -- Skip over l sites.
    skipStretch :: Int -> Lump (Fix Lump) -> Fix Lump
    skipStretch !l           _ | l < 0 || l > 300000000 = error "WTF?!"

    skipStretch !l (Del1 _ s) = skipStretch l (unFix s)
    skipStretch !l (Del2 _ s) = skipStretch l (unFix s)
    skipStretch !l (DelH _ s) = skipStretch l (unFix s)

    skipStretch !l (Ins1 _ s) = skipStretch l (unFix s)
    skipStretch !l (Ins2 _ s) = skipStretch l (unFix s)
    skipStretch !l (InsH _ s) = skipStretch l (unFix s)

    skipStretch _ (Break    s) = Fix $ Break s
    skipStretch _  Done        = Fix $ Done

    skipStretch !l (Ns   !n s) | l == n = s
                               | l <  n = Fix $ Ns (n-l) s
                               | otherwise = skipStretch (l-n) (unFix s)
    skipStretch !l (Eqs1 !n s) | l == n = s
                               | l <  n = Fix $ Eqs1 (n-l) s
                               | otherwise = skipStretch (l-n) (unFix s)
    skipStretch !l (Eqs2 !n s) | l == n = s
                               | l <  n = Fix $ Eqs2 (n-l) s
                               | otherwise = skipStretch (l-n) (unFix s)

    skipStretch  0          s  = Fix s
    skipStretch !l (Trans1  s) = skipStretch (l-1) (unFix s)
    skipStretch !l (Compl1  s) = skipStretch (l-1) (unFix s)
    skipStretch !l (TCompl1 s) = skipStretch (l-1) (unFix s)

    skipStretch !l (RefTrans    s) = skipStretch (l-1) (unFix s)
    skipStretch !l (Trans2      s) = skipStretch (l-1) (unFix s)
    skipStretch !l (RefCompl    s) = skipStretch (l-1) (unFix s)
    skipStretch !l (TransCompl  s) = skipStretch (l-1) (unFix s)
    skipStretch !l (Compl2      s) = skipStretch (l-1) (unFix s)
    skipStretch !l (RefTCompl   s) = skipStretch (l-1) (unFix s)
    skipStretch !l (TransTCompl s) = skipStretch (l-1) (unFix s)
    skipStretch !l (ComplTCompl s) = skipStretch (l-1) (unFix s)
    skipStretch !l (TCompl2     s) = skipStretch (l-1) (unFix s)


    longestStretch :: V.Vector (Fix Lump) -> Int
    longestStretch = V.minimum . V.map (get_stretch_len . unFix)
      where
        get_stretch_len (Ns   n _) = n
        get_stretch_len (Eqs1 n _) = n
        get_stretch_len (Eqs2 n _) = n
        get_stretch_len (Break  _) = maxBound
        get_stretch_len  Done      = maxBound
        get_stretch_len  _         = 0


    -- Find all the variants, anchored on the reference allele, and
    -- split them.  Misfitting alleles are not counted.
    mkVar :: Int -> Int -> V.Vector (Fix Lump) -> [[Variant]] -> [[Variant]]
    mkVar ix pos ss = (:)
            [ Variant ix pos alt calls
            | (alt, ct) <- zip [1..3] [ct_trans, ct_compl, ct_tcompl]
            , let calls = V.convert $ V.map (ct . unFix) ss
            -- it's only a variant if at least one alt called
            , U.any (\c -> c .&. 0xC /= 0) (U.drop noutgroups calls) ]


    -- Variant codes:  #ref + 4 * #alt
    ct_trans :: Lump a -> Word8
    ct_trans (Eqs1      _ _) = 1
    ct_trans (Trans1      _) = 4

    ct_trans (Eqs2      _ _) = 2
    ct_trans (RefTrans    _) = 5
    ct_trans (Trans2      _) = 8
    ct_trans              _  = 0

    ct_compl :: Lump a -> Word8
    ct_compl (Eqs1      _ _) = 1
    ct_compl (Compl1      _) = 4

    ct_compl (Eqs2      _ _) = 2
    ct_compl (RefCompl    _) = 5
    ct_compl (Compl2      _) = 8
    ct_compl              _  = 0

    ct_tcompl :: Lump a -> Word8
    ct_tcompl (Eqs1      _ _) = 1
    ct_tcompl (TCompl1     _) = 4

    ct_tcompl (Eqs2      _ _) = 2
    ct_tcompl (RefTCompl   _) = 5
    ct_tcompl (TCompl2     _) = 8
    ct_tcompl              _  = 0


-- | This gunk is need to make a map over a 'Vector' strict.  Looks
-- ugly, but is perfectly servicable.
newtype S a = S { unS :: a }

instance Functor S where
    {-# INLINE fmap #-}
    fmap f (S a) = S (f a)

instance Applicative S where
    {-# INLINE pure #-}
    pure = S
    {-# INLINE (<*>) #-}
    S f <*> S a = S (f a)

instance Monad S where
    {-# INLINE return #-}
    return = S
    {-# INLINE (>>=) #-}
    S !a >>= k = k a

