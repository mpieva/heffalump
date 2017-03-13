module Lump where

import Data.Fix

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
-- length and a number of bases.
-- 
-- Needed Codes:
--   (Complex indels have to be broken into multiple variants.)
-- 

data Lump a
    = NNs !Int a                        -- ^ uncalled stretch

    | Eqs1 !Int a                       -- ^ stretch of monoallelic matches
    | Trans1 a                          -- ^ monoallelic transversion
    | Compl1 a                          -- ^ monoallelic complement
    | TCompl1 a                         -- ^ monoallelic trans-complement

    | Eqs2 !Int a                       -- ^ stretch of diallelic matches
    | RefTrans a                        -- ^ nine more diallelic SNP codes (*groan)
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

encodeLump :: Fix Lump -> Builder
encodeLump = cata encode1
  where
    encode1 (NNs n b)                          = stretch_of 0x40 n <> b
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
        | n < 0x100                            = word8 (k .|. 0x3F) <> word8 (fromIntegral n)
        | n < 0x10000                          = word8 (k .|. 0x3E) <> word8 (fromIntegral (n .&. 0xff)) 
                                                                    <> word8 (fromIntegral (n `shiftR`  8 .&. 0xff)
        | n < 0x1000000                        = word8 (k .|. 0x3D) <> word8 (fromIntegral (n .&. 0xff))
                                                                    <> word8 (fromIntegral (n `shiftR`  8 .&. 0xff)
                                                                    <> word8 (fromIntegral (n `shiftR` 16 .&. 0xff)
        | otherwise                            = word8 (k .|. 0x3C) <> word8 (fromIntegral (n .&. 0xff))
                                                                    <> word8 (fromIntegral (n `shiftR`  8 .&. 0xff)
                                                                    <> word8 (fromIntegral (n `shiftR` 16 .&. 0xff)
                                                                    <> word8 (fromIntegral (n `shiftR` 24 .&. 0xff)
    indel_of k n 
        | n < 8         = word8 (k .|. fromIntegral n)
        | n < 0x100     = word8 k <> word8 (fromIntegral n)

    seq_of s 
        | U.length s == 0 = mempty
        | U.length s == 1 = word8 (U.unsafeIndex s 0)
        | U.length s == 2 = word8 (U.unsafeIndex s 0 .|.
                                  (U.unsafeIndex s 1 `shiftL` 2)
        | U.length s == 3 = word8 (U.unsafeIndex s 0 .|.
                                  (U.unsafeIndex s 1 `shiftL` 2 .|.
                                  (U.unsafeIndex s 2 `shiftL` 4)
        | U.length s == 4 = word8 (U.unsafeIndex s 0 .|.
                                  (U.unsafeIndex s 1 `shiftL` 2 .|.
                                  (U.unsafeIndex s 2 `shiftL` 4 .|.
                                  (U.unsafeIndex s 3 `shiftL` 6)
                            <> seq_of (U.unsafeDrop 4 s)

--     if y == 1, it's an insertion and xxx bases follow (2 bits per base,
--                ATGC, padded to full bytes)

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
-- We are left with four "reserved" codes (0x04, 0x0D, 0x0E, 0x0F) and three
-- nonsensical ones (0x40, 0x80, 0xC0).  Should the integers maybe be coded
-- similarly to UTF8?  Loses some codes for short stretches, increases the
-- range for the longer codes.  Complicates the decoder, I think.
-- 
-- Note the bases: bit 0 codes for transition, bit 1 for complement.
-- Possible codes when outputting SNPs without knowledge of the reference
-- base:  I (Identical), O (transitiOn), P (comPlement), X
-- (trans-complement)
