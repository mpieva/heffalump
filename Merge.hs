module Merge where

import BasePrelude

import qualified Data.ByteString.Char8          as B
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.ByteString.Internal       as B
import qualified Data.ByteString.Unsafe         as B
import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U

import Stretch
import Util ( Reference(..) )

{-# INLINE tr' #-}
tr' :: Alleles -> Char
tr' (Alleles w) = B.w2c . B.unsafeIndex chars_for_alleles . fromIntegral $ w

chars_for_alleles :: B.ByteString
chars_for_alleles = "NACMGRSVTWYHKDBN"

newtype Alleles = Alleles Word8 deriving (Bits, Eq, Show)

{-# INLINE alleles #-}
alleles :: NucCode -> Alleles
alleles (NucCode w) = Alleles . B.unsafeIndex vector_alleles . fromIntegral $ w

vector_alleles :: B.ByteString
vector_alleles = "\0\1\2\4\8\3\5\9\6\10\12\1\2\4\8\0"

-- A variant call.  Has a genomic position, ref and alt alleles, and a
-- bunch of calls.
data Variant = Variant { v_chr :: !Int                  -- index into 'chroms'
                       , v_pos :: !Int                  -- 0-based
                       , v_ref :: !Char                 -- A,C,G or T
                       , v_alt :: !Char                 -- A,C,G or T
                       , v_calls :: !(U.Vector Word8) }
  deriving Show

merge_hefs :: Bool -> Int -> Reference -> [Stretch] -> [Variant]
merge_hefs nosplit nougroups ref = concat . merge_hefs' nosplit nougroups ref

merge_hefs' :: Bool -> Int -> Reference -> [Stretch] -> [[Variant]]
merge_hefs'       !_          !_ (Reference [          ]) = const [ ]
merge_hefs' !nosplit !noutgroups (Reference (ref0:refs0)) = go refs0 0 0 ref0 . V.fromList
    -- Merging stretches.  We define a 'Variant' as anything that is
    -- different from the reference.  Therefore, 'Eqs' ('Eqs1') and 'Ns'
    -- never create a 'Variant' and we can skip forwards.  A 'Done' is
    -- an error.  Else, if at least one sample has 'Chrs', we need to
    -- collect the alleles, check if we found exactly two, then output
    -- and skip forward by two bases.
  where
    go :: [L.ByteString] -> Int -> Int -> L.ByteString -> V.Vector Stretch -> [[Variant]]
    go refs !ix !pos !ref !smps
        | L.null ref = case refs of (r:rs) -> go rs (succ ix) 0 r (V.map skipBreaks smps)
                                    [    ] -> []
        | L.null (L.tail ref) =
                 tryVar  ix  pos    (code $ L.index ref 0)  firstVars smps $
                 case refs of (r:rs) -> go rs (succ ix) 0 r (V.map skipBreaks smps)
                              [    ] -> []

        -- We ignore outgroups in computing the longest stretch.  That
        -- way, outgroups don't get to define variants, put participate
        -- in variants found in another way.
        | otherwise  = case longestStretch $ V.drop noutgroups smps of
            -- no stretch, have to look for vars
            0 -> tryVar  ix  pos    (code_ref $ L.index ref 0)  firstVars smps $
                 tryVar  ix (pos+1) (code_ref $ L.index ref 1) secondVars smps $
                 go refs ix (pos+2) (L.drop 2 ref) (V.map (skipStretch 1) smps)
            -- maxBound == stretches to infinity ==> we a re done with this ref
            l | l == maxBound ->
                 case refs of
                        [    ] -> []
                        (r:rs) -> go rs (succ ix) 0 r (V.map (skipBreaks . skipStretch lref) smps)
                            where lref = succ (L.length ref) `div` 2

            -- a stretch of l non-variants can be skipped over
            l -> go refs ix (pos + 2 * fromIntegral l) (L.drop (2*l) ref) (unS $ V.mapM (S . skipStretch l) smps)

    code_ref 'a' = NucCode 11
    code_ref 'c' = NucCode 12
    code_ref 'g' = NucCode 13
    code_ref 't' = NucCode 14
    code_ref 'A' = NucCode 11
    code_ref 'C' = NucCode 12
    code_ref 'G' = NucCode 13
    code_ref 'T' = NucCode 14
    code_ref  _  = NucCode 0

    -- Skip over 2*l(!) sites.
    skipStretch :: Int64 -> Stretch -> Stretch
    skipStretch !l           _ | l <= 0 || l > 300000000 = error "WTF?!"
    skipStretch !l (Ns   !n s) | l == n = s
                               | l <  n = Ns (n-l) s
                               | otherwise = skipStretch (l-n) s
    skipStretch !l (Eqs  !n s) | l == n = s
                               | l <  n = Eqs (n-l) s
                               | otherwise = skipStretch (l-n) s
    skipStretch !l (Eqs1 !n s) | l == n = s
                               | l <  n = Eqs1 (n-l) s
                               | otherwise = skipStretch (l-n) s
    skipStretch  1 (Chrs _ _ s) = s
    skipStretch !n (Chrs _ _ s) = skipStretch (n-1) s

    skipStretch _ (Break    s) = Break s
    skipStretch _  Done        = Done

    skipBreaks (Break s) = s
    skipBreaks        s  = s

    longestStretch :: V.Vector Stretch -> Int64
    longestStretch = V.minimum . V.map (\s -> case s of Eqs n _ -> n ; Eqs1 n _ -> n ; Ns n _ -> n ; Chrs _ _ _ -> 0 ; _ -> maxBound)

    firstVars, secondVars :: NucCode -> V.Vector Stretch -> V.Vector NucCode
    firstVars  r = V.map $ \s -> case s of Eqs _ _ -> twice r ; Eqs1 _ _ -> once r ; Chrs x _ _ -> x ; _ -> NucCode 0
    secondVars r = V.map $ \s -> case s of Eqs _ _ -> twice r ; Eqs1 _ _ -> once r ; Chrs _ y _ -> y ; _ -> NucCode 0

    twice (NucCode x) |    x == 0 = NucCode 0
                      |    x > 10 = NucCode (x-10)
                      | otherwise = NucCode x

    once (NucCode x) |    x == 0 = NucCode 0
                     |    x <= 4 = NucCode (x+10)
                     | otherwise = NucCode x

    -- Find all the variants, anchored on the reference allele, and
    -- split them.  Misfitting alleles are not counted.
    tryVar :: Int -> Int -> NucCode -> (NucCode -> V.Vector Stretch -> V.Vector NucCode) -> V.Vector Stretch -> [[Variant]] -> [[Variant]]
    tryVar ix pos r which ss
        | r == NucCode 0      = id
        | nosplit && not good = id
        | otherwise           = (:) [ Variant ix pos (tr r) (tr' v) (V.convert $ V.map (ct (alleles r) v) vs)
                                    | v <- map Alleles [1,2,4,8], vacc .&. v /= Alleles 0 ]
      where
        vs   = which r ss
        -- collect variant alleles, ref and outgroups don't count
        vacc = V.foldl' (\a c -> a .|. alleles c) (Alleles 0) (V.drop noutgroups vs) .&. complement (alleles r)
        -- good == exactly one variant
        good = vacc `elem` map Alleles [1,2,4,8]

-- Variant codes:  #ref + 4 * #alt
ct :: Alleles -> Alleles -> NucCode -> Word8
ct r v n | alleles n == r                       = if isHap n then 1 else 2
         | alleles n == v                       = if isHap n then 4 else 8
         | alleles n == r .|. v                 = 5
         | otherwise                            = 0

isHap :: NucCode -> Bool
isHap (NucCode n) = n >= 11

-- | This gunk is need to make a map over a 'Vector' strict.  Looks
-- ugly, but is perfectly servicable.
data S a = S { unS :: !a }

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
    S a >>= k = k a

