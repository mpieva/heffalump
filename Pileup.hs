module Pileup where

-- Simple and fast Pileup... iff we can pull it off.

import Bio.Base
import Bio.Bam
import Control.Monad.Primitive

import qualified Data.Vector.Generic as V
import qualified Data.Vector.Generic.Mutable as VM
import qualified Data.Vector as VV
import qualified Data.Vector.Unboxed as VU

-- | An iterator over a BamRec.
-- (I'm inclined to make this more compact and copy compacted shit in
-- from the BamRec.  We'll see.)
data BamIter = BamIter
    { pos_seq :: Int           -- current position into the sequence
    , pos_cig :: Int           -- position into CIGAR
    , pos_cig_op :: Int        -- position into CIGAR op
    , bam_raw :: BamRaw }
  deriving Show

instance Show BamRaw where show = show . unpackBam

data VecIters a = VecIters
    { vec_seq :: !(VU.Vector Int)
    , vec_cig :: !(VU.Vector Int)
    , vec_cig_op :: !(VU.Vector Int)
    , vec_recs :: !(VV.Vector BamRaw) }

data MVecIters m a = MVecIters
    { mvec_seq :: !(VU.MVector m Int)
    , mvec_cig :: !(VU.MVector m Int)
    , mvec_cig_op :: !(VU.MVector m Int)
    , mvec_recs :: !(VV.MVector m BamRaw) }

type instance V.Mutable VecIters = MVecIters

instance V.Vector VecIters BamIter where
    basicUnsafeFreeze (MVecIters x y z w) =
        VecIters <$> V.basicUnsafeFreeze x <*> V.basicUnsafeFreeze y
                 <*> V.basicUnsafeFreeze z <*> V.basicUnsafeFreeze w

    basicUnsafeThaw (VecIters x y z w) =
        MVecIters <$> V.basicUnsafeThaw x <*> V.basicUnsafeThaw y
                  <*> V.basicUnsafeThaw z <*> V.basicUnsafeThaw w

    basicLength (VecIters x _ _ _) = V.basicLength x

    basicUnsafeSlice i l (VecIters x y z w) =
        VecIters (V.basicUnsafeSlice i l x) (V.basicUnsafeSlice i l y)
                 (V.basicUnsafeSlice i l z) (V.basicUnsafeSlice i l w)

    basicUnsafeIndexM (VecIters x y z w) i =
        BamIter <$> V.basicUnsafeIndexM x i <*> V.basicUnsafeIndexM y i
                <*> V.basicUnsafeIndexM z i <*> V.basicUnsafeIndexM w i

instance VM.MVector MVecIters BamIter where
    basicLength (MVecIters x _ _ _) = VM.basicLength x

    basicUnsafeSlice i l (MVecIters x y z w) =
        MVecIters (VM.basicUnsafeSlice i l x) (VM.basicUnsafeSlice i l y)
                  (VM.basicUnsafeSlice i l z) (VM.basicUnsafeSlice i l w)

    basicOverlaps (MVecIters a _ _ _) (MVecIters b _ _ _) = VM.basicOverlaps a b

    basicUnsafeNew l =
        MVecIters <$> VM.basicUnsafeNew l <*> VM.basicUnsafeNew l
                  <*> VM.basicUnsafeNew l <*> VM.basicUnsafeNew l

    basicUnsafeRead (MVecIters x y z w) i =
        BamIter <$> VM.basicUnsafeRead x i <*> VM.basicUnsafeRead y i
                <*> VM.basicUnsafeRead z i <*> VM.basicUnsafeRead w i

    basicUnsafeWrite (MVecIters x y z w) i BamIter{..} =
        VM.basicUnsafeWrite x i pos_seq >>
        VM.basicUnsafeWrite y i pos_cig >>
        VM.basicUnsafeWrite z i pos_cig_op >>
        VM.basicUnsafeWrite w i bam_raw

newBuf :: PrimMonad m => Int -> m (MVecIters (PrimState m) BamIter)
newBuf = VM.new

-- Makes an iterator from a BamRec.  Returns 'Nothing' if the CIGAR
-- doesn't contain anything valuable.
newBamIter :: BamRaw -> Maybe BamIter
newBamIter = scanBamIter . BamIter 0 0 0

-- | Steps an iterator to the next genomic location.  Returns 'Nothing'
-- instead of stepping off the end.
stepBamIter :: BamIter -> Maybe BamIter
stepBamIter bi
    = case b_cigar (unpackBam (bam_raw bi)) V.! pos_cig bi of
            Mat :* _ -> scanBamIter $ bi { pos_seq    = succ $ pos_seq    bi
                                         , pos_cig_op = succ $ pos_cig_op bi }
            Del :* _ -> scanBamIter $ bi { pos_cig_op = succ $ pos_cig_op bi }

-- | Steps an iterator until it validly points at a CIGAR op that would
-- advance along the reference.
scanBamIter :: BamIter -> Maybe BamIter
scanBamIter bi@BamIter{..}
    | pos_cig == V.length (b_cigar (unpackBam (bam_raw)))                   = Nothing

    | otherwise
        = case (b_cigar (unpackBam (bam_raw))) V.! pos_cig of
            Mat :* n | pos_cig_op /= n              -> Just bi
            Del :* n | pos_cig_op /= n              -> Just bi
            _   :* _ -> scanBamIter bi { pos_cig = succ pos_cig, pos_cig_op = 0 }

type SomeBase = (Nucleotides, Qual, Bool)

-- | Reads the current 'Nucleotides' base off the iterator.  Returns
-- 'Nothing' if we're looking at a deletion.
baseFromBamIter :: BamIter -> Maybe SomeBase
baseFromBamIter BamIter{..} | pos_cig >= V.length (b_cigar (unpackBam (bam_raw))) = error "Huh?"
baseFromBamIter BamIter{..} =
    case (b_cigar (unpackBam (bam_raw))) V.! pos_cig of
        Mat :* _ | pos_seq >= V.length (b_seq (unpackBam (bam_raw))) -> error $ "Haeh?" ++ show BamIter{..}
        Mat :* _ -> Just ( b_seq (unpackBam bam_raw)  V.! pos_seq
                         , b_qual (unpackBam bam_raw) V.! pos_seq
                         , isReversed (unpackBam bam_raw) )
        _   :* _ -> Nothing

-- Fold over the extracted stuff is implicit.
-- (We bail as soon as we hit an invalid refseq, we skip over unaligned
-- reads.)
pileup :: MonadIO m => (a -> SomeBase -> m a) -> (Refseq -> Int -> a) -> Enumeratee [BamRaw] [a] m b
pileup cons nil iter = do buf <- liftIO $ newBuf 256
                          eneeCheckIfDone (mSkip (Refseq 0) 0 0 buf) iter
  where
    mSkip !rs !po !nbuf !buf out
        | nbuf /= 0 = mFeed rs po nbuf buf out
        | otherwise -- nothing in buffer, skip to next read
            = do mbr <- peekStream
                 case fmap unpackBam mbr of
                    Nothing -> return (liftI out)
                    Just br
                        | not (isValidRefseq (b_rname br))
                            -> return (liftI out)

                        | isUnmapped br
                            -> mSkip rs po 0 buf out

                        | otherwise
                            -> mFeed (b_rname br) (b_pos br) 0 buf out


    mFeed !rs !po !nbuf !buf out = do
        mbr <- peekStream
        case mbr of
            Just br
                | isValidRefseq (b_rname (unpackBam br)) && isUnmapped (unpackBam br)
                    -> headStream >>
                       mFeed rs po nbuf buf out

                | b_rname (unpackBam br) == rs && b_pos (unpackBam br) == po
                    -> headStream >>
                       case newBamIter br of
                            Nothing -> mFeed rs po nbuf buf out
                            Just bi -> do
                                buf' <- if VM.length buf > nbuf then return buf
                                        else liftIO $ VM.grow buf (VM.length buf)
                                liftIO $ VM.write buf' nbuf bi
                                mFeed rs po (succ nbuf) buf' out

            _ -> mPile rs po nbuf buf out

    mPile !rs !po !nbuf !buf out = do
        !pile <- lift $ pile_loop nbuf buf 0 (nil rs po)
        !nbuf' <- liftIO $ step_loop nbuf buf 0 0
        eneeCheckIfDone (mSkip rs (succ po) nbuf' buf) . out $ Chunk [pile]

    pile_loop !nbuf !buf !i !acc
        | i == nbuf = return acc
        | otherwise = do bi <- liftIO $ VM.read buf i
                         case baseFromBamIter bi of
                            Nothing ->                pile_loop nbuf buf (succ i) acc
                            Just  b -> cons acc b >>= pile_loop nbuf buf (succ i)

    step_loop !nbuf !buf !i !j
        | i == nbuf = return j
        | otherwise = do bi <- VM.read buf i
                         case stepBamIter bi of
                            Nothing -> step_loop nbuf buf (succ i) j
                            Just bj -> do VM.write buf j bj
                                          step_loop nbuf buf (succ i) (succ j)


