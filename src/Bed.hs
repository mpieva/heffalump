module Bed where

import Bio.Prelude
import Streaming

import qualified Data.ByteString                as B
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.Vector.Algorithms.Intro   as U ( sort )
import qualified Data.Vector.Unboxed            as U
import qualified Streaming.Prelude              as Q

import NewRef ( Variant(..) )

newtype Bed = Bed ( U.Vector (Int32,Int32,Int32) )

mkBedFilter :: Monad m
            => Maybe FilePath -> [ B.ByteString ]
            -> IO ( Stream (Of Variant) m r -> Stream (Of Variant) m r )
mkBedFilter Nothing _ = return id
mkBedFilter (Just fp) refs = do
        bed <- L.readFile fp >>= parseBed refs
        return $ filterWithBed bed

parseBed :: [ B.ByteString ] -> L.ByteString -> IO Bed
parseBed chrms raw = do
    vau <- U.unsafeThaw . U.fromList . mapMaybe parse1 . L.lines $ raw
    U.sort vau
    Bed <$> U.unsafeFreeze vau
  where
    parse1 :: L.ByteString -> Maybe (Int32, Int32, Int32)
    parse1 ln = do sq:frm:tho:_ <- Just $ L.words ln
                   ci <- findIndex ((==) sq . L.fromStrict) chrms
                   (start,"") <- L.readInt frm
                   (end,"") <- L.readInt tho
                   return (fromIntegral ci, fromIntegral start, fromIntegral end)

filterWithBed :: Monad m => Bed -> Stream (Of Variant) m r -> Stream (Of Variant) m r
filterWithBed (Bed vau) = go (U.toList vau)
  where
    go [                 ] = lift . Q.effects >=> pure
    go ((ch, ps, pe) : rs) = lift . Q.next >=> \case
        Left    r    -> pure r
        Right (v,vs)
            -- variant before interval, drop the variant
            | (v_chr v, v_pos v) <  (fromIntegral ch, fromIntegral ps) -> go ((ch,ps,pe):rs) vs

            -- variant after interval, drop the interval
            | (v_chr v, v_pos v) >= (fromIntegral ch, fromIntegral pe) -> go rs (Q.cons v vs)

            -- must be a good variant
            | otherwise                                                -> v `Q.cons` go ((ch,ps,pe):rs) vs

