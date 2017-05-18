module Bed where

import BasePrelude

import qualified Data.ByteString                as B
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Data.Vector.Algorithms.Intro   as U ( sort )
import qualified Data.Vector.Unboxed            as U

import NewRef ( Variant(..) )

newtype Bed = Bed ( U.Vector (Int32,Int32,Int32) )

mkBedFilter :: Maybe FilePath -> [ B.ByteString ] -> IO ( [Variant] -> [Variant] )
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

filterWithBed :: Bed -> [ Variant ] -> [ Variant ]
filterWithBed (Bed vau) = go (U.toList vau)
  where
    go [                 ] _ = []
    go                  _ [] = []
    go ((ch, ps, pe) : rs) (v:vs)

        -- variant before interval, drop the variant
        | (v_chr v, v_pos v) < (fromIntegral ch, fromIntegral ps)  = go ((ch,ps,pe):rs) vs

        -- variant after interval, drop the interval
        | (v_chr v, v_pos v) >= (fromIntegral ch, fromIntegral pe) = go rs (v:vs)

        -- must be a good variant
        | otherwise                                                = v : go ((ch,ps,pe):rs) vs

