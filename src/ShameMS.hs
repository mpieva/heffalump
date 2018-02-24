module ShameMS ( mainShameMSout ) where

import Bio.Prelude
import System.Console.GetOpt
import Streaming

import qualified Data.ByteString.Builder         as B
import qualified Data.Vector                     as W
import qualified Data.Vector.Generic             as V
import qualified Data.Vector.Generic.Mutable     as M
import qualified Data.Vector.Storable            as U
import qualified Streaming.Prelude               as Q

import Bed ( mkBedFilter)
import Genome
import Lump
import Util ( parseFile1Opts )

data ConfShameMS = ConfShameMS {
    conf_noutgroups      :: Maybe Int,
    conf_blocklength     :: Maybe Int,
    conf_min_informative :: Double,
    conf_all             :: Bool,
    conf_split           :: Bool,
    conf_reference       :: Maybe FilePath,
    conf_regions         :: Maybe FilePath }
  deriving Show

defaultMsConf :: ConfShameMS
defaultMsConf = ConfShameMS (Just 0) Nothing 90.0 True True Nothing Nothing

optsVcfout :: [ OptDescr (ConfShameMS -> IO ConfShameMS) ]
optsVcfout =
    [ Option "r" ["reference"]          (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "n" ["numoutgroups"]       (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (0)"
    , Option "l" ["blocklength"]      (ReqArg set_lblock "NUM") "The length of each block (0)"
    , Option "m" ["min-informative"] (ReqArg set_tblock "PERC") "The minimum percentage of informative sites per block (90)"
    , Option "D" ["dense"]                    (NoArg set_dense) "Output invariant sites, too"
    , Option "t" ["only-transversions"]      (NoArg set_no_all) "Output only transversion sites"
    , Option "b" ["only-biallelic"]        (NoArg set_no_split) "Discard, don't split, polyallelic sites"
    , Option "R" ["regions"]           (ReqArg set_rgns "FILE") "Restrict to regions in bed-file FILE" ]
  where
    set_ref    a c = return $ c { conf_reference  =  Just a }
    set_nout   a c = (\n ->   c { conf_noutgroups =  Just n }) <$> readIO a
    set_lblock a c = (\n ->   c { conf_blocklength = Just n }) <$> readIO a
    set_dense    c = return $ c { conf_noutgroups = Nothing }
    set_no_all   c = return $ c { conf_all        =   False }
    set_no_split c = return $ c { conf_split      =   False }
    set_rgns   a c = return $ c { conf_regions    =  Just a }

    set_tblock a c = do n <- readIO a
                        if n < 0 || n > 100
                          then error "percentage must be between 0 and 100"
                          else return $ c { conf_min_informative = n / 100 }


mainShameMSout :: [String] -> IO ()
mainShameMSout args = do
    ( hefs, ConfShameMS{..} ) <- parseFile1Opts defaultMsConf "pseudo_ms"
                                                "[hef-file...]" optsVcfout args

    decodeMany conf_reference hefs $ \refs inps -> do
        region_filter <- mkBedFilter conf_regions (either error rss_chroms refs)

        mapsM_ (blocktoShame conf_blocklength conf_min_informative (V.length inps)) $
            maybe yields blocks conf_blocklength $
            region_filter $
            bool singles_only Q.concat conf_split $
            addRef (either error id refs) $
            maybe mergeLumpsDense mergeLumps conf_noutgroups inps
  where
    singles_only = Q.mapMaybe (\case [x] -> Just x ; _ -> Nothing)

    blocktoShame :: Maybe Int -> Double -> Int -> Stream (Of Variant) IO r -> IO r
    blocktoShame ml t m s = do
        ff :> r <- Q.mapOf (Block m) <$> vconcats (Q.map smashVariants s)
        unless (nullBlock ff) $ do
            forM_ ml $ \l -> hPutStrLn stdout $ "\n//\nblockSize_" ++ show l

            if enoughInfo t ff
                then B.hPutBuilder stdout $ foldMap
                            (\ind -> oneLine fstW ind <> oneLine sndW ind)
                            (individuals ff)
                else hPutStrLn stdout $ "# Not enough of info in this block"
        return r

    oneLine :: (Word8 -> Nuc2b) -> U.Vector Word8 -> B.Builder
    oneLine which = U.foldr ((<>) . B.char7 . toRefCode . which) (B.char7 '\n')

    enoughInfo :: Double -> Block -> Bool
    enoughInfo t blk = fromIntegral n_valid >= t/100 * fromIntegral n_variants
      where
        n_variants = W.length $ variants blk
        n_valid    = W.length $ W.filter (U.all $ \w -> isKnown (fstW w) && isKnown (sndW w)) $ variants blk

    smashVariants :: Variant -> U.Vector Word8
    smashVariants Variant{..} =
        U.map (smashVariant v_ref v_alt) v_calls

    smashVariant :: Nuc2b -> Var2b -> AlleleCounts -> Word8
    smashVariant (N2b r) (V2b a) = \case
        AC 0 0 ->      15 `pack2` 15
        AC _ 0 ->       r `pack2` r
        AC 0 _ -> xor r a `pack2` xor r a
        AC _ _ ->       r `pack2` xor r a

    pack2 :: Word8 -> Word8 -> Word8
    pack2 a b = a .|. shiftL b 4

    fstW, sndW :: Word8 -> Nuc2b
    fstW x = N2b $        x   .&. 0xF
    sndW x = N2b $ shiftR x 4 .&. 0xF


blocks :: Monad m => Int -> Stream (Of Variant) m r -> Stream (Stream (Of Variant) m) m r
blocks ln = go (-1) 0
  where
    go c p = lift . Q.next >=> \case
        Left    r                   -> pure r
        Right (v,vs) | v_chr v == c -> yields (Q.span (before     c (p+ln)) (Q.cons v vs)) >>= go     c (p+ln)
                     | otherwise    -> yields (Q.span (before (v_chr v) ln) (Q.cons v vs)) >>= go (v_chr v) ln

    before c p v = v_chr v == c && v_pos v < p


-- | Our variant data comes in as a stream of variants, but for output
-- we need a stream of individuals.  Put another way, the rows of our
-- data matrix are variants, and we need to transpose it so the variants
-- are in the columns.  This means we need to buffer at least a block in
-- memory.
--
-- To keep storage compact, we encode two alleles into a single 'Word8',
-- then store the whole matrix in a 'U.Vector' 'Word8'.  Use 'variants'
-- to access the 'Block' as a list of variants (traverse row-by-row) or
-- 'individuals' to access it as a list of individuals (traverse
-- column-by-column).

data Block = Block !Int                 -- ^ the stride (number of individuals)
                   !(U.Vector Word8)    -- ^ the data matrix

nullBlock :: Block -> Bool
nullBlock (Block _ v) = U.null v

variants :: Block -> W.Vector (U.Vector Word8)
variants (Block w ff) =
    W.map (\i -> U.slice i w ff) $
    W.enumFromStepN 0 w h
  where
    h = U.length ff `div` w

individuals :: Block -> W.Vector (U.Vector Word8)
individuals (Block w ff) =
    W.map (\j -> U.map (ff U.!) $ U.enumFromStepN j w h) $
    W.enumFromN 0 w
  where
    h = U.length ff `div` w

vconcats :: (V.Vector v a, MonadIO m) => Stream (Of (v a)) m r -> m (Of (v a) r)
vconcats s = do v <- liftIO $ M.unsafeNew 1024 ; go v 0 s
  where
    go !v !i = Q.next >=> \case
        Left    r    -> liftIO $ (:> r) <$> V.unsafeFreeze (M.take i v)
        Right (x,xs) -> do
            let lx = V.length x
            v' <- if i + lx <= M.length v
                    then return v
                    else liftIO $ M.grow v (M.length v `max` lx)
            liftIO $ M.move (M.slice i lx v') =<< V.unsafeThaw x
            go v' (i + lx) xs
