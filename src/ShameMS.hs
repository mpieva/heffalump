module ShameMS ( mainShameMSout ) where

import Bio.Prelude
import System.Console.GetOpt
import Streaming


import Data.List.Split                           as LS
import qualified Data.Vector.Generic             as V
import qualified Data.Vector.Generic.Mutable     as M
import qualified Data.Vector.Unboxed             as U
import qualified Streaming.Prelude               as Q

import Bed ( mkBedFilter)
import Genome
import Lump
import Util ( parseFile1Opts )

data ConfShameMS = ConfShameMS {
    conf_noutgroups     :: Maybe Int,
    conf_blocklength    :: Maybe Int,
    conf_minInformative :: Double,
    conf_all            :: Bool,
    conf_split          :: Bool,
    conf_reference      :: Maybe FilePath,
    conf_regions        :: Maybe FilePath }
  deriving Show

defaultMsConf :: ConfShameMS
defaultMsConf = ConfShameMS (Just 0) Nothing 90.0 True True Nothing Nothing

optsVcfout :: [ OptDescr (ConfShameMS -> IO ConfShameMS) ]
optsVcfout =
    [ Option "r" ["reference"]     (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "n" ["numoutgroups"]  (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (0)"
    , Option "l" ["blocklength"] (ReqArg set_lblock "NUM") "The length of each block (0)"
    , Option "m" ["minInformative"] (ReqArg set_tblock "NUM") "The minimum informative threshold of each block (0.0-100.0)"
    , Option "D" ["dense"]               (NoArg set_dense) "Output invariant sites, too"
    , Option "t" ["only-transversions"] (NoArg set_no_all) "Output only transversion sites"
    , Option "b" ["only-biallelic"]   (NoArg set_no_split) "Discard, don't split, polyallelic sites"
    , Option "R" ["regions"]      (ReqArg set_rgns "FILE") "Restrict to regions in bed-file FILE" ]
  where
    set_ref    a c = return $ c { conf_reference  =  Just a }
    set_nout   a c = (\n ->   c { conf_noutgroups =  Just n }) <$> readIO a
    set_lblock a c = (\n ->   c { conf_blocklength = Just n }) <$> readIO a
    set_tblock a c = (\n ->   c { conf_minInformative =   n }) <$> readIO a
    set_dense    c = return $ c { conf_noutgroups = Nothing }
    set_no_all   c = return $ c { conf_all        =   False }
    set_no_split c = return $ c { conf_split      =   False }
    set_rgns   a c = return $ c { conf_regions    =  Just a }


mainShameMSout :: [String] -> IO ()
mainShameMSout args = do
    ( hefs, ConfShameMS{..} ) <- parseFile1Opts defaultMsConf "pseudo_ms"
                                                "[hef-file...]" optsVcfout args

    decodeMany conf_reference hefs $ \refs inps -> do
        region_filter <- mkBedFilter conf_regions (either error rss_chroms refs)

        mapsM_ (blocktoShame conf_blocklength conf_minInformative (V.length inps)) $
            maybe yields blocks conf_blocklength $
            region_filter $
            bool singles_only Q.concat conf_split $
            addRef (either error id refs) $
            maybe mergeLumpsDense mergeLumps conf_noutgroups inps
  where
    singles_only = Q.mapMaybe (\case [x] -> Just x ; _ -> Nothing)

    blocktoShame :: Maybe Int -> Double -> Int -> Stream (Of Variant) IO r -> IO r
    blocktoShame ml t m s = do
        ff :> r <- vconcats $ Q.map smashVariants s
        unless (U.null ff) $ do
            
            let tt = [k (ff U.! z) | z <- [ 0 .. U.length ff -1 ] , k <- [ fstW, sndW ]]

            forM_ ml $ \l -> hPutStrLn stdout $ "\n//\nblockSize_" ++ show l
           
            case enoughInfo t m tt of
               Nothing -> return ()
               Just False -> hPutStrLn stdout $ "# Not enough of info in this block"
               Just True  -> hPutStr stdout . unlines $
                          [ [ toRefCode . N2b $ (tt !! i) 
                            | i <- [ j, j+(m*2) .. length tt -1 ] ]
                            | j <- [ 0 .. (m*2)-1 ]]
                            
        return r

    enoughInfo :: Double -> Int -> [Word8] -> Maybe Bool       
    enoughInfo t s list        
        | t < 0 || t> 100   = Nothing
        | otherwise         = Just (containNs * 100 <= (100 - t))
            where containNs = fromIntegral (length . filter (== True) $ map (elem 0xF ) $ LS.chunksOf (s*2) list) / fromIntegral (length list)

    smashVariants :: Variant -> U.Vector Word8
    smashVariants Variant{..} =
        U.map (smashVariant v_ref v_alt) v_calls

    smashVariant :: Nuc2b -> Var2b -> Word8 -> Word8
    smashVariant (N2b r) (V2b a) x
        | x .&. 3 > 0 && x .&. 12 > 0 =       r `pack2` xor r a
        | x .&. 3 > 0                 =       r `pack2` r
        | x .&. 12 > 0                = xor r a `pack2` xor r a
        | otherwise                   =      15 `pack2` 15

    pack2 :: Word8 -> Word8 -> Word8
    pack2 a b = a .|. shiftL b 4

    fstW, sndW :: Word8 -> Word8
    fstW x =        x   .&. 0xF
    sndW x = shiftR x 4 .&. 0xF


blocks :: Monad m => Int -> Stream (Of Variant) m r -> Stream (Stream (Of Variant) m) m r
blocks ln = go (-1) 0
  where
    go c p = lift . Q.next >=> \case
        Left    r                   -> pure r
        Right (v,vs) | v_chr v == c -> yields (Q.span (before     c (p+ln)) (Q.cons v vs)) >>= go     c (p+ln)
                     | otherwise    -> yields (Q.span (before (v_chr v) ln) (Q.cons v vs)) >>= go (v_chr v) ln

    before c p v = v_chr v == c && v_pos v < p


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
