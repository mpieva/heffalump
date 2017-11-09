module ShameMS ( main_shameMSout ) where

import Bio.Prelude
import System.FilePath                  ( takeBaseName )
import System.Console.GetOpt
import Streaming

import qualified Data.ByteString.Builder         as B
import qualified Data.ByteString.Char8           as B
import qualified Data.Vector                     as V
import qualified Data.Vector.Split               as S
import qualified Data.Maybe                      as M
import qualified Data.Vector.Unboxed             as U
import qualified Streaming.Prelude               as Q

import Bed ( mkBedFilter)
import Genome
import Lump
import Util ( parseFileOpts, mk_opts )

data ConfShameMS = ConfShameMS {
    conf_noutgroups :: Maybe Int,
    conf_blocklength:: Maybe Int,
    conf_all        :: Bool,
    conf_split      :: Bool,
    conf_reference  :: Maybe FilePath,
    conf_regions    :: Maybe FilePath }
  deriving Show

defaultMsConf :: ConfShameMS
defaultMsConf = ConfShameMS (Just 0) (Just 0) True True Nothing Nothing

opts_vcfout :: [ OptDescr (ConfShameMS -> IO ConfShameMS) ]
opts_vcfout =
    [ Option "r" ["reference"]     (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "n" ["numoutgroups"]  (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (0)"
    , Option "l" ["blocklength"]   (ReqArg set_lblock "NUM") "The length of each block (0)"
    , Option "D" ["dense"]               (NoArg set_dense) "Output invariant sites, too"
    , Option "t" ["only-transversions"] (NoArg set_no_all) "Output only transversion sites"
    , Option "b" ["only-biallelic"]   (NoArg set_no_split) "Discard, don't split, polyallelic sites"
    , Option "R" ["regions"]      (ReqArg set_rgns "FILE") "Restrict to regions in bed-file FILE" ]
  where
    set_ref    a c = return $ c { conf_reference  =  Just a }
    set_nout   a c = (\n ->   c { conf_noutgroups =  Just n }) <$> readIO a
    set_lblock a c = (\n ->   c { conf_blocklength = Just n }) <$> readIO a
    set_dense    c = return $ c { conf_noutgroups = Nothing }
    set_no_all   c = return $ c { conf_all        =   False }
    set_no_split c = return $ c { conf_split      =   False }
    set_rgns   a c = return $ c { conf_regions    =  Just a }


main_shameMSout :: [String] -> IO ()
main_shameMSout args = do
    ( hefs, ConfShameMS{..} ) <- parseFileOpts defaultMsConf
                                              (mk_opts "pseudo_ms" "[hef-file...]" opts_vcfout) args

    decodeMany conf_reference hefs $ \refs inps -> do
        region_filter <- mkBedFilter conf_regions (either error rss_chroms refs)
        --window_filter <- mkwindowFilter (Just 1) conf_blocklength 
        let chrs = either error (V.fromList . rss_chroms) refs

        let the_vars = addRef (either error id refs) $
                       blocks (Just 1) conf_blocklength $ --window_filter $
                       region_filter $
                       bool singles_only Q.concat conf_split $
                       maybe mergeLumpsDense mergeLumps conf_noutgroups inps
        
        my_vec <- U.fromList . concat <$> Q.toList_ (Q.map smashVariants the_vars)


        let m = length inps
        let res = [[ my_vec U.! i | i <- [ 2 * j + o, 2 * (j+m) + o .. U.length my_vec-1] ] | j <- [ 0 .. m-1 ] , o <- [0,1] ]
        

        hPutStr stdout . unlines $ map (map (toRefCode . N2b)) res

  where
    singles_only = Q.concat . Q.map (\case [x] -> Just x ; _ -> Nothing)

    smashVariants :: Variant -> [Word8]
    smashVariants Variant{..} = 
        concatMap (smashVariant v_ref v_alt) (U.toList v_calls)

    smashVariant :: Nuc2b -> Var2b -> Word8 -> [Word8]
    smashVariant (N2b r) (V2b a) x
  		| x .&. 3 > 0 && x .&. 12 > 0 = [r,xor r a]
  		| x .&. 3 > 0 = [r,r]
  		| x .&. 12 > 0 = [xor r a,xor r a]
  		| x == 0 = [255,255]

myGroupBy :: Monad m => Maybe Int -> Maybe Int -> Stream (Of Variant) m r -> Stream (Stream (Of Variant) m ) m r
myGroupBy s w = groupBy (\l g -> v_pos l >=  fromJust s ) -- $ S.each $ V.fromList [(1,'a'),(2,'b'),(3,'c'),(2,'b'),(3,'e'),(5,'f'),(4,'g'),(6,'h')]

--mblocks :: Monad m => Maybe Int -> Int -> Bool
--mblocks =
--  where
--    go _      Nothing = lift . Q.effects >=> pure
--    go Nothing _      = lift . Q.effects >=> pure
--    go s       w      = lift . Q.next >=> \case

blocks :: Monad m => Maybe Int -> Maybe Int -> Stream (Of Variant) m r -> Stream (Of Variant) m r --Stream (Stream (Of Variant) m ) m r
blocks = go 
  where
    go _      Nothing = lift . Q.effects >=> pure
    go Nothing _      = lift . Q.effects >=> pure
    go s       w      = lift . Q.next >=> \case
        Left   r    -> pure r
        Right (v,vs)
            -- variant belong to this window keep it
            | v_pos v <  fromJust s -> v `Q.cons` go s w vs

            -- variant passes the window size, make next boundry
            | v_pos v >= fromJust((+) <$> s <*> w) -> v `Q.cons` go ((+) <$> s <*> w) w vs

            -- when logic fails!
            | otherwise                     -> v `Q.cons` go Nothing Nothing vs


