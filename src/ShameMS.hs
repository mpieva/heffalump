module ShameMS ( mainShameMSout ) where

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
import qualified Data.Typeable                    as D

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

optsVcfout :: [ OptDescr (ConfShameMS -> IO ConfShameMS) ]
optsVcfout =
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


mainShameMSout :: [String] -> IO ()
mainShameMSout args = do
    ( hefs, ConfShameMS{..} ) <- parseFileOpts defaultMsConf
                                              (mk_opts "pseudo_ms" "[hef-file...]" optsVcfout) args

    decodeMany conf_reference hefs $ \refs inps -> do
        region_filter <- mkBedFilter conf_regions (either error rss_chroms refs)
        let chrs = either error (V.fromList . rss_chroms) refs

        let the_vars = blocks conf_blocklength $ 
                       addRef (either error id refs) $
                       region_filter $
                       bool singles_only Q.concat conf_split $
                       maybe mergeLumpsDense mergeLumps conf_noutgroups inps
      
        let m = length inps
        --my_vec <-  U.fromList . concat <$> Q.toList_ (Q.map smashVariantss $ Q.mapped Q.toList the_vars)
        my_vec <-  Q.toList_ $ Q.map smashVariantss $ Q.mapped Q.toList  the_vars
        
        mapM_ (blocktoShame m . U.fromList) my_vec

  where
    singles_only = Q.concat . Q.map (\case [x] -> Just x ; _ -> Nothing)

    blocktoShame :: Int -> U.Vector Word8 -> IO ()
    blocktoShame m ff = hPutStr stdout . unlines $ map (map (toRefCode . N2b)) [[ ff U.! i | i <- [ 2 * j + o, 2 * (j+m) + o .. U.length ff -1] ] | j <- [ 0 .. m-1 ] , o <-[0,1] ]  
   
    smashVariantss :: [Variant] -> [Word8] 
    smashVariantss = concatMap smashVariants

    smashVariants :: Variant -> [Word8]
    smashVariants Variant{..} = 
        concatMap (smashVariant v_ref v_alt) (U.toList v_calls)

    smashVariant :: Nuc2b -> Var2b -> Word8 -> [Word8]
    smashVariant (N2b r) (V2b a) x
        | x .&. 3 > 0 && x .&. 12 > 0 = [r,xor r a]
        | x .&. 3 > 0                 = [r,r]
        | x .&. 12 > 0                = [xor r a,xor r a]
        | x == 0                      = [255,255]
        | otherwise                   = [255,255]

blocks :: Monad m => Maybe Int -> Stream (Of Variant) m r -> Stream (Stream (Of Variant) m) m r
blocks ln = go (-1) (Just 0)
  where
    go c p = lift . Q.next >=> \case
        Left    r                   -> pure r
        Right (v,vs) | v_chr v == c -> yields (Q.span (before     c ((+) <$> ln <*> p)) (Q.cons v vs)) >>= go     c  ((+) <$> ln <*> p)
                     | otherwise    -> yields (Q.span (before (v_chr v) ln) (Q.cons v vs)) >>= go (v_chr v) ln

    before c p v = v_chr v == c && v_pos v < fromJust p
