module ShameMS ( mainShameMSout ) where

import Bio.Prelude
import System.Console.GetOpt
import Streaming

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
defaultMsConf = ConfShameMS (Just 0) Nothing True True Nothing Nothing

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

        let the_vars = maybe yields blocks conf_blocklength $
                       addRef (either error id refs) $
                       region_filter $
                       bool singles_only Q.concat conf_split $
                       maybe mergeLumpsDense mergeLumps conf_noutgroups inps

        Q.mapM_ (blocktoShame (length inps) . U.fromList . concatMap smashVariants) $
                Q.mapped Q.toList the_vars
  where
    singles_only = Q.concat . Q.map (\case [x] -> Just x ; _ -> Nothing)

    blocktoShame :: Int -> U.Vector Word8 -> IO ()
    blocktoShame m ff = unless (U.null ff) $
                            hPutStr stdout $ unlines
                                [ [ toRefCode . N2b $ k (ff U.! i)
                                  | i <- [ j, j+m .. U.length ff -1 ] ]
                                | j <- [ 0 .. m-1 ], k <- [ fstW, sndW ] ]

    smashVariants :: Variant -> [Word8]
    smashVariants Variant{..} =
        map (smashVariant v_ref v_alt) (U.toList v_calls)

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
