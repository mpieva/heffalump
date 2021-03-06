module VcfOut ( main_vcfout ) where

import Bio.Prelude
import System.FilePath                  ( takeBaseName )
import System.Console.GetOpt

import qualified Data.ByteString.Builder         as B
import qualified Data.ByteString.Char8           as B
import qualified Data.Vector                     as V
import qualified Data.Vector.Storable            as U
import qualified Streaming.Prelude               as Q

import Bed ( mkBedFilter )
import Genome
import Lump
import Util ( parseFile1Opts )

data ConfVcfOut = ConfVcfOut {
    conf_noutgroups :: Maybe Int,
    conf_all        :: Bool,
    conf_split      :: Bool,
    conf_reference  :: Maybe FilePath,
    conf_regions    :: Maybe FilePath }
  deriving Show

defaultVcfConf :: ConfVcfOut
defaultVcfConf = ConfVcfOut (Just 0) True True Nothing Nothing

opts_vcfout :: [ OptDescr (ConfVcfOut -> IO ConfVcfOut) ]
opts_vcfout =
    [ Option "r" ["reference"]     (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "n" ["numoutgroups"]  (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (0)"
    , Option "D" ["dense"]               (NoArg set_dense) "Output invariant sites, too"
    , Option "t" ["only-transversions"] (NoArg set_no_all) "Output only transversion sites"
    , Option "b" ["only-biallelic"]   (NoArg set_no_split) "Discard, don't split, polyallelic sites"
    , Option "R" ["regions"]      (ReqArg set_rgns "FILE") "Restrict to regions in bed-file FILE" ]
  where
    set_ref    a c = return $ c { conf_reference  =  Just a }
    set_nout   a c = (\n ->   c { conf_noutgroups =  Just n }) <$> readIO a
    set_dense    c = return $ c { conf_noutgroups = Nothing }
    set_no_all   c = return $ c { conf_all        =   False }
    set_no_split c = return $ c { conf_split      =   False }
    set_rgns   a c = return $ c { conf_regions    =  Just a }


main_vcfout :: [String] -> IO ()
main_vcfout args = do
    ( hefs, ConfVcfOut{..} ) <- parseFile1Opts defaultVcfConf "vcfexport"
                                               "[hef-file...]" opts_vcfout args

    decodeMany conf_reference hefs $ \refs inps -> do
        region_filter <- mkBedFilter conf_regions (either error rss_chroms refs)
        let chrs = either error (V.fromList . rss_chroms) refs
        B.putStr $ "##fileformat=VCFv4.1\n" <>
                   "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" <>
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" <>
                   mconcat [ '\t' `B.cons` B.pack (takeBaseName f) | f <- hefs ] <>
                   B.singleton '\n'

        let the_vars = region_filter $
                       bool singles_only Q.concat conf_split $
                       maybe (mergeLumpsRef (either error id refs))
                             (\n -> addRef (either error id refs) . mergeLumps n)
                             conf_noutgroups inps

        flip Q.mapM_ the_vars $ \Variant{..} ->
            -- samples (not outgroups) must show alt allele at least once
            when (conf_all || isTransversion v_alt) $
                B.hPutBuilder stdout $
                    B.byteString (chrs V.! v_chr) <> B.char8 '\t' <>
                    B.intDec (v_pos+1) <> B.string8 "\t.\t" <>
                    B.char8 (toRefCode v_ref) <> B.char8 '\t' <>
                    B.char8 (toAltCode v_alt v_ref) <> B.string8 "\t.\t.\t.\tGT" <>
                    U.foldr ((<>) . gts) mempty v_calls <>
                    B.char8 '\n'
  where
    singles_only = Q.concat . Q.map (\case [x] -> Just x ; _ -> Nothing)

    gts :: AlleleCounts -> B.Builder
    gts (AC 0 0) = "\t./."      -- N
    gts (AC 1 0) = "\t0"        -- 1x ref
    gts (AC 0 1) = "\t1"        -- 1x alt
    gts (AC _ 0) = "\t0/0"      -- >1x ref
    gts (AC 0 _) = "\t1/1"      -- >1x alt
    gts (AC _ _) = "\t0/1"      -- ref+alt

