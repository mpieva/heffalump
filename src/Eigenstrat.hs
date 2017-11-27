module Eigenstrat ( main_eigenstrat ) where

import Bio.Prelude
import System.Console.GetOpt
import System.IO

import qualified Data.ByteString.Char8          as B
import qualified Data.Vector.Unboxed            as U
import qualified Streaming.Prelude              as Q

import Bed
import Genome
import Lump
import Util

-- ^ The Eigenstrat and Ancestrymap exporters.  Includes generation of
-- meaningless names and hash functions.


-- We wrote Eigenstrat in the first iteration.  Eigenstrat is one
-- character per genotype.  Alternatively, we could write packed
-- Ancestrymap, then it's one char per four genotypes.  It starts with a
-- header record, which contains text and is then padded to the record
-- size as defined below:
-- printf("GENO %d %d %x %x", n_individuals, n_snps, individual_hash, snp_hash);
--
-- That's of course retarded, we can't know n_snps or snp_hash before we
-- generate the snp file.  So all we can do is reserve space and
-- overwrite the header later.  m(
--
-- To generate individual_hash, we apply nick_hasharr to the vector of
-- individual names; to generate snp_hash we apply nick_hasharr to the
-- vector of SNP identifiers.
--
-- The remainder is binary:  one record per snp, each record uses 2 bits
-- per individual, padded to a full byte, padded to at least 48 bytes.


-- Not used anywhere
-- Nickhash for strings.  (It's not clear if we got the signs right?)
-- nick_hashit :: L.ByteString -> Int32
-- nick_hashit = L.foldl (\h c -> 23 * h + fromIntegral c) 0

-- Not used anywhere
-- Nickhash for vectors of strings.  (This will need incremental
-- updates, but that's easy.)
-- nick_hasharr :: F.Foldable v => v L.ByteString -> Int32
-- nick_hasharr = F.foldl (\h s -> 17 * h `xor` nick_hashit s) 0

mkname :: Int -> Int -> Char -> String
mkname x y z = enc (4*y + numOf z) [ B.index chars x ]
  where
    numOf 'A' = 0 ; numOf 'C' = 1 ; numOf 'G' = 2 ; numOf 'T' = 3
    numOf 'a' = 0 ; numOf 'c' = 1 ; numOf 'g' = 2 ; numOf  _  = 3

    chars = "0123456789ABCDEFGHKLMNPQRSTUWXYZ"

    enc 0 = id
    enc u = (:) (B.index chars (u .&. 31)) . enc (u `shiftR` 5)

data ConfEigen = ConfEigen {
    conf_noutgroups :: Int,
    conf_blocksize  :: Int,
    conf_all        :: Bool,
    conf_split      :: Bool,
    conf_reference  :: Maybe FilePath,
    conf_regions    :: Maybe FilePath,
    conf_nrefpanel  :: Int,
    conf_output     :: FilePath,
    conf_sample     :: FilePath }
  deriving Show

defaultConfEigen :: ConfEigen
defaultConfEigen = ConfEigen 0 5000000 True True Nothing Nothing
                                (error "size of reference panel not known")
                                (error "no output file specified")
                                (error "no sample file specified")

opts_eigen :: [ OptDescr (ConfEigen -> IO ConfEigen) ]
opts_eigen =
    [ Option "o" ["output"]     (ReqArg set_output "FILE") "Write output to FILE.geno and FILE.snp"
    , Option "r" ["reference"]     (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "n" ["numoutgroups"]  (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (0)"
    , Option "t" ["only-transversions"] (NoArg set_no_all) "Output only transversion sites"
    , Option "b" ["only-biallelic"]   (NoArg set_no_split) "Discard, don't split, polyallelic sites"
    , Option "R" ["regions"]      (ReqArg set_rgns "FILE") "Restrict to regions in bed-file FILE" ]
  where
    set_output a c = return $ c { conf_output     =      a }
    set_ref    a c = return $ c { conf_reference  = Just a }
    set_nout   a c = (\n ->   c { conf_noutgroups =      n }) <$> readIO a
    set_no_all   c = return $ c { conf_all        =  False }
    set_no_split c = return $ c { conf_split      =  False }
    set_rgns   a c = return $ c { conf_regions    = Just a }

-- merge multiple files with the reference, write Eigenstrat format (geno & snp files)
main_eigenstrat :: [String] -> IO ()
main_eigenstrat args = do
    ( hefs, ConfEigen{..} ) <- parseFile1Opts defaultConfEigen "eigenstrat"
                                              "[hef-file...]" opts_eigen args

    decodeMany conf_reference hefs $ \refs inps -> do
      region_filter <- mkBedFilter conf_regions (either error rss_chroms refs)
      withFile (conf_output ++ ".snp") WriteMode $ \hsnp ->
        withFile (conf_output ++ ".geno") WriteMode $ \hgeno -> do
            let vars = either (const id) addRef refs $
                       region_filter $
                       bool singles_only Q.concat conf_split $
                       mergeLumps conf_noutgroups inps
            flip Q.mapM_ vars $ \Variant{..} ->
                -- samples (not outgroups) must show ref and alt allele at least once
                let ve = U.foldl' (.|.) 0 $ U.drop conf_noutgroups v_calls
                    is_ti = conf_all || isTransversion v_alt in
                when (ve .&. 3 /= 0 && ve .&. 12 /= 0 && is_ti) $ do
                    hPutStrLn hgeno $ map (B.index "9222011101110111" . fromIntegral) $ U.toList v_calls
                    hPutStrLn hsnp $ intercalate "\t"
                        -- 1st column is SNP name
                        [ mkname v_chr v_pos (toAltCode v_alt v_ref)
                        -- "2nd column is chromosome.  X chromosome is encoded as 23.
                        -- Also, Y is encoded as 24, mtDNA is encoded as 90, ..."
                        , show $ if v_chr == 24 then 90 else v_chr + 1
                        -- "3rd column is genetic position (in Morgans).
                        -- If unknown, ok to set to 0.0"
                        , "0.0"
                        -- "4th column is physical position (in bases)"
                        , show (v_pos+1)
                        -- "Optional 5th and 6th columns are reference and variant alleles"
                        , [toRefCode v_ref], [toAltCode v_alt v_ref] ]
  where
    singles_only = Q.concat . Q.map (\case [x] -> Just x ; _ -> Nothing)


