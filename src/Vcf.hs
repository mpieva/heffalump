module Vcf ( main_xcf, conf_vcf, conf_bcf ) where

-- ^ Stuff related to Vcf and Bcf

import Bio.Prelude hiding ( Ns )
import Streaming
import System.Console.GetOpt
import System.IO

import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Streaming       as S
import qualified Streaming.Prelude               as Q

import BcfScan
import Lump
import NewRef
import Util
import VcfScan

type LumpXform = (Stream (Of Lump) IO () -> Stream (Of Lump) IO ())
type GapCons = Int -> Lump

data ConfXcf = ConfXcf
    { conf_output  :: FilePath
    , conf_ref     :: FilePath
    , conf_density :: GapCons
    , conf_clean   :: Stream (Of RawVariant) IO () -> Stream (Of RawVariant) IO ()
    , conf_ploidy  :: LumpXform
    , conf_key     :: String
    , conf_ext     :: String
    , conf_reader  :: FilePath -> (Stream (Of RawVariant) IO () -> IO ()) -> IO () }

conf_vcf :: ConfXcf
conf_vcf = ConfXcf (error "no output file specified")
                   (error "no reference specified")
                   (error "specify either --dense or --sparse")
                   id id "vcfin" "vcf" readVcf

conf_bcf :: ConfXcf
conf_bcf = conf_vcf { conf_key = "bcfin", conf_ext = "bcf", conf_reader = readBcf }


opts_xcf :: [ OptDescr ( ConfXcf -> IO ConfXcf ) ]
opts_xcf =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "D" ["dense"]          (NoArg  set_dense) "Input stores all sites"
    , Option "S" ["sparse"]         (NoArg set_sparse) "Input stores only variants"
    , Option "L" ["low"]            (NoArg    set_low) "Low coverage, stores haploid calls"
    , Option "H" ["high"]           (NoArg   set_high) "High coverage, stores diploid calls" ]
  where
    set_output a c = return $ c { conf_output  = a }
    set_ref    a c = return $ c { conf_ref     = a }

    set_dense    c = return $ c { conf_density = Ns,   conf_clean = cleanMissing }
    set_sparse   c = return $ c { conf_density = Eqs2, conf_clean = id           }
    set_low      c = return $ c { conf_ploidy  = make_hap }
    set_high     c = return $ c { conf_ploidy  = id       }


main_xcf :: ConfXcf -> [String] -> IO ()
main_xcf conf0 args = do
    ( vcfs, ConfXcf{..} ) <- let opts = mk_opts (conf_key conf0) fspec opts_xcf
                                 fspec = "["++conf_ext conf0++"-file...]"
                             in parseFileOpts conf0 opts args
    ref <- readTwoBit conf_ref
    withMany conf_reader vcfs $ \inp ->
        withFile conf_output WriteMode $ \hdl ->
            S.hPut hdl . encode ref . conf_ploidy
            . importVcf conf_density (nrss_chroms ref)
            . progress conf_output . dedupVcf
            . conf_clean . cleanVcf $ inp
  where
    progress :: MonadIO m => String -> Stream (Of RawVariant) m r -> Stream (Of RawVariant) m r
    progress fp = go 0 0
      where
        go rs po = lift . Q.next >=> \case
          Left r -> pure r
          Right (v,vs)
            | rs /= rv_chrom v || po + 10000000 <= rv_pos v -> do
                    liftIO $ hPutStrLn stderr $ fp ++ "@" ++ show (rv_chrom v) ++ ":" ++ show (rv_pos v)
                    v `Q.cons` go (rv_chrom v) (rv_pos v) vs
            | otherwise ->  v `Q.cons` go rs po vs

    withMany _ [ ] k = k $ pure ()
    withMany r (fp:fps) k = r fp $ \s -> withMany r fps $ k . (>>) s

-- | Some idiot decided to output multiple records for the same position
-- into some VCF files.  If we hit that, we take the first.  (Einmal mit
-- Profis arbeiten!)
dedupVcf :: MonadIO m => Stream (Of RawVariant) m r -> Stream (Of RawVariant) m r
dedupVcf = lift . Q.next >=> \case Left       r  -> pure r
                                   Right (v1,vs) -> go v1 vs
  where
    go v1 = lift . Q.next >=> \case Left r                      -> v1 `Q.cons` pure r
                                    Right (v2,vs) | match v1 v2 ->             go v1 vs
                                                  | otherwise   -> v1 `Q.cons` go v2 vs

    match v1 v2 = rv_chrom v1 == rv_chrom v2 && rv_pos v1 == rv_pos v2


-- | Remove indel variants, since we can't very well use them.
cleanVcf :: Monad m => Stream (Of RawVariant) m r -> Stream (Of RawVariant) m r
cleanVcf = Q.filter $ \RawVariant{..} ->
    B.length rv_vars  == B.length (B.filter (== ',') rv_vars) * 2 + 1

-- | Removes "no call" variants.  When reading dense files, they are
-- equivalent to missing entries.
cleanMissing :: Monad m => Stream (Of RawVariant) m r -> Stream (Of RawVariant) m r
cleanMissing = Q.filter $ \RawVariant{..} ->
    rv_gt /= 0x0000 && rv_gt /= 0xff00

-- Reads a VCF file and returns 'RawVariant's, which is not exactly
-- practical.
readVcf :: FilePath -> (Stream (Of RawVariant) IO () -> IO r) -> IO r
readVcf fp k = do
    sc <- initVcf fp
    k $ Q.untilRight (getVariant sc)

importVcf :: Monad m => GapCons -> [ B.ByteString ] -> Stream (Of RawVariant) m r -> Stream (Of Lump) m r
importVcf ns = (.) normalizeLump . go
  where
    -- We start the coordinate at one(!), for VCF is one-based.
    -- Pseudo-variants of the "no call" type must be filtered
    -- beforehand, see 'cleanVcf'.
    go [    ] = Q.cons Break . lift . Q.effects
    go (c:cs) = generic (hashChrom c) 1
      where
        generic !hs pos = lift . Q.next >=> \case
            Left r                    -> Q.cons Break $ pure r
            Right (var1,vars)
                | hs /= rv_chrom var1 -> Q.cons Break $ go cs (var1 `Q.cons` vars)

                -- long gap, creates Ns or Eqs2
                | rv_pos var1 > pos   -> Q.cons (ns (rv_pos var1 - pos)) $
                                         generic hs (rv_pos var1) (var1 `Q.cons` vars)

                -- positions must match now
                | rv_pos var1 < pos   -> error $ "Got variant position " ++ show (rv_pos var1)
                                              ++ " when expecting " ++ show pos ++ " or higher."

                | isVar var1          -> Q.cons (get_var_code var1) $ generic hs (succ pos) vars
                | otherwise           -> Q.cons (Eqs2     1       ) $ generic hs (succ pos) vars


    -- *sigh*  Have to turn a numeric genotype into a 'NucCode'.  We
    -- have characters for the variants, and we need to map a pair of
    -- them to a code.
    get_var_code RawVariant{..}
        | B.any (== 'N') rv_vars || B.null rv_vars = Ns 1

        -- missing call
        | rv_gt == 0xFF00 || rv_gt == 0x0000 = Ns 1

        -- haploid call or one allele missing
        | rv_gt .&. 0xFF00 == 0xFF00 || rv_gt .&. 0xFF00 == 0x0000 =
                encTwoVars $ 16 + (char_to_2b c1 `xor` char_to_2b n0)

        -- diploid call; one called allele must be the reference
        | otherwise = encTwoVars $ (char_to_2b c1 `xor` char_to_2b n0)
                                 + (char_to_2b c2 `xor` char_to_2b n0) * 4

        -- diploid call, but unrepresentable
        | otherwise = Ns 1
      where
        v1 = fromIntegral $ rv_gt            .&. 0x00FE - 2
        v2 = fromIntegral $ rv_gt `shiftR` 8 .&. 0x00FE - 2

        n0 = toUpper $ B.head rv_vars
        c1 = toUpper $ safeIndex "c1" rv_vars v1
        c2 = toUpper $ safeIndex ("c2 "++shows rv_pos " " ++ showHex rv_gt " ") rv_vars v2

    safeIndex m s i | B.length s > i = B.index s i
                    | otherwise = error $ "Attempted to index " ++ shows i " in " ++ shows s " (" ++ m ++ ")."

    char_to_2b 'T' =  0
    char_to_2b 'C' =  1
    char_to_2b 'A' =  2
    char_to_2b 'G' =  3
    char_to_2b  c  = error $ "What's a " ++ shows c "?"

    isVar RawVariant{..} | rv_gt == 0xFF02            = False     -- "0"
                         | rv_gt .&. 0xFCFE == 0x0002 = False     -- "0|.", "0/.", "0|0", "0/0"
                         | otherwise                  = True

