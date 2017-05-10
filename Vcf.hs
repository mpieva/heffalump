module Vcf ( main_xcf, conf_vcf, conf_bcf ) where

-- ^ Stuff related to Vcf and Bcf

import BasePrelude
import Streaming
import System.Console.GetOpt
import System.IO

import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Streaming       as S

import BcfScan
import Lump
import NewRef
import Util
import VcfScan

type LumpXform = (Stream Lump IO () -> Stream Lump IO ())
type GapCons = Int -> () -> Lump ()

data ConfXcf = ConfXcf
    { conf_output  :: FilePath
    , conf_ref     :: FilePath
    , conf_density :: GapCons
    , conf_ploidy  :: LumpXform
    , conf_key     :: String
    , conf_ext     :: String
    , conf_reader  :: FilePath -> IO [RawVariant] }

conf_vcf :: ConfXcf
conf_vcf = ConfXcf (error "no output file specified")
                   (error "no reference specified")
                   Ns id "vcfin" "vcf" readVcf

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

    set_dense    c = return $ c { conf_density = Ns       }
    set_sparse   c = return $ c { conf_density = Eqs2     }
    set_low      c = return $ c { conf_ploidy  = make_hap }
    set_high     c = return $ c { conf_ploidy  = id       }


main_xcf :: ConfXcf -> [String] -> IO ()
main_xcf conf0 args = do
    ( vcfs, ConfXcf{..} ) <- let opts = mk_opts (conf_key conf0) fspec opts_xcf
                                 fspec = "["++conf_ext conf0++"-file...]"
                             in parseOpts True conf0 opts args
    ref <- readTwoBit conf_ref
    withFile conf_output WriteMode $ \hdl ->
        S.hPut hdl . encodeLump' . conf_ploidy
        . importVcf conf_density (nrss_chroms ref)
        . progress conf_output . dedupVcf . cleanVcf
        . concat =<< mapM conf_reader vcfs
  where
    progress :: String -> [RawVariant] -> [RawVariant]
    progress fp = go 0 0
      where
        go  _  _ [    ] = []
        go rs po (v:vs)
            | rs /= rv_chrom v || po + 10000000 <= rv_pos v
                = unsafePerformIO $ do
                    hPutStrLn stderr $ fp ++ "@" ++ show (rv_chrom v) ++ ":" ++ show (rv_pos v)
                    return $ v : go (rv_chrom v) (rv_pos v) vs
            | otherwise =  v : go rs po vs

-- | Some idiot decided to output multiple records for the same position
-- into some VCF files.  If we hit that, we take the first.  (Einmal mit
-- Profis arbeiten!)
dedupVcf :: [ RawVariant ] -> [ RawVariant ]
dedupVcf (v1:v2:vs)
    | rv_chrom v1 == rv_chrom v2 && rv_pos v1 == rv_pos v2  =      dedupVcf (v1:vs)
    | otherwise                                             = v1 : dedupVcf (v2:vs)
dedupVcf [v1] = [v1]
dedupVcf [  ] = [  ]

-- | Remove indel variants, since we can't very well use them.  Also
-- removes "no call" variants, which are equivalent to missing entries.
cleanVcf :: [ RawVariant ] -> [ RawVariant ]
cleanVcf = filter $ \RawVariant{..} ->
    B.length rv_vars  == B.length (B.filter (== ',') rv_vars) * 2 + 1
    && rv_gt /= 0x0000 && rv_gt /= 0xff00

-- Reads a VCF file and returns 'RawVariant's, which is not exactly
-- practical.
readVcf :: FilePath -> IO [ RawVariant ]
readVcf fp = do
    sc <- initVcf fp
    fix $ \self -> unsafeInterleaveIO $ do
        mv <- getVariant sc
        case mv of
            Nothing -> return []
            Just  v -> (:) v <$> self


importVcf :: Monad m => GapCons -> [ B.ByteString ] -> [ RawVariant ] -> Stream Lump m ()
importVcf ns = (.) normalizeLump' . go
  where
    -- We start the coordinate at one(!), for VCF is one-based.
    -- Pseudo-variants of the "no call" type must be filtered
    -- beforehand, see 'cleanVcf'.
    go [    ] = const $ yields (Break ())
    go (c:cs) = generic (hashChrom c) 1
      where
        generic !_    _ [         ] = yields (Break ())
        generic !hs pos (var1:vars)
            | hs /= rv_chrom var1   = yields (Break ()) >> go cs (var1:vars)

            -- long gap, creates Ns or Eqs2
            | rv_pos var1 > pos     = yields (ns (rv_pos var1 - pos) ()) >> generic hs (rv_pos var1) (var1:vars)

            -- positions must match now
            | rv_pos var1 == pos =
                        if isVar var1
                          then yields (get_var_code var1 ()) >> generic hs (succ pos) vars
                          else yields (Eqs2     1        ()) >> generic hs (succ pos) vars

            | otherwise = error $ "Got variant position " ++ show (rv_pos var1)
                               ++ " when expecting " ++ show pos ++ " or higher."

    -- *sigh*  Have to turn a numeric genotype into a 'NucCode'.  We
    -- have characters for the variants, and we need to map a pair of
    -- them to a code.
    get_var_code RawVariant{..}
        | B.any (== 'N') rv_vars || B.null rv_vars = Ns 1

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
    char_to_2b 'A' =  1
    char_to_2b 'C' =  2
    char_to_2b 'G' =  3
    char_to_2b  c  = error $ "What's a " ++ shows c "?"

    isVar RawVariant{..} | rv_gt == 0xFF02            = False     -- "0"
                         | rv_gt .&. 0xFCFE == 0x0002 = False     -- "0|.", "0/.", "0|0", "0/0"
                         | otherwise                  = True

