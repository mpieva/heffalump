import Bio.Prelude
import Control.Monad.Trans.Resource
import Paths_heffalump                  ( version )
import Streaming
import System.Console.GetOpt
import System.IO

import qualified Data.ByteString.Builder         as B
import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.ByteString.Streaming       as S
import qualified Data.IntMap                     as I
import qualified Streaming.Prelude               as Q

import Bamin
import Eigenstrat
import Emf
import Genome
import Lump
import SillyStats
import Stretch ( main_dumppatch )
import Treemix
import Util
import Vcf
import VcfOut
import ShameMS

main :: IO ()
main = do
    as0 <- getArgs
    case as0 of
        [] -> usage
        a:as -> case [ m | (k,(m,_)) <- mains
                     , map toLower a `isPrefixOf` map toLower k ] of
            [m] -> m as
            _   -> usage
  where
    usage = do
        pn <- getProgName
        let kl = maximum $ map (length . fst) mains
        hPutStrLn stderr . unlines $ "Usage:" :
              [ "  "++pn++' ':k++replicate (kl-length k) ' '++"  -- "++h
              | (k,(_,h)) <- mains ]

    mains = let z a b c = (a,(b,c)) in
        [ z "hetfa"        main_hetfa             "Import hetfa file"
        , z "maf"          main_maf               "Import one genome from maf files"
        , z "emf"          main_emf               "Import one genome from emf (Compara)"
        , z "bcfin"       (main_xcf conf_bcf)     "Import one genome from bcf files"
        , z "vcfin"       (main_xcf conf_vcf)     "Import one genome from vcf files"
        , z "bamin"        main_bam               "Import crudely from bam files"
        , z "patch"        main_patch             "Make a hetfa file by patching the reference"
        , z "treemix"      main_treemix           "Merge heffalumps into Treemix format"
        , z "eigenstrat"   main_eigenstrat        "Merge heffalumps into Eigenstrat format"
        , z "vcfexport"    main_vcfout            "Merge heffalumps into vcf format"
        , z "pseudo_ms"    mainShameMSout         "Merge heffalumps into pseudo_ms format"
        , z "twobitinfo"   main_2bitinfo          "List reference sequences"
        , z "twobittofa"   main_2bittofa          "Extract Fasta from 2bit"
        , z "fatotwobit"   main_fato2bit          "Convert Fasta to 2bit"
        , z "kayvergence"  main_kayvergence       "Compute Kayvergence ratios"
        , z "f4statistics"(main_patterson WantF4) "Compute Patterson's F4 and D"
        , z "dstatistics" (main_patterson NoF4)   "Compute Patterson's D"
        , z "yaddayadda"   main_yaddayadda        "Compute Yadda-Yadda-counts"
        , z "dumppatch"    main_dumppatch         "(debugging aid)"
        , z "dumplump"     main_dumplump          "(debugging aid)"
        , z "--help"      (const usage)           "list commands and exit"
        , z "--version"    main_version           "print version and exit" ]

main_version :: [String] -> IO ()
main_version _ = hPutStrLn stderr $ showVersion version

data ConfImportGen = ConfImportGen {
    conf_imp_reference  :: FilePath,
    conf_imp_output     :: FilePath,
    conf_imp_sample     :: FilePath }
  deriving Show

defaultImportConf :: ConfImportGen
defaultImportConf = ConfImportGen (error "no reference specified")
                                  (error "no output file specified")
                                  (error "no sample file specified")


opts_hetfa :: [ OptDescr (ConfImportGen -> IO ConfImportGen) ]
opts_hetfa =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "s" ["sample"] (ReqArg set_sample "FILE") "Read sample from FILE (.hetfa)" ]
  where
    set_output a c = return $ c { conf_imp_output    = a }
    set_ref    a c = return $ c { conf_imp_reference = a }
    set_sample a c = return $ c { conf_imp_sample    = a }


main_hetfa :: [String] -> IO ()
main_hetfa args = do
    ConfImportGen{..} <- parseOpts defaultImportConf (mk_opts "hetfa" [] opts_hetfa) args
    refs <- readTwoBit conf_imp_reference

    runResourceT $ S.writeFile conf_imp_output $
            importHetfa refs $ readSampleFa conf_imp_sample


-- | Import Hetfa by diffing it against a reference.  We must read the
-- Hetfa in the order of the file (Fasta doesn't have an index), but
-- write it in the order of the reference.  So we parse a chromosome at
-- a time, sort the diffs in memory and dump it at the end.
--
-- Chromosomes in input that don't have a match in the reference are
-- silently ignored.  Chromsomes that don't show up in the input are
-- emitted as empty sequences.  Only if no sequence matches do we
-- produce an error message.
importHetfa :: MonadIO m => RefSeqs -> Stream (FastaSeq m) m r -> S.ByteString m r
importHetfa ref smps = S.mwrap $ do
    (map1,r) <- fold_1 I.empty smps

    when (I.null map1) $ error
            "Found only unexpected sequences.  Is this the right reference?"

    return $ do S.fromLazy $ B.toLazyByteString $ encodeHeader ref
                forM_ [0 .. length (rss_chroms ref) - 1] $
                        \i -> S.fromLazy $ unpackLump $ I.findWithDefault noLump i map1
                pure r
  where
    enc2 :: MonadIO m => Int -> S.ByteString m r -> m (Of PackedLump r)
    enc2 i sq = encodeLumpToMem $ diff2 (rss_seqs ref !! i) sq <* Q.yield Break

    fold_1 :: MonadIO m => I.IntMap PackedLump -> Stream (FastaSeq m) m r -> m (I.IntMap PackedLump, r)
    fold_1 !acc s = inspect s >>= \case
        Left r -> return (acc,r)
        Right (FastaSeq nm sq) -> case elemIndex nm (rss_chroms ref) of
            Nothing -> S.effects sq >>= fold_1 acc
            Just  i -> do lump :> s' <- enc2 i sq
                          liftIO $ hPrint stderr (i, L.length $ unpackLump lump)
                          fold_1 (I.insert i lump acc) s'

data ConfMaf = ConfMaf
    { conf_maf_output :: FilePath
    , conf_maf_reference :: FilePath
    , conf_maf_ref_species :: B.ByteString
    , conf_maf_oth_species :: B.ByteString }

opts_maf :: [ OptDescr ( ConfMaf -> IO ConfMaf ) ]
opts_maf =
    [ Option "o" ["output"]     (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"]     (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "R" ["ref-species"]  (ReqArg set_from "NAME") "Set reference species to NAME"
    , Option "S" ["sample-species"] (ReqArg set_to "NAME") "Set sample species to NAME" ]
  where
    set_output a c = return $ c { conf_maf_output = a }
    set_ref    a c = return $ c { conf_maf_reference = a }
    set_from   a c = return $ c { conf_maf_ref_species = fromString a }
    set_to     a c = return $ c { conf_maf_oth_species = fromString a }

conf_maf :: ConfMaf
conf_maf = ConfMaf (error "no output file specified")
                   (error "no reference specified")
                   (error "no reference species specified")
                   (error "no sample species specified")

main_maf :: [String] -> IO ()
main_maf args = do
    (maffs,ConfMaf{..}) <- parseFileOpts conf_maf "maf" "[maf-file...]" opts_maf args
    ref <- readTwoBit conf_maf_reference
    withFile conf_maf_output WriteMode $ \ohdl ->
        L.hPut ohdl . encodeGenome =<<
           foldM (\g f -> withInputFile f $
                            parseMaf (conf_maf_ref_species, conf_maf_oth_species)
                                     (rss_chroms ref) g . decomp . S.fromHandle)
                 emptyGenome (if null maffs then ["-"] else maffs)


data ConfPatch = ConfPatch {
    conf_patch_width      :: Int64,
    conf_patch_output     :: (Handle -> IO ()) -> IO (),
    conf_patch_reference  :: Either String FilePath,
    conf_patch_sample     :: FilePath }

defaultPatchConf :: ConfPatch
defaultPatchConf = ConfPatch 50 ($ stdout)
                             (Left    "no reference specified")
                             (error "no sample file specified")


opts_patch :: [ OptDescr (ConfPatch -> IO ConfPatch) ]
opts_patch =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hetfa)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "s" ["sample"] (ReqArg set_sample "FILE") "Read sample from FILE (.hef)"
    , Option "w" ["width"]  (ReqArg  set_width  "NUM") "Set width of FastA to NUM (50)" ]
  where
    set_output "-" c = return $ c { conf_patch_output    = ($ stdout) }
    set_output  a  c = return $ c { conf_patch_output    = withFile a WriteMode }
    set_ref     a  c = return $ c { conf_patch_reference = Right a }
    set_sample  a  c = return $ c { conf_patch_sample    =       a }
    set_width   a  c = (\n   -> c { conf_patch_width     =       n }) <$> readIO a

main_patch :: [String] -> IO ()
main_patch args = do
    ConfPatch{..} <- parseOpts defaultPatchConf (mk_opts "patch" [] opts_patch) args
    withFile conf_patch_sample ReadMode $ \hdli -> do
        (mref, raw) <- getRefPath $ decomp $ S.fromHandle hdli
        ref <- readTwoBit $ either (`fromMaybe` mref) id conf_patch_reference

        conf_patch_output $ \hdlo ->
            patchFasta hdlo conf_patch_width
                (rss_chroms ref) (rss_seqs ref) (decode (Right ref) raw)

main_dumplump :: [String] -> IO ()
main_dumplump [ref,inf] = do rs <- readTwoBit ref
                             withFile inf ReadMode $ debugLump . decode (Right rs) . decomp . S.fromHandle
main_dumplump [  inf  ] =    withFile inf ReadMode $ debugLump . decode (Left "no reference given") . decomp . S.fromHandle
main_dumplump     _     =    hPutStrLn stderr "Usage: dumplump [foo.hef]"


patchFasta :: Handle                    -- ^ output handle
           -> Int64                     -- ^ line width
           -> [B.ByteString]            -- ^ chromosome names
           -> [() -> RefSeq]            -- ^ reference sequences
           -> Stream (Of Lump) IO r     -- ^ input heffalump
           -> IO r
patchFasta hdl wd = p1
  where
    p1 [    ]     _  p = Q.effects p
    p1      _ [    ] p = Q.effects p
    p1 (c:cs) (r:rs) p = do hPutStrLn hdl $ '>' : B.unpack c
                            p2 (p1 cs rs) 0 (patch (r ()) p)

    p2 k l | l == wd = \f -> hPutChar hdl '\n' >> p2 k 0 f
    p2 k l = inspect >=> \case
        Left p            -> when (l>0) (L.hPutStrLn hdl L.empty) >> k p
        Right (Short c f) -> hPutChar hdl c >> p2 k (succ l) f
        Right (Long  s f)
            | L.null s    -> p2 k l f
            | otherwise   -> do let (u,v) = L.splitAt (wd-l) s
                                L.hPutStr hdl u
                                p2 k (l + L.length u) (wrap $ Long v f)

