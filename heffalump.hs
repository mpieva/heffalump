{-# LANGUAGE LambdaCase #-}
import BasePrelude
import Paths_heffalump                  ( version )
import Streaming
import System.Console.GetOpt
import System.FilePath                  ( takeBaseName )
import System.IO

import qualified Data.ByteString.Builder         as B
import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.ByteString.Streaming       as S
import qualified Data.IntMap                     as I
import qualified Data.Vector                     as V
import qualified Data.Vector.Unboxed             as U

import Bamin
import Bed
import Eigenstrat
import Emf
import Lump
import NewRef
import SillyStats
import Stretch ( main_dumppatch )
import Treemix
import Util
import Vcf

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
        hPutStrLn stderr $ "Usage:"
        let kl = maximum $ map (length . fst) mains
        forM_ mains $ \(k,(_,h)) -> hPutStrLn stderr $
                "  "++pn++' ':k++replicate (kl-length k) ' '++"  -- "++h

    mains = let z a b c = (a,(b,c)) in
        [ z "eigenstrat"  main_eigenstrat         "Merge heffalumps into Eigenstrat format"
        , z "vcfexport"   main_vcfout             "Merge heffalumps into VCF"
        , z "hetfa"       main_hetfa              "Import hetfa file"
        , z "bamin"       main_bam                "Import BAM file"
        , z "maf"         main_maf                "Import two-species maf"
        , z "emf"         main_emf                "Import one genome from emf (Compara)"
        , z "patch"       main_patch              "Make a hetfa file by patching the reference"
        , z "treemix"     main_treemix            "Merge heffalumps into Treemix format"
        , z "kayvergence" main_kayvergence        "Compute Kayvergence ratios"
        , z "dstatistics" main_patterson          "Compute Patterson's D"
        , z "yaddayadda"  main_yaddayadda         "Compute Yadda-yadda-counts"
        , z "bcfin"      (main_xcf conf_bcf)      "Import BCF"
        , z "vcfin"      (main_xcf conf_vcf)      "Import VCF"
        , z "dumppatch"   main_dumppatch          "(debugging aid)"
        , z "dumplump"    main_dumplump           "(debugging aid)"
        , z "twobitinfo"  main_2bitinfo           "list reference sequences"
        , z "twobittofa"  main_2bittofa           "extract Fasta from 2bit"
        , z "fatotwobit"  main_fato2bit           "convert Fasta to 2bit"
        , z "--version"   main_version            "print version and exit" ]

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
    ( _, ConfImportGen{..} ) <- parseOpts False defaultImportConf (mk_opts "hetfa" [] opts_hetfa) args
    refs <- readTwoBit conf_imp_reference

    runResourceT $ S.writeFile conf_imp_output $
            importHetfa refs $ readSampleFa conf_imp_sample


-- | Import Hetfa by diffing it against a reference.  We must read the
-- Hetfa in the order of the file (Fasta doesn't have an index), but
-- write it in the order of the reference.  So we parse a chromosome at
-- a time, sort in memory and dump it at the end.
--
-- Chromosomes in input that don't have match in the reference are
-- silently ignored.  Chromsomes that don't show up in the input are
-- emitted as empty sequences.  Only f no sequence matches, we produce
-- an error message.
importHetfa :: MonadIO m => NewRefSeqs -> Stream (FastaSeq m) m r -> S.ByteString m r
importHetfa ref smps = S.mwrap $ do
    (map1,r) <- fold_1 I.empty smps

    when (I.null map1) $ error
            "Found only unexpected sequences.  Is this the right reference?"

    return $ do S.fromLazy $ B.toLazyByteString $ encodeHeader ref
                forM_ [0 .. length (nrss_chroms ref) - 1] $
                        \i -> S.fromLazy $ unpackLump $ I.findWithDefault noLump i map1
                pure r
  where
    enc2 :: MonadIO m => Int -> S.ByteString m r -> m (Of PackedLump r)
    enc2 i sq = encodeLumpToMem $ diff2' (nrss_seqs ref !! i) sq >>= yields . Break

    fold_1 :: MonadIO m => I.IntMap PackedLump -> Stream (FastaSeq m) m r -> m (I.IntMap PackedLump, r)
    fold_1 !acc s = inspect s >>= \case
        Left r -> return $ (acc,r)
        Right (FastaSeq nm sq) -> case findIndex (nm ==) (nrss_chroms ref) of
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
    (maffs,ConfMaf{..}) <- parseOpts True conf_maf (mk_opts "maf" "[maf-file...]" opts_maf) args
    ref <- readTwoBit conf_maf_reference
    withFile conf_maf_output WriteMode $ \hdl ->
        L.hPut hdl . encodeGenome =<<
           foldM (\g f -> parseMaf (conf_maf_ref_species, conf_maf_oth_species)
                                   (nrss_chroms ref) g . decomp =<< L.readFile f)
                 emptyGenome maffs


data ConfPatch = ConfPatch {
    conf_patch_width      :: Int64,
    conf_patch_reference  :: Either String FilePath,
    conf_patch_output     :: FilePath,
    conf_patch_sample     :: FilePath }
  deriving Show

defaultPatchConf :: ConfPatch
defaultPatchConf = ConfPatch 50 (Left    "no reference specified")
                                (error "no output file specified")
                                (error "no sample file specified")


opts_patch :: [ OptDescr (ConfPatch -> IO ConfPatch) ]
opts_patch =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hetfa)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "s" ["sample"] (ReqArg set_sample "FILE") "Read sample from FILE (.hef)"
    , Option "w" ["width"]  (ReqArg  set_width  "NUM") "Set width of FastA to NUM (50)" ]
  where
    set_output a c = return $ c { conf_patch_output    =       a }
    set_ref    a c = return $ c { conf_patch_reference = Right a }
    set_sample a c = return $ c { conf_patch_sample    =       a }
    set_width  a c = (\n   -> c { conf_patch_width     =       n }) <$> readIO a

main_patch :: [String] -> IO ()
main_patch args = do
    ( _, ConfPatch{..} ) <- parseOpts False defaultPatchConf (mk_opts "patch" [] opts_patch) args
    raw <- decomp <$> L.readFile conf_patch_sample
    ref <- readTwoBit $ either (\e -> fromMaybe e $ getRefPath raw) id conf_patch_reference

    withFile conf_patch_output WriteMode $ \hdl ->
        patchFasta hdl conf_patch_width (nrss_chroms ref)
                   (nrss_seqs ref) (decode (Right ref) raw)

main_dumplump :: [String] -> IO ()
main_dumplump [ref,inf] = do rs <- readTwoBit ref
                             debugLump . decode (Right rs) . decomp =<< L.readFile inf
main_dumplump [  inf  ] =    debugLump . decode (Left "no reference given") . decomp =<< L.readFile inf
main_dumplump     _     =    hPutStrLn stderr "Usage: dumplump [foo.hef]"

-- VCF output, this time going through 'Lump' instead of 'Stretch'
main_vcfout :: [String] -> IO ()
main_vcfout args = do
    ( hefs, ConfMergeGen{..} ) <- parseOpts True defaultMergeConf { conf_noutgroups = 0 }
                                            (mk_opts "vcfexport" "[hef-file...]" (tail opts_eigen)) args
    (refs, inps) <- decodeMany conf_reference hefs
    region_filter <- mkBedFilter conf_regions (either error nrss_chroms refs)
    let chrs = either error (V.fromList . nrss_chroms) refs

    B.putStr $ "##fileformat=VCFv4.1\n" <>
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" <>
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" <>
               mconcat [ '\t' `B.cons` B.pack (takeBaseName f) | f <- hefs ] <>
               B.singleton '\n'

    let the_vars = addRef (either error id refs) $ region_filter $ concat $ mergeLumps conf_noutgroups inps

    forM_ the_vars $ \Variant{..} ->
        -- samples (not outgroups) must show alt allele at least once
        let ve    = U.foldl' (.|.) 0 $ U.drop conf_noutgroups v_calls
            is_ti = conf_all || isTransversion v_alt in
        when (ve .&. 12 /= 0 && is_ti) $ do
            B.hPutBuilder stdout $
                B.byteString (chrs V.! v_chr) <> B.char8 '\t' <>
                B.intDec (v_pos+1) <> B.string8 "\t.\t" <>
                B.char8 (toRefCode v_ref) <> B.char8 '\t' <>
                B.char8 (toAltCode v_alt v_ref) <> B.string8 "\t.\t.\t.\tGT" <>
                U.foldr ((<>) . B.byteString . (V.!) gts . fromIntegral) mempty v_calls <>
                B.char8 '\n'
  where
    gts :: V.Vector B.ByteString
    gts = V.fromList [ "\t./."      -- 0, N
                     , "\t0"        -- 1, 1xref
                     , "\t0/0"      -- 2, 2xref
                     , "\t0/0"      -- 3, 3xref (whatever)
                     , "\t1"        -- 4, 1xalt
                     , "\t0/1"      -- 5, ref+alt
                     , "\t0/1"      -- 6, 2xref+alt (whatever)
                     , "\t0/1"      -- 7, 3xref+alt (whatever)
                     , "\t1/1"      -- 8, 2xalt
                     , "\t0/1"      -- 9, ref+2xalt (whatever)
                     , "\t0/1"      -- 10, 2xref+2xalt (whatever)
                     , "\t0/1"      -- 11, 3xref+2xalt (whatever)
                     , "\t1/1"      -- 12, 3xalt (whatever)
                     , "\t0/1"      -- 13, ref+3xalt (whatever)
                     , "\t0/1"      -- 14, 2xref+3xalt (whatever)
                     , "\t0/1" ]    -- 15, 3xref+3xalt (whatever)


patchFasta :: Handle -> Int64 -> [B.ByteString] -> [() -> NewRefSeq] -> Fix Lump -> IO ()
patchFasta hdl wd = p1
  where
    p1 :: [B.ByteString] -> [() -> NewRefSeq] -> Fix Lump -> IO ()
    p1 [    ]     _  _ = return ()
    p1      _ [    ] _ = return ()
    p1 (c:cs) (r:rs) p = do hPutStrLn hdl $ '>' : B.unpack c
                            p2 (p1 cs rs) 0 (patch (r ()) p)

    p2 :: (Fix Lump -> IO ()) -> Int64 -> Frag -> IO ()
    p2 k l f | l == wd = L.hPutStrLn hdl L.empty >> p2 k 0 f
    p2 k l (Term p)    = when (l>0) (L.hPutStrLn hdl L.empty) >> k p
    p2 k l (Short c f) = hPutChar hdl c >> p2 k (succ l) f
    p2 k l (Long  s f) = case L.splitAt (wd-l) s of
            _ | L.null s -> p2 k l f
            (u,v)        -> L.hPutStr hdl u >> p2 k (l + L.length u) (Long v f)

