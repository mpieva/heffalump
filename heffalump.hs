{-# LANGUAGE LambdaCase #-}
import BasePrelude
import Data.ByteString.Builder          ( hPutBuilder )
import Data.Fix
import Streaming
import System.Console.GetOpt
import System.FilePath                  ( takeBaseName )
import System.IO

import qualified Data.ByteString.Builder         as B
import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.IntMap                     as I
import qualified Data.Vector                     as V
import qualified Data.Vector.Unboxed             as U
import qualified Data.ByteString.Streaming as S

import Bamin
import Bcf
import Eigenstrat
import Emf
import Lump
import NewRef
import SillyStats
import Stretch ( main_dumppatch )
import Treemix
import Util
import VcfScan

main :: IO ()
main = do
    as0 <- getArgs
    case as0 of
        [] -> usage
        a:as -> case [ m | (k,(m,_)) <- mains, a `isPrefixOf` k ] of
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
        , z "het2fa"      main_hetfa'             "Import hetfa file (new)"
        , z "bamin"       main_bam                "Import BAM file"
        , z "maf"         main_maf                "Import two-species maf"
        , z "emf"         main_emf                "Import one genome from emf (Compara)"
        , z "patch"       main_patch              "Make a hetfa file by patching the reference"
        , z "treemix"     main_treemix            "Merge heffalumps into Treemix format"
        , z "kayvergence" main_kayvergence        "Compute Kayvergence ratios"
        , z "dstatistics" main_patterson          "Compute Patterson's D"
        , z "yaddayadda"  main_yaddayadda         "Compute Yadda-yadda-counts"
        , z "bcfhc"   (main_bcf "bcfhc" id)       "Import high coverage bcf (diploid)"
        , z "bcflc"   (main_bcf "bcflc" make_hap) "Import low coverage bcf (haploid)"
        , z "vcfhc"   (main_vcf "vcfhc" id)       "Import high coverage vcf (diploid)"
        , z "vcflc"   (main_vcf "vcflc" make_hap) "Import low coverage vcf (haploid)"
        , z "debhetfa"    main_debhetfa           "(debugging aid)"
        , z "debmaf"      main_debmaf             "(debugging aid)"
        , z "dumppatch"   main_dumppatch          "(debugging aid)"
        , z "dumplump"    main_dumplump           "(debugging aid)" ]


data ConfImportGen = ConfImportGen {
    conf_imp_width      :: Int64,           -- only for FastA output
    conf_imp_reference  :: FilePath,
    conf_imp_output     :: FilePath,
    conf_imp_sample     :: FilePath }
  deriving Show

defaultImportConf :: ConfImportGen
defaultImportConf = ConfImportGen 50 (error "no reference specified")
                                     (error "no output file specified")
                                     (error "no sample file specified")


data ConfMergeGen = ConfMergeGen {
    conf_noutgroups :: Int,
    conf_blocksize  :: Int,
    conf_all        :: Bool,
    conf_split      :: Bool,
    conf_reference  :: Maybe FilePath,
    conf_nrefpanel  :: Int,
    conf_output     :: FilePath,
    conf_sample     :: FilePath }
  deriving Show

defaultMergeConf :: ConfMergeGen
defaultMergeConf = ConfMergeGen 1 5000000 True True Nothing
                                (error "size of reference panel not known")
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
    smp  <- readSampleFa conf_imp_sample

    withFile conf_imp_output WriteMode $ \hdl ->
        hPutBuilder hdl $ importHetfa refs smp

main_hetfa' :: [String] -> IO ()
main_hetfa' args = do
    ( _, ConfImportGen{..} ) <- parseOpts False defaultImportConf (mk_opts "hetfa" [] opts_hetfa) args
    refs <- readTwoBit conf_imp_reference

    runResourceT $ S.writeFile conf_imp_output $
            importHetfa' refs $ readSampleFa' conf_imp_sample

-- Import Hetfa by diffing it against a reference.  We must read the
-- Hetfa in the order of the file (Fasta doesn't have an index), but
-- write it in the order of the reference.  So we parse a chromosome at
-- a time, sort in memory and dump it at the end.
--
-- Chromosomes in input that don't have match in the reference are
-- silently ignored.  Chromsomes that don't show up in the input are
-- emitted as empty sequences.  If no sequence matches, we produce a
-- warning message.
importHetfa :: NewRefSeqs -> [( B.ByteString, L.ByteString )] -> B.Builder
importHetfa ref smps
    | I.null map1 = error "Found only unexpected sequences.  Is this the right reference?"
    | otherwise   = encodeHeader ref <>
                    foldMap (\i -> maybe noLump B.lazyByteString  $ I.lookup i map1)
                            [0 .. length (nrss_chroms ref) - 1]
  where
    noLump = encodeLump $ Fix (Break (Fix Done))

    map1 :: I.IntMap L.ByteString
    map1 = foldl' enc1 I.empty smps

    enc1 :: I.IntMap L.ByteString -> (B.ByteString, L.ByteString) -> I.IntMap L.ByteString
    enc1 im (nm,sq) = maybe im (\i -> let !l = enc2 i sq in trace (show (i, L.length l)) (I.insert i l im))
                      $ findIndex (nm ==) (nrss_chroms ref)

    enc2 :: Int -> L.ByteString -> L.ByteString
    enc2 i sq = B.toLazyByteString . encodeLump . normalizeLump .
                diff2 (nrss_seqs ref !! i) sq $ Fix (Break (Fix Done))

importHetfa' :: Monad m => NewRefSeqs -> Stream (FastaSeq m) m r -> S.ByteString m r
importHetfa' ref smps = S.mwrap $ do -- ?!
    (map1,r) <- fold_1 I.empty smps

    when (I.null map1) $ error
            "Found only unexpected sequences.  Is this the right reference?"

    return $ do S.fromLazy $ B.toLazyByteString $ encodeHeader ref
                forM_ [0 .. length (nrss_chroms ref) - 1] $
                        \i -> S.fromLazy $ I.findWithDefault noLump i map1
                pure r
  where
    noLump = B.toLazyByteString $ encodeLump $ Fix (Break (Fix Done))

    enc2 :: Int -> S.ByteString IO r -> L.ByteString
    enc2 i sq = B.toLazyByteString .
                S.concatBuilders . encodeLump' . normalizeLump' $
                diff2' (nrss_seqs ref !! i) sq >>= yields . Break

    fold_1 :: Monad m => I.IntMap L.ByteString -> Stream (FastaSeq m) m r -> m (I.IntMap L.ByteString, r)
    fold_1 !acc s = inspect s >>= \case
        Left r -> return $ (acc,r)
        Right (FastaSeq nm sq) -> case findIndex (nm ==) (nrss_chroms ref) of
            Nothing -> S.effects sq >>= fold_1 acc
            Just  i -> do (!lump, s') <- enc2 i sq
                          fold_1 (I.insert i lump acc) s'


opts_maf :: [ OptDescr ( (FilePath,FilePath) -> IO (FilePath,FilePath) ) ]
opts_maf =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)" ]
  where
    set_output a (_,b) = return (a,b)
    set_ref    b (a,_) = return (a,b)

main_maf :: [String] -> IO ()
main_maf args = do
    ( maffs, (conf_output, conf_ref) ) <- parseOpts True (error "no output file specified")
                                             (mk_opts "maf" "[maf-file...]" opts_maf) args
    ref <- readTwoBit conf_ref
    withFile conf_output WriteMode $ \hdl ->
        hPutBuilder hdl . encode ref . concatLumps =<<
           mapM (\f -> parseMaf . decomp <$> L.readFile f) maffs


opts_patch :: [ OptDescr (ConfImportGen -> IO ConfImportGen) ]
opts_patch =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hetfa)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "s" ["sample"] (ReqArg set_sample "FILE") "Read sample from FILE (.hef)"
    , Option "w" ["width"]  (ReqArg  set_width  "NUM") "Set width of FastA to NUM (50)" ]
  where
    set_output a c = return $ c { conf_imp_output    = a }
    set_ref    a c = return $ c { conf_imp_reference = a }
    set_sample a c = return $ c { conf_imp_sample    = a }
    set_width  a c = (\n   -> c { conf_imp_width     = n }) <$> readIO a

main_patch :: [String] -> IO ()
main_patch args = do
    ( _, ConfImportGen{..} ) <- parseOpts False defaultImportConf (mk_opts "patch" [] opts_patch) args
    ref <- readTwoBit conf_imp_reference
    pat <- decode (Right ref) . decomp <$> L.readFile conf_imp_sample

    withFile conf_imp_output WriteMode $ \hdl ->
        patchFasta hdl conf_imp_width (nrss_chroms ref) (nrss_seqs ref) pat


main_debmaf :: [String] -> IO ()
main_debmaf [maff] = debugLump . ($ Fix Done) . parseMaf . decomp =<< L.readFile maff
main_debmaf      _ = hPutStrLn stderr $ "Usage: debmaf [ref.2bit] [smp.hetfa]"

-- This is broken.  Have to do it differently, maybe match the
-- incoming chromosomes against the reference, print the index, dump
-- them one at a time.  Clone some code from importHetfa.
main_debhetfa :: [String] -> IO ()
main_debhetfa args = case args of
    [reff, smpf] -> do ref <- readTwoBit reff
                       readSampleFa smpf >>= mapM_ (enc1 ref)
    _            -> hPutStrLn stderr $ "Usage: debhetfa [ref.2bit] [smp.hetfa]"
  where
    -- enc1 :: NewRefSeqs -> (B.ByteString, L.ByteString) -> L.ByteString
    enc1 ref (nm,sq) = {- L.hPut stdout . B.toLazyByteString . encodeLump .-}
                       debugLump . normalizeLump .
                       ($ Fix (Break (Fix Done))) .
                       maybe id (\i -> diff2 (nrss_seqs ref !! i) sq) $
                       findIndex (nm ==) (nrss_chroms ref)

main_dumplump :: [String] -> IO ()
main_dumplump [ref,inf] = do rs <- readTwoBit ref
                             debugLump . decode (Right rs) . decomp =<< L.readFile inf
main_dumplump [  inf  ] =    debugLump . decode (Left "no reference given") . decomp =<< L.readFile inf
main_dumplump     _     =    hPutStrLn stderr "Usage: dumplump [foo.hef]"

opts_eigen :: [ OptDescr (ConfMergeGen -> IO ConfMergeGen) ]
opts_eigen =
    [ Option "o" ["output"]     (ReqArg set_output "FILE") "Write output to FILE.geno and FILE.snp"
    , Option "r" ["reference"]     (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "n" ["numoutgroups"]  (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (1)"
    , Option "t" ["only-transversions"] (NoArg set_no_all) "Output only transversion sites"
    , Option "b" ["only-biallelic"]   (NoArg set_no_split) "Discard, don't split, polyallelic sites" ]
  where
    set_output a c = return $ c { conf_output     =      a }
    set_ref    a c = return $ c { conf_reference  = Just a }
    set_nout   a c = (\n ->   c { conf_noutgroups =      n }) <$> readIO a
    set_no_all   c = return $ c { conf_all        =  False }
    set_no_split c = return $ c { conf_split      =  False }

-- merge multiple files with the reference, write Eigenstrat format (geno & snp files)
main_eigenstrat :: [String] -> IO ()
main_eigenstrat args = do
    ( hefs, ConfMergeGen{..} ) <- parseOpts True defaultMergeConf (mk_opts "eigenstrat" "[hef-file...]" opts_eigen) args
    (refs, inps) <- decodeMany conf_reference hefs

    withFile (conf_output ++ ".snp") WriteMode $ \hsnp ->
        withFile (conf_output ++ ".geno") WriteMode $ \hgeno -> do
            let vars = either (const id) addRef refs $
                       bool singles_only concat conf_split $
                       mergeLumps conf_noutgroups inps
            forM_ vars $ \Variant{..} ->
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
    singles_only = foldr (\xs xss -> case xs of [x] -> x : xss ; _ -> xss) []


-- VCF output, this time going through 'Lump' instead of 'Stretch'
main_vcfout :: [String] -> IO ()
main_vcfout args = do
    ( hefs, ConfMergeGen{..} ) <- parseOpts True defaultMergeConf { conf_noutgroups = 0 }
                                            (mk_opts "vcfexport" "[hef-file...]" (tail opts_eigen)) args
    (refs, inps) <- decodeMany conf_reference hefs
    let chrs = either error (V.fromList . nrss_chroms) refs

    B.putStr $ "##fileformat=VCFv4.1\n" <>
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" <>
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" <>
               mconcat [ '\t' `B.cons` B.pack (takeBaseName f) | f <- hefs ] <>
               B.singleton '\n'

    let the_vars = addRef (either error id refs) $ concat $ mergeLumps conf_noutgroups inps

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


main_vcf :: String -> (Fix Lump -> Fix Lump) -> [String] -> IO ()
main_vcf = main_xcf "vcf" readVcf

main_bcf :: String -> (Fix Lump -> Fix Lump) -> [String] -> IO ()
main_bcf = main_xcf "bcf" readBcf

main_xcf :: String -> (FilePath -> IO [RawVariant]) -> String -> (Fix Lump -> Fix Lump) -> [String] -> IO ()
main_xcf ext reader key trans args = do
    ( vcfs, (conf_output, conf_ref) ) <- parseOpts True (error "no output file specified")
                                                   (mk_opts key ("["++ext++"-file...]") opts_maf) args
    ref <- readTwoBit conf_ref
    withFile conf_output WriteMode $ \hdl ->
        hPutBuilder hdl . encode ref
        . trans . importVcf (nrss_chroms ref)
        . progress conf_output . dedupVcf . cleanVcf
        . concat =<< mapM reader vcfs
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

-- Lowers all calls to haploid.
make_hap :: Fix Lump -> Fix Lump
make_hap = id

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


-- Our FastA files don't run exactly parallel (some samples are
-- truncated, it appears).  Therefore, we have to parse the FastA input
-- to some extent.  We shall also guard against messed up chromosome
-- order, just to be safe.  To this end, we pass the list of desired
-- chromosomes in.  The parser looks for the correct header, then
-- extracts the sequence (concatenated, no line feeds).  If the header
-- isn't found (this may happen after a while if the sorting is messed
-- up), a hard error is caused.  Unexpected sequences are skipped.  This
-- function doesn't care about the alphabet, so it would work for hetfa
-- and even protein FastA.  Input is ordered list of chromosomes and
-- lines of input.
--
-- Note that we turn small letters into capital letters in the
-- sequences.  Saves us from dealing with them later.

-- We also have two use cases:  dense files (more important), where
-- missing parts are considered "no call", and sparse files (less
-- important), where missing parts are considered "matches".

-- We start the coordinate at one(!), for VCF is one-based.
-- Pseudo-variants of the "no call" type must be filtered beforehand,
-- see 'cleanVcf'.
importVcf :: [ B.ByteString ] -> [ RawVariant ] -> Fix Lump
importVcf = (.) normalizeLump . go
  where
    go [    ] = const $ Fix $ Break $ Fix Done
    go (c:cs) = generic (hashChrom c) 1
      where
        generic !_    _ [         ] = Fix $ Break $ Fix $ Done
        generic !hs pos (var1:vars)
            | hs /= rv_chrom var1   = Fix $ Break $ go cs (var1:vars)

            -- long gap, creates Ns
            | rv_pos var1 > pos     = Fix $ Ns (rv_pos var1 - pos) $ generic hs (rv_pos var1) (var1:vars)

            -- positions must match now
            | rv_pos var1 == pos =
                        if isVar var1
                          then Fix $ get_var_code var1 $ generic hs (succ pos) vars
                          else Fix $ Eqs2 1 $ generic hs (succ pos) vars

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

-- | Remove indel variants, since we can't very well use them.  Also
-- removes "no call" variants, which are equivalent to missing entries.
cleanVcf :: [ RawVariant ] -> [ RawVariant ]
cleanVcf = filter $ \RawVariant{..} ->
    B.length rv_vars  == B.length (B.filter (== ',') rv_vars) * 2 + 1
    && rv_gt /= 0x0000 && rv_gt /= 0xff00

-- | Some idiot decided to output multiple records for the same position
-- into some VCF files.  If we hit that, we take the first.  (Einmal mit
-- Profis arbeiten!)
dedupVcf :: [ RawVariant ] -> [ RawVariant ]
dedupVcf (v1:v2:vs)
    | rv_chrom v1 == rv_chrom v2 && rv_pos v1 == rv_pos v2  =      dedupVcf (v1:vs)
    | otherwise                                             = v1 : dedupVcf (v2:vs)
dedupVcf [v1] = [v1]
dedupVcf [  ] = [  ]

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

