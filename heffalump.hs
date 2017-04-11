-- TODO:
-- - switch everything to the 'Lump' representation (and acknowledge
--   that the 'Stretch' was a bad choice)
-- - allow arbitrary reference genomes (use 2bit everywhere)
-- - import from multi-MAF (either create multiple heffalumps or allow
--   selection of one species)
-- - allow import from unsorted MAF or EMF (make 'Lump's from each block
--   in the file, reorder them in memory)


import BasePrelude
import Data.ByteString.Builder          ( hPutBuilder, Builder )
import System.Console.GetOpt
import System.FilePath                  ( takeBaseName )
import System.IO

import qualified Data.ByteString.Builder         as B
import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.Vector                     as V
import qualified Data.Vector.Unboxed             as U

import Bamin
import Bcf
import Eigenstrat
import Emf
import qualified Lump
import Merge
import NewRef ( Nuc2b(..), Var2b(..), addRef, readTwoBit )
import qualified NewRef
import SillyStats
import Stretch
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
        , z "vcfexport"   main_vcfout             "Merge heffalumps into VCF (Stretch)"
        , z "vcflexport"  main_vcfout_2           "Merge heffalumps into VCF (Lump)"
        , z "hetfa"       main_hetfa              "Import hetfa file"
        , z "bamin"       main_bam                "Import BAM file"
        , z "maf"         main_maf                "Import two-species maf"
        , z "emf"         main_emf                "Import one genome from emf (Compara)"
        , z "patch"       main_patch              "Make a hetfa file by patching the reference"
        , z "treemix"     main_treemix            "Merge heffalumps into Treemix format"
        , z "kayvergence" main_kayvergence        "Compute Kayvergence ratios"
        , z "dstatistics" main_patterson          "Compute Patterson's D"
        , z "yaddayadda"  main_yaddayadda         "Compute Yadda-yadda-counts"
        , z "bcfhc" (main_bcf "bcfhc" encode_dip) "Import high coverage bcf (diploid)"
        , z "bcflc" (main_bcf "bcflc" encode_hap) "Import low coverage bcf (haploid)"
        , z "vcfhc" (main_vcf "vcfhc" encode_dip) "Import high coverage vcf (diploid)"
        , z "vcflc" (main_vcf "vcflc" encode_hap) "Import low coverage vcf (haploid)"
        , z "cdstatistics" main_patterson_corr    "(failed experiment)"
        , z "debhetfa"    main_debhetfa           "(debugging aid)"
        , z "debmaf"      main_debmaf             "(debugging aid)"
        , z "dumppatch"   main_dumppatch          "(debugging aid)"
        , z "test"        main_test               "(test)" ]


data ConfGen = ConfGen {
    conf_noutgroups :: Int,
    conf_blocksize  :: Int,
    conf_width      :: Int64,
    conf_all        :: Bool,
    conf_nosplit    :: Bool,
    conf_nrefpanel  :: Int,
    conf_output     :: FilePath,
    conf_reference  :: FilePath,
    conf_sample     :: FilePath }
  deriving Show

defaultConfGen :: ConfGen
defaultConfGen = ConfGen 1 5000000 50 True False
                         (error "size of reference panel not known")
                         (error "no output file specified")
                         (error "no reference file specified")
                         (error "no sample file specified")


opts_hetfa :: [ OptDescr (ConfGen -> IO ConfGen) ]
opts_hetfa =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "s" ["sample"] (ReqArg set_sample "FILE") "Read sample from FILE (.hetfa)" ]
  where
    set_output a c = return $ c { conf_output    = a }
    set_ref    a c = return $ c { conf_reference = a }
    set_sample a c = return $ c { conf_sample    = a }


main_hetfa :: [String] -> IO ()
main_hetfa args = do
    ( _, ConfGen{..} ) <- parseOpts False defaultConfGen (mk_opts "hetfa" [] opts_hetfa) args
    -- (_, Reference ref) <- readReference conf_reference
    refs <- map (\(_,_,s) -> s) <$> readTwoBit conf_reference
    smp <- readSampleFa conf_sample

    withFile conf_output WriteMode $ \hdl ->
        hPutBuilder hdl . encode_dip . catStretches $ zipWith diff2 refs smp


opts_maf :: [ OptDescr ( FilePath -> IO FilePath ) ]
opts_maf =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hef)" ]
  where
    set_output a _ = return a

main_maf :: [String] -> IO ()
main_maf args = do
    ( maffs, conf_output ) <- parseOpts True (error "no output file specified")
                                             (mk_opts "maf" "[maf-file...]" opts_maf) args
    withFile conf_output WriteMode $ \hdl ->
        hPutBuilder hdl . encode_hap . catStretches =<<
           mapM (\f -> parseMaf . decomp <$> L.readFile f) maffs


opts_patch :: [ OptDescr (ConfGen -> IO ConfGen) ]
opts_patch =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hetfa)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "s" ["sample"] (ReqArg set_sample "FILE") "Read sample from FILE (.hef)"
    , Option "w" ["width"]  (ReqArg  set_width  "NUM") "Set width of FastA to NUM (50)" ]
  where
    set_output a c = return $ c { conf_output    = a }
    set_ref    a c = return $ c { conf_reference = a }
    set_sample a c = return $ c { conf_sample    = a }
    set_width  a c = (\n   -> c { conf_width     = n }) <$> readIO a

main_patch :: [String] -> IO ()
main_patch args = do
    ( _, ConfGen{..} ) <- parseOpts False defaultConfGen (mk_opts "patch" [] opts_patch) args
    (_,ref) <- readReference conf_reference
    pat <- decode . decomp <$> L.readFile conf_sample

    withFile conf_output WriteMode $ \hdl ->
            patchFasta hdl conf_width chroms ref pat


main_debmaf :: [String] -> IO ()
main_debmaf [maff] = debugStretch . ($ Done) . parseMaf . decomp =<< L.readFile maff
main_debmaf _ =  hPutStrLn stderr $ "Usage: debmaf [ref.fa] [smp.hetfa]"

main_debhetfa :: [String] -> IO ()
main_debhetfa [reff, smpf] = do
    (_,Reference ref) <- readReference reff
    smp <- readSampleFa smpf
    debugStretch . ($ Done) . head $ zipWith diff ref smp
main_debhetfa _ = hPutStrLn stderr $ "Usage: debhetfa [ref.fa] [smp.hetfa]"

main_dumppatch :: [String] -> IO ()
main_dumppatch [inf] = debugStretch . decode . decomp =<< L.readFile inf
main_dumppatch     _ = hPutStrLn stderr "Usage: dumppatch [foo.hef]"


opts_eigen :: [ OptDescr (ConfGen -> IO ConfGen) ]
opts_eigen =
    [ Option "o" ["output"]     (ReqArg set_output "FILE") "Write output to FILE.geno and FILE.snp"
    , Option "r" ["reference"]     (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "n" ["numoutgroups"]  (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (1)"
    , Option "t" ["only-transversions"] (NoArg set_no_all) "Output only transversion sites"
    , Option "b" ["only-biallelic"]   (NoArg set_no_split) "Discard, don't split, polyallelic sites" ]
  where
    set_output a c = return $ c { conf_output    = a }
    set_ref    a c = return $ c { conf_reference = a }
    set_nout   a c = readIO a >>= \n -> return $ c { conf_noutgroups = n }
    set_no_all   c = return $ c { conf_all = False }
    set_no_split c = return $ c { conf_nosplit = True }

-- merge multiple files with the reference, write Eigenstrat format (geno & snp files)
main_eigenstrat :: [String] -> IO ()
main_eigenstrat args = do
    ( hefs, ConfGen{..} ) <- parseOpts True defaultConfGen (mk_opts "eigenstrat" "[hef-file...]" opts_eigen) args
    (_,ref) <- readReference conf_reference
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs

    withFile (conf_output ++ ".snp") WriteMode $ \hsnp ->
        withFile (conf_output ++ ".geno") WriteMode $ \hgeno ->
            forM_ (merge_hefs conf_nosplit conf_noutgroups ref inps) $ \Variant{..} ->
                -- samples (not outgroups) must show ref and alt allele at least once
                let ve = U.foldl' (.|.) 0 $ U.drop conf_noutgroups v_calls
                    is_ti = conf_all || is_transversion (toUpper v_ref) (toUpper v_alt) in
                when (ve .&. 3 /= 0 && ve .&. 12 /= 0 && is_ti) $ do
                    hPutStrLn hgeno $ map (B.index "9222011101110111" . fromIntegral) $ U.toList v_calls
                    hPutStrLn hsnp $ intercalate "\t"
                        -- 1st column is SNP name
                        [ mkname v_chr v_pos v_alt
                        -- "2nd column is chromosome.  X chromosome is encoded as 23.
                        -- Also, Y is encoded as 24, mtDNA is encoded as 90, ..."
                        , show $ if v_chr == 24 then 90 else v_chr + 1
                        -- "3rd column is genetic position (in Morgans).
                        -- If unknown, ok to set to 0.0"
                        , "0.0"
                        -- "4th column is physical position (in bases)"
                        , show (v_pos+1)
                        -- "Optional 5th and 6th columns are reference and variant alleles"
                        , [toUpper v_ref], [toUpper v_alt] ]
  where
    is_transversion 'C' 'T' = False
    is_transversion 'T' 'C' = False
    is_transversion 'G' 'A' = False
    is_transversion 'A' 'G' = False
    is_transversion  _   _  = True

main_vcfout :: [String] -> IO ()
main_vcfout args = do
    ( hefs, ConfGen{..} ) <- parseOpts True defaultConfGen { conf_noutgroups = 0 }
                                            (mk_opts "vcfexport" "[hef-file...]" (tail opts_eigen)) args
    (chrs,ref) <- readReference conf_reference
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs

    B.putStr $ "##fileformat=VCFv4.1\n" <>
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" <>
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" <>
               mconcat [ '\t' `B.cons` B.pack (takeBaseName f) | f <- hefs ] <>
               B.singleton '\n'

    forM_ (merge_hefs False conf_noutgroups ref inps) $ \Variant{..} ->
        -- samples (not outgroups) must show alt allele at least once
        let ve    = U.foldl' (.|.) 0 $ U.drop conf_noutgroups v_calls
            is_ti = conf_all || is_transversion (toUpper v_ref) (toUpper v_alt) in
        when (ve .&. 12 /= 0 && is_ti) $ do
            B.hPutBuilder stdout $
                B.byteString (chrs !! v_chr) <> B.char8 '\t' <>
                B.intDec (v_pos+1) <> B.string8 "\t.\t" <>
                B.char8 (toUpper v_ref) <> B.char8 '\t' <>
                B.char8 (toUpper v_alt) <> B.string8 "\t.\t.\t.\tGT" <>
                U.foldr ((<>) . B.byteString . (V.!) gts . fromIntegral) mempty v_calls <>
                B.char8 '\n'
  where
    is_transversion 'C' 'T' = False
    is_transversion 'T' 'C' = False
    is_transversion 'G' 'A' = False
    is_transversion 'A' 'G' = False
    is_transversion  _   _  = True

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



-- VCF output, this time going through 'Lump' instead of 'Stretch'
main_vcfout_2 :: [String] -> IO ()
main_vcfout_2 args = do
    ( hefs, ConfGen{..} ) <- parseOpts True defaultConfGen { conf_noutgroups = 0 }
                                            (mk_opts "vcfexport" "[hef-file...]" (tail opts_eigen)) args
    refs <- readTwoBit conf_reference
    let chrs = V.fromList [ c | (c,_,_) <- refs ]
    inps <- V.fromList <$> mapM (fmap (Lump.decode refs . decomp) . L.readFile) hefs

    B.putStr $ "##fileformat=VCFv4.1\n" <>
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" <>
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" <>
               mconcat [ '\t' `B.cons` B.pack (takeBaseName f) | f <- hefs ] <>
               B.singleton '\n'

    let the_vars = addRef refs $ concat $ Lump.merge_lumps conf_noutgroups inps

    forM_ the_vars $ \NewRef.Variant{..} ->
        -- samples (not outgroups) must show alt allele at least once
        let ve    = U.foldl' (.|.) 0 $ U.drop conf_noutgroups v_calls
            is_ti = conf_all || is_transversion v_alt in
        when (ve .&. 12 /= 0 && is_ti) $ do
            B.hPutBuilder stdout $
                B.byteString (chrs V.! v_chr) <> B.char8 '\t' <>
                B.intDec (v_pos+1) <> B.string8 "\t.\t" <>
                B.char8 (toRefCode v_ref) <> B.char8 '\t' <>
                B.char8 (toAltCode v_alt v_ref) <> B.string8 "\t.\t.\t.\tGT" <>
                U.foldr ((<>) . B.byteString . (V.!) gts . fromIntegral) mempty v_calls <>
                B.char8 '\n'
  where
    toAltCode (V2b v) (N2b r) = B.index "TCAGXPOI" $ fromIntegral (xor r v)
    toRefCode = toAltCode (V2b 0)

    is_transversion (V2b v) = testBit v 1

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


main_vcf :: String -> (Stretch -> Builder) -> [String] -> IO ()
main_vcf = main_xcf "vcf" readVcf

main_bcf :: String -> (Stretch -> Builder) -> [String] -> IO ()
main_bcf = main_xcf "bcf" readBcf

main_xcf :: String -> (FilePath -> IO [RawVariant]) -> String
         -> (Stretch -> Builder) -> [String] -> IO ()
main_xcf ext reader key enc args = do
    ( vcfs, conf_output ) <- parseOpts True (error "no output file specified")
                                            (mk_opts key ("["++ext++"-file...]") opts_maf) args
    withFile conf_output WriteMode $ \hdl ->
        hPutBuilder hdl . enc .
            importVcf chroms . progress conf_output . dedupVcf . cleanVcf . concat
                =<< mapM reader vcfs
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


patchFasta :: Handle -> Int64 -> [L.ByteString] -> Reference -> Stretch -> IO ()
patchFasta hdl wd = p1
  where
    p1 :: [L.ByteString] -> Reference -> Stretch -> IO ()
    p1 [    ]                 _  _ = return ()
    p1      _ (Reference [    ]) _ = return ()
    p1 (c:cs) (Reference (r:rs)) p = do
        hPutStrLn hdl $ '>' : L.unpack c
        p2 (p1 cs (Reference rs)) 0 (patch r p)

    p2 :: (Stretch -> IO ()) -> Int64 -> Frag -> IO ()
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


-- Reads and parses a VCF.  We assume one VCF per chromosome, so we go
-- ahead and ignore the chromosome completely.  The rest has to be
-- turned into a 'Stretch'.
--
-- Two options...  either convert to two sequences and pass to diff, or
-- scan the list of variants and produce the stretch directly.
--
-- We also have two use cases:  dense files (more important), where
-- missing parts are considered "no call", and sparse files (less
-- important), where missing parts are considered "matches".

-- We start the coordinate at one(!), for VCF is one-based.
-- Pseudo-variants of the "no call" type must be filtered beforehand,
-- see 'cleanVcf'.
importVcf :: [ L.ByteString ] -> [ RawVariant ] -> Stretch
importVcf [    ] = const $ Break Done
importVcf (c:cs) = generic (hashChrom $ B.concat $ L.toChunks c) 1
  where
    generic !_    _ [         ] = Break Done
    generic !hs pos (var1:vars)
        | hs /= rv_chrom var1 = Break $ importVcf cs (var1:vars)

        -- long gap, creates Ns
        | rv_pos var1 >= pos + 2 = let l = (rv_pos var1 - pos) `div` 2
                                   in Ns (fromIntegral l) $ generic hs (pos + 2*l) (var1:vars)

        -- small gap, have to emit codes
        | rv_pos var1 == pos + 1 = Chrs (NucCode 0) (get_nuc_code var1) $ generic hs (pos+2) vars

        -- positions must match now
        | rv_pos var1 == pos = case vars of
                -- two variant calls next to each other
                -- if both are reference, we can try and build a stretch
                var2:vars'' | rv_pos var2 == pos+1 ->
                    if isVar var1 || isVar var2
                      then Chrs (get_nuc_code var1) (get_nuc_code var2) $ generic hs (pos+2) vars''
                      else matches hs 1 (pos+2) vars''

                -- one followed by gap, will become an N
                _ -> Chrs (get_nuc_code var1) (NucCode 0) $ generic hs (pos+2) vars

        | otherwise = error $ "Got variant position " ++ show (rv_pos var1)
                           ++ " when expecting " ++ show pos ++ " or higher."

    -- To extend a stretch of matches, we need
    -- two non-vars at the next two positions
    matches hs num pos (var1:var2:vars)
        | rv_pos var1 == pos && rv_pos var2 == pos+1 && not (isVar var1) && not (isVar var2)
            = matches hs (num+1) (pos+2) vars

    -- anything else, we dump the stretch out and pass the buck
    matches hs num pos vars = Eqs num $ generic hs pos vars

    -- *sigh*  Have to turn a numeric genotype into a 'NucCode'.  We
    -- have characters for the variants, and we need to map a pair of
    -- them to a code.
    get_nuc_code RawVariant{..}
        -- haploid call or one allele missing
        | rv_gt .&. 0xFF00 == 0xFF00 || rv_gt .&. 0xFF00 == 0x0000
            = let z = toUpper $ safeIndex "z" rv_vars ( fromIntegral (rv_gt .&. 0x00FE - 2) )
              in NucCode $ maybe (error $ "What's a " ++ shows z "?") fromIntegral
                         $ B.elemIndex z nuc_chars
        | otherwise                         -- diploid
            = let c1 = toUpper $ safeIndex "c1" rv_vars ( fromIntegral (rv_gt            .&. 0x00FE - 2) )
                  c2 = toUpper $ safeIndex ("c2 "++shows rv_pos " " ++ showHex rv_gt " ") rv_vars ( fromIntegral (rv_gt `shiftR` 8 .&. 0x00FE - 2) )

                  n1 = maybe (error $ "What's a " ++ shows c1 "?") id $ B.elemIndex c1 nuc_chars
                  n2 = maybe (error $ "What's a " ++ shows c2 "?") id $ B.elemIndex c2 nuc_chars
              in NucCode $ two_to_code U.! (n1+5*n2)

    two_to_code = U.fromList [0,1,2,3,4,1,1,5,6,7,2,5,2,8,9,3,6,8,3,10,4,7,9,10,4]

    safeIndex m s i | B.length s > i = B.index s i
                    | otherwise = error $ "Attempted to index " ++ shows i " in " ++ shows s " (" ++ m ++ ")."

    nuc_chars :: B.ByteString
    nuc_chars = "NACGT"

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

main_test :: [String] -> IO ()
main_test hefs = do
    (chrs,ref) <- readReference "/var/tmp/ref.fa"
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs
    forM_ (merge_hefs' False 0 ref inps) $ \vs ->
        -- We need an ABBA+ pattern.  So all samples except the second
        -- and third must be alt alleles only while the second and third
        -- are ref only, or the other way around.  Any of the set of
        -- variants is fine, we print all.
        when (any is_good vs) $
            forM_ vs $ \Variant{..} ->
                B.hPutBuilder stdout $
                    B.byteString (chrs !! v_chr) <> B.char8 '\t' <>
                    B.intDec (v_pos+1) <> B.string8 "\t.\t" <>
                    B.char8 (toUpper v_ref) <> B.char8 '\t' <>
                    B.char8 (toUpper v_alt) <> B.string8 "\t.\t.\t.\tGT" <>
                    U.foldr ((<>) . B.byteString . (V.!) gts . fromIntegral) mempty v_calls <>
                    B.char8 '\n'
  where
    is_good Variant{..} = is_ti && defnd && (adda || daad)
      where
        c1    = v_calls U.! 0
        c2    = v_calls U.! 1
        c3    = v_calls U.! 2
        cs    = U.drop 3 v_calls
        is_ti = is_transversion (toUpper v_ref) (toUpper v_alt)

        adda  = U.all (not . isRef) cs && isRef c2 && isRef c3 && isAlt c1
        daad  = U.all (not . isAlt) cs && isAlt c2 && isAlt c3 && isRef c1
        defnd = U.length (U.filter (== 0) cs) <= 3

    is_transversion 'C' 'T' = False
    is_transversion 'T' 'C' = False
    is_transversion 'G' 'A' = False
    is_transversion 'A' 'G' = False
    is_transversion  _   _  = True

    isRef v = v /= 0 && v .&. 12 == 0
    isAlt v = v /= 0 && v .&.  3 == 0

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


