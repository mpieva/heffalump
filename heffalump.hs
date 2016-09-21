import BasePrelude
import Data.ByteString.Builder          ( hPutBuilder, Builder, intDec, char7 )
import System.Console.GetOpt
import System.Directory
import System.FilePath
import System.IO

import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Lazy            as LB
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.HashMap.Strict             as H
import qualified Data.Vector.Unboxed             as U

import Bamin
import Eigenstrat
import Merge
import SillyStats
import Stretch
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
        , z "hetfa"       main_hetfa              "Import hetfa file"
        , z "bamin"       main_bam                "Import BAM file"
        , z "maf"         main_maf                "Import two-species maf"
        , z "patch"       main_patch              "Make a hetfa file by patching the reference"
        , z "treemix"     main_treemix            "Merge heffalumps into Treemix format"
        , z "kayvergence" main_kayvergence        "Compute Kayvergence ratios"
        , z "dstatistics" main_patterson          "Compute Patterson's D"
        , z "yaddayadda"  main_yaddayadda         "Compute Yadda-yadda-counts"
        , z "vcfhc" (main_vcf "vcfhc" encode_dip) "Import high coverage vcf (diploid)"
        , z "vcflc" (main_vcf "vcflc" encode_hap) "Import low coverage vcf (haploid)"
        , z "debhetfa"    main_debhetfa           "(debugging aid)"
        , z "debmaf"      main_debmaf             "(debugging aid)"
        , z "dumppatch"   main_dumppatch          "(debugging aid)" ]


data ConfGen = ConfGen {
    conf_noutgroups :: Int,
    conf_blocksize  :: Int,
    conf_width      :: Int64,
    conf_all        :: Bool,
    conf_nrefpanel  :: Int,
    conf_output     :: FilePath,
    conf_reference  :: FilePath,
    conf_sample     :: FilePath }
  deriving Show

defaultConfGen :: ConfGen
defaultConfGen = ConfGen 1 5000000 50 True
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
    Reference ref <- readReference conf_reference
    smp <- readSampleFa conf_sample

    withFile conf_output WriteMode $ \hdl ->
        hPutBuilder hdl . encode_dip . catStretches $ zipWith diff ref smp


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
    ref <- readReference conf_reference
    pat <- decode . decomp <$> L.readFile conf_sample

    withFile conf_output WriteMode $ \hdl ->
            patchFasta hdl conf_width chroms ref pat


main_debmaf :: [String] -> IO ()
main_debmaf [maff] = debugStretch . ($ Done) . parseMaf . decomp =<< L.readFile maff
main_debmaf _ =  hPutStrLn stderr $ "Usage: debmaf [ref.fa] [smp.hetfa]"

main_debhetfa :: [String] -> IO ()
main_debhetfa [reff, smpf] = do
    Reference ref <- readReference reff
    smp <- readSampleFa smpf
    debugStretch . ($ Done) . head $ zipWith diff ref smp
main_debhetfa _ = hPutStrLn stderr $ "Usage: debhetfa [ref.fa] [smp.hetfa]"

main_dumppatch :: [String] -> IO ()
main_dumppatch [inf] = debugStretch . decode . decomp =<< L.readFile inf
main_dumppatch     _ = hPutStrLn stderr "Usage: dumppatch [foo.hef]"


opts_treemix :: [ OptDescr (ConfGen -> IO ConfGen) ]
opts_treemix =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.tmx)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "i" ["individuals"] (ReqArg set_indiv "FILE") "Read individuals from FILE (.ind)"
    , Option "n" ["numoutgroups"] (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (1)" ]
  where
    set_output a c = return $ c { conf_output    = a }
    set_ref    a c = return $ c { conf_reference = a }
    set_indiv  a c = return $ c { conf_sample    = a }
    set_nout   a c = readIO a >>= \n -> return $ c { conf_noutgroups = n }

main_treemix :: [String] -> IO ()
main_treemix args = do
    ( _, ConfGen{..} ) <- parseOpts False defaultConfGen (mk_opts "treemix" [] opts_treemix) args
    treemix conf_noutgroups conf_output conf_sample conf_reference


opts_eigen :: [ OptDescr (ConfGen -> IO ConfGen) ]
opts_eigen =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE.geno and FILE.snp"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "n" ["numoutgroups"] (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (1)"
    , Option "t" ["only-transversions"] (NoArg set_no_all) "Output only transversion sites" ]
  where
    set_output a c = return $ c { conf_output    = a }
    set_ref    a c = return $ c { conf_reference = a }
    set_nout   a c = readIO a >>= \n -> return $ c { conf_noutgroups = n }
    set_no_all   c = return $ c { conf_all = False }

-- merge multiple files with the reference, write Eigenstrat format (geno & snp files)
main_eigenstrat :: [String] -> IO ()
main_eigenstrat args = do
    ( hefs, ConfGen{..} ) <- parseOpts True defaultConfGen (mk_opts "eigenstrat" "[hef-file...]" opts_eigen) args
    ref <- readReference conf_reference
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs

    withFile (conf_output ++ ".snp") WriteMode $ \hsnp ->
        withFile (conf_output ++ ".geno") WriteMode $ \hgeno ->
            forM_ (merge_hefs conf_noutgroups ref inps) $ \Variant{..} ->
                -- samples (not outgroups) must show ref and alt allele at least once
                let ve = U.foldl' (.|.) 0 $ U.drop conf_noutgroups v_calls
                    is_ti = conf_all || is_transversion (toUpper v_ref) (toUpper v_alt) in
                when (ve .&. 3 /= 0 && ve .&. 12 /= 0 && is_ti) $ do
                    hPutStrLn hgeno $ map (B.index "9222011101110111" . fromIntegral) $ U.toList v_calls
                    hPutStrLn hsnp $ intercalate "\t"
                        -- 1st column is SNP name
                        [ mkname v_chr v_pos
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



main_vcf :: String -> (Stretch -> Builder) -> [String] -> IO ()
main_vcf key enc args = do
    ( vcfs, conf_output ) <- parseOpts True (error "no output file specified")
                                            (mk_opts key "[vcf-file...]" opts_maf) args
    withFile conf_output WriteMode $ \hdl ->
        hPutBuilder hdl . enc .
            importVcf chroms . dedupVcf . cleanVcf . concat
                =<< mapM readVcf vcfs


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

-- Read MAF.  We discard the header lines starting with '#'.  The rest
-- must be alignment blocks.  Each a-block must contain its own header
-- line and exactly two s-lines.  The first line must have human
-- coordinates, and the file must be sorted by human coordinates.   We
-- ignore the chromosome name.
--
-- We wish to skip about between blocks.  But sometimes blocks join
-- immediately, and then we get into trouble because 'diff' only works
-- correctly for stretches of even length.  So instead we return a new
-- data type to make that explicit.

data MafStretch = Skip !Int MafStretch | Aln !L.ByteString !L.ByteString MafStretch | MafDone

parseMaf :: L.ByteString -> Stretch -> Stretch
parseMaf inp done = norm . parse1 0 $ L.lines inp
  where
    parse1 oldpos lns = case dropWhile (\l -> L.null l || L.all isSpace l || L.head l == '#') lns of
        [] -> MafDone
        (ahdr:shum:sother:lns')
            | "a " `L.isPrefixOf` ahdr
            , ["s", _, hpos_, _hlen, "+",  _clen, hum_seq] <- L.words shum
            , ["s", _, _opos, _olen,   _, _oclen, oth_seq] <- L.words sother
            , Just (hpos, "") <- L.readInt hpos_
            , hpos >= oldpos ->
                -- We have to deal with gaps.  Gaps in the reference can
                -- simply be removed.  If the sample is gapped, we want
                -- to encode an N.  Thankfully, that's precisely what
                -- diff already does when it sees a gap.
                let (ref', smp') = LB.unzip $ map (up *** up) $ filter ((/=) 45 . fst) $ LB.zip hum_seq oth_seq
                in Skip (hpos - oldpos) $ Aln ref' smp' $
                   parse1 (fromIntegral (L.length ref') + hpos) lns'

        ls -> error $ "WTF?! (" ++ show oldpos ++ ") near \n" ++ L.unpack (L.unlines $ take 3 ls)

    -- remove degenerate stretches
    norm :: MafStretch -> Stretch
    norm (Skip  0 ms)            = norm ms
    norm (Aln r _ ms) | L.null r = norm ms

    -- join adjacent similar stretches
    norm (Skip n1   (Skip n2   ms)) = norm $ Skip (n1+n2) ms
    norm (Aln r1 s1 (Aln r2 s2 ms)) = norm $ Aln (L.append r1 r2) (L.append s1 s2) ms

    -- emit 'Skip' or 'Aln', but transfer an N to make the length even
    norm (Skip n1 (Aln ref smp ms))
        | L.null ref = norm $ Skip n1 ms
        | even n1    = Ns (fromIntegral $ n1 `div` 2) $ norm $ Aln             ref              smp  ms
        | n1 == 1    =                                  norm $ Aln (L.cons 'X' ref) (L.cons 'N' smp) ms
        | otherwise  = Ns (fromIntegral $ n1 `div` 2) $ norm $ Aln (L.cons 'X' ref) (L.cons 'N' smp) ms

    norm (Aln r1 s1 (Skip n2 ms))
        | n2 == 0            = norm $ Aln r1 s1 ms
        | even (L.length r1) = diff         r1              s1      (norm $ Skip       n2  ms)
        | otherwise          = diff (L.snoc r1 'X') (L.snoc s1 'N') (norm $ Skip (pred n2) ms)

    -- anything else runs into the end
    norm (Aln r1 s1 MafDone) = diff r1 s1 done
    norm (Skip n1   MafDone) = Ns (fromIntegral $ succ n1 `div` 2) done
    norm            MafDone  = done


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
        | rv_gt .&. 0xFF00 == 0xFF00        -- haploid call
            = let z = toUpper $ safeIndex "z" rv_vars ( fromIntegral (rv_gt .&. 0x00FE - 2) )
              in NucCode $ maybe (error $ "What's a " ++ shows z "?") fromIntegral
                         $ B.elemIndex z nuc_chars
        | otherwise                         -- diploid
            = let c1 = toUpper $ safeIndex "c1" rv_vars ( fromIntegral (rv_gt            .&. 0x00FE - 2) )
                  c2 = toUpper $ safeIndex "c2" rv_vars ( fromIntegral (rv_gt `shiftR` 8 .&. 0x00FE - 2) )

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

-- read individual file, then read HEF files for each individual,
-- group them according to population, print allele counts for each
-- variable site
treemix :: Int -> FilePath -> FilePath -> FilePath -> IO ()
treemix noutgr outfile indfile reff = do
    (indivs, popixs, pops) <- foldIndFile <$> B.readFile indfile
    heffalumps <- scanDirectory "."
    let hefs = fmap (\i -> case H.lookup i heffalumps of
                        Nothing -> error $ "no .hef file found for " ++ show i
                        Just hf -> hf) indivs
        npops = length pops

    ref <- readReference reff
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs

    withFile outfile WriteMode $ \hout -> do
        B.hPut hout $ B.unwords pops `B.snoc` '\n'
        forM_ (merge_hefs noutgr ref inps) $ \Variant{..} -> do
            -- XXX  Some variants are probably useless.  Should we
            -- attempt to remove private variants?  Variants where
            -- everyone is different from the reference?

            let refcounts = U.accumulate (+) (U.replicate npops 0) $
                            U.zip popixs $ U.map (fromIntegral . (3 .&.)) v_calls
                altcounts = U.accumulate (+) (U.replicate npops 0) $
                            U.zip popixs $ U.map (fromIntegral . (`shiftR` 2)) v_calls

            let show1 (a,b) k = intDec a <> char7 ',' <> intDec b <> char7 ' ' <> k
            hPutBuilder hout $ U.foldr show1 (char7 '\n') $ U.zip refcounts altcounts


-- | Recursively scans a directory for .hef files, makes a hashmap.
-- This assumes unique names, which have to match the indivividual
-- names.
scanDirectory :: FilePath -> IO (H.HashMap B.ByteString FilePath)
scanDirectory fp = do
    nms <- filter (\n -> n /= "." && n /= "..") <$> getDirectoryContents fp
    foldM (\t n -> do isdir <- doesDirectoryExist (fp </> n)
                      t' <- if isdir then H.union t <$> scanDirectory (fp </> n) else return t
                      return $ if takeExtension n == ".hef" then
                            H.insert (B.pack $ dropExtension n) (fp </> n) t' else t') H.empty nms



-- Folds over an individual file, returns list of individuals, indices
-- into the population list, and the population list.
foldIndFile :: B.ByteString -> ([B.ByteString], U.Vector Int, [B.ByteString])
foldIndFile = go [] [] [] H.empty . filter (not . null) . map B.words
            . filter (not . B.isPrefixOf "#") . B.lines
  where
    go indivs popixs pops _hpops []
        = ( reverse indivs, U.reverse $ U.fromList popixs, reverse pops )

    go indivs popixs pops  hpops ([ ind, _, pop ] : xs)
        = case H.lookup pop hpops of
            Nothing -> go (ind:indivs) (H.size hpops:popixs) (pop:pops)
                          (H.insert pop (H.size hpops) hpops) xs
            Just ix -> go (ind:indivs) (ix:popixs) pops hpops xs

    go _ _ _ _ (ws:_) = error $ "Parse error in individual file at:\n" ++ show (B.unwords ws)


