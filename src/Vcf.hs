{-# LANGUAGE RankNTypes #-}
module Vcf ( main_xcf, conf_vcf, conf_bcf, toDensity ) where

-- ^ Stuff related to Vcf and Bcf

import Bio.Prelude               hiding ( Ns )
import Data.ByteString.Builder          ( hPutBuilder )
import Data.ByteString.Lazy             ( hPut )
import Streaming
import System.Console.GetOpt
import System.IO

import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Lazy            as L
import qualified Data.IntMap.Strict              as I
import qualified Data.Vector.Unboxed             as U
import qualified Streaming.Prelude               as Q

import BcfScan
import Bed
import Genome
import Lump
import Util
import VcfScan

-- | Position and length to list of 'Lump's.  Used to create either 'Ns'
-- or 'Eqs2' to enable dense and sparse input, also used to create
-- possibly multiple 'Ns' and 'Eqs2' when unbreaking piecewise broken
-- input.
type GapCons = Int -> Int -> [Lump]

-- | One 'GapCons' function for each chromosome.  Usually a const
-- function, but table driven when unbreaking piecewise broken input.
type Density = Int -> GapCons

type StreamRV = Stream (Of RawVariant) IO ()

dense, sparse :: Density
dense  _ _ = pure . Ns
sparse _ _ = pure . Eqs2

data ConfXcf = ConfXcf
    { conf_output  :: FilePath
    , conf_ref     :: FilePath
    , conf_density :: RefSeqs -> StreamRV -> IO (Of Density StreamRV)
    , conf_ploidy  :: Lump -> Lump
    , conf_key     :: String
    , conf_ext     :: String
    , conf_reader  :: [Bytes] -> FilePath -> (StreamRV -> IO ()) -> IO () }

conf_vcf :: ConfXcf
conf_vcf = ConfXcf (error "no output file specified")
                   (error "no reference specified")
                   detect_density id "vcfin" "vcf" readVcf

conf_bcf :: ConfXcf
conf_bcf = conf_vcf { conf_key = "bcfin", conf_ext = "bcf", conf_reader = readBcf }


opts_xcf :: [ OptDescr ( ConfXcf -> IO ConfXcf ) ]
opts_xcf =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.2bit)"
    , Option "D" ["dense"]          (NoArg  set_dense) "Input stores all sites"
    , Option "S" ["sparse"]         (NoArg set_sparse) "Input stores only variants"
    , Option ""  ["patchwork"] (ReqArg set_bed "FILE") "Read HQ regions from FILE (.bed)"
    , Option "L" ["low"]            (NoArg    set_low) "Low coverage, stores haploid calls"
    , Option "H" ["high"]           (NoArg   set_high) "High coverage, stores diploid calls" ]
  where
    set_output a c = return $ c { conf_output  = a }
    set_ref    a c = return $ c { conf_ref     = a }

    set_bed f c = do raw <- L.readFile f
                     let mkdens refs str = do bed <- (parseBed (rss_chroms refs) raw)
                                              return $ toDensity bed :> str
                     return $ c { conf_density = mkdens }
    set_dense    c = return $ c { conf_density = const $ pure . (:>)  dense }
    set_sparse   c = return $ c { conf_density = const $ pure . (:>) sparse }
    set_low      c = return $ c { conf_ploidy  = make_hap }
    set_high     c = return $ c { conf_ploidy  = id       }


main_xcf :: ConfXcf -> [String] -> IO ()
main_xcf conf0 args = do
    ( vcfs, ConfXcf{..} ) <- let opts = mk_opts (conf_key conf0) fspec opts_xcf
                                 fspec = "[" ++ conf_ext conf0 ++ "-file...]"
                             in parseFileOpts conf0 opts args
    ref <- readTwoBit conf_ref
    withMany (conf_reader $ rss_chroms ref) vcfs $ \inp0 -> do
        density :> inp <- conf_density ref inp0
        withFile conf_output WriteMode $ \hdl -> do
            yenome <- Q.fold_ (\y (k,v) -> I.insertWith srsly k v y) I.empty id
                      . mapsM (\s0 -> do
                              -- get one variant to know the chromsome, store it
                              -- together with the packed form of the whole chromsome
                              Right (v1,s1) <- Q.next s0
                              pl :> r <- encodeLumpToMem . Q.map conf_ploidy $
                                         importVcf (density $ rv_chrom v1) (v1 `Q.cons` s1)
                              return $! (rv_chrom v1, pl) :> r )
                      . Q.groupBy ((==) `on` rv_chrom)
                      . progress conf_output . dedupVcf
                      . cleanVcf $ inp

            hPutBuilder hdl $ encodeHeader ref
            forM_ [0 .. length (rss_chroms ref) - 1] $ \i ->
                hPut hdl . unpackLump $ I.findWithDefault noLump i yenome
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

    srsly _ _ = error "Repeated chromosome.  Is the input really sorted by coordinate?"


-- | Some idiot decided to output multiple records for the same position
-- into some VCF files.  If we hit that, we take the first.  (Einmal mit
-- Profis arbeiten!)
-- XXX  The idiocy applies to multiple SNPs at the same location.  Once
-- we support Indels, the logic in here has to change, because Indels
-- usually share the coordinate of a nearby SNP.
dedupVcf :: MonadIO m => Stream (Of RawVariant) m r -> Stream (Of RawVariant) m r
dedupVcf = lift . Q.next >=> \case Left       r  -> pure r
                                   Right (v1,vs) -> go v1 vs
  where
    go v1 = lift . Q.next >=> \case Left r                      -> v1 `Q.cons` pure r
                                    Right (v2,vs) | match v1 v2 ->             go v1 vs
                                                  | otherwise   -> v1 `Q.cons` go v2 vs

    match v1 v2 = rv_chrom v1 == rv_chrom v2 && rv_pos v1 == rv_pos v2


-- | Remove indel variants, since we can't very well use them, and
-- \"variants\" with invalid chromosome numbers, since we clearly don't
-- want them.
-- XXX  If we want to support Indel variants, this has to go.
cleanVcf :: Monad m => Stream (Of RawVariant) m r -> Stream (Of RawVariant) m r
cleanVcf = Q.filter $ \RawVariant{..} ->
    rv_chrom >= 0 && B.length rv_vars == B.length (B.filter (== ',') rv_vars) * 2 + 1

-- Reads a VCF file and returns 'RawVariant's, which is not exactly
-- practical.
readVcf :: [Bytes] -> FilePath -> (Stream (Of RawVariant) IO () -> IO r) -> IO r
readVcf cs fp k = do
    sc <- initVcf fp cs
    k $ Q.untilRight (getVariant sc)

-- | Imports one chromosome worth of VCF style 'RawVariant's into a
-- stream of 'Lump's with a 'Break' tacked onto the end.  This assumes
-- that all variants have the same chromsome, it's best preceeded by
-- 'Streaming.Prelude.groupBy'.  Pseudo-variants of the "no call" type
-- must be filtered beforehand, see 'cleanVcf'.
importVcf :: Monad m => GapCons -> Stream (Of RawVariant) m r -> Stream (Of Lump) m r
importVcf ns = normalizeLump . generic 1 -- Start the coordinate at one, for VCF is one-based.
  where
    generic pos = lift . Q.next >=> \case
        Left r                    -> Q.cons Break $ pure r
        Right (var1,vars)
            -- long gap, creates Ns or Eqs2
            | rv_pos var1 > pos   -> Q.each (ns pos (rv_pos var1 - pos)) >>
                                     generic (rv_pos var1) (var1 `Q.cons` vars)

            -- positions must match now
            | rv_pos var1 < pos   -> error $ "Got variant position " ++ show (rv_pos var1)
                                          ++ " when expecting " ++ show pos ++ " or higher."

            | isVar var1          -> Q.cons (get_var_code var1) $ generic (succ pos) vars
            | otherwise           -> Q.cons (Eqs2     1       ) $ generic (succ pos) vars

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

        -- diploid call
        | otherwise = encTwoVars $ (char_to_2b c1 `xor` char_to_2b n0)
                                 + (char_to_2b c2 `xor` char_to_2b n0) * 4
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


-- | Trying to detect density.  A sequence of entries for, say, 128
-- consecutive positions is a pretty strong tell that the stream is
-- dense.  Otherwise, a gap of, say, 32 positions where the reference
-- does not also have a gap, is a pretty strong tell it is sparse.  Just
-- by looking at the first couple thousand records, we know.  Else, the
-- command line options are still there.
detect_density :: RefSeqs -> StreamRV -> IO (Of Density StreamRV)
detect_density ref =
    Q.toList . Q.splitAt 2048 >=> \(vars :> rest) ->
    if is_dense vars then do
        hPutStrLn stderr "Hmm, this looks like a dense data set."
        return (dense :> (Q.each vars >> rest))
    else if is_sparse vars then do
        hPutStrLn stderr "Hmm, this looks like a sparse data set."
        return (sparse :> (Q.each vars >> rest))
    else do
        hPutStrLn stderr "Hmm, this data set does not make sense to me."
        hPutStrLn stderr "Specify either --dense or --sparse to import it correctly."
        exitFailure
  where
    is_dense = is_dense' 0 0 (0::Int)

    is_dense' !_ !_ 128 _ = True
    is_dense' !_ !_ !_ [] = False
    is_dense' !c !p !n (v:vs) | rv_chrom v /= c = is_dense' (rv_chrom v) (rv_pos v)   1   vs
                              | rv_pos v == p   = is_dense' (rv_chrom v) (rv_pos v)   n   vs
                              | rv_pos v == p+1 = is_dense' (rv_chrom v) (rv_pos v) (n+1) vs
                              | otherwise       = is_dense' (rv_chrom v) (rv_pos v)   1   vs

    is_sparse [    ] = False
    is_sparse (v:vs) = is_sparse' (find_ref $ rv_chrom v) (rv_chrom v) 0 (0::Int) (v:vs)

    is_sparse' _ !_ !_ 32   _         = True
    is_sparse' _ !_ !_ !_ [    ]      = False
    is_sparse' r !c !p !n (v:vs)
        | rv_chrom v /= c             = is_sparse (v:vs)
        | rv_pos v < p                = is_sparse' r c p n vs
        | otherwise = case viewRS r of
            NilRef                   -> is_sparse' r  c (rv_pos v)  n  vs
            l :== r'                 -> is_sparse' r' c (p+l)   0   (v:vs)
            _ :^  r' | rv_pos v == p -> is_sparse' r' c (p+1)   0      vs
                     | otherwise     -> is_sparse' r' c (p+1) (n+1) (v:vs)

    find_ref ix = case drop ix $ rss_seqs ref of s:_ -> s () ; [] -> RefEnd


toDensity :: Bed -> Int -> Int -> Int -> [Lump]
toDensity (Bed vee) rs p0 l0
    = go (fromIntegral p0) (fromIntegral $ p0+l0)
    . U.toList . U.drop (binSearch 0 (U.length vee -1)) $ vee
  where
--    go s e                [] = Ns   (fromIntegral $  e-s) : []                                -- done, default to LQ
--    go s e rgns@((_,s1,e1):rgns')
--        | s1 <= s && s1 <= e && s <= e1 && e <= e1 = Eqs2 (fromIntegral $  e-s) : []            -- done, contained in HQ rgn
--        | s1 <= s && s1 <= e && s <= e1            = Eqs2 (fromIntegral $  e1-s): go e1 e rgns' -- overlap left
--        | s1 <= s && s1 <= e            && e <= e1 = error "We don't support strandedness"      -- TODO
--        | s1 <= s && s1 <= e                       = Ns (fromIntegral $  e-s): []              -- done, contained in LQ rgn
--        | s1 <= s && s1 >= e            && e1 >= e = error "We don't support strandedness"      -- TODO
--        | s1 <= s && s1 >= e                       = error "logic fail"                         -- shouldn't happen
--       
--
--        | s1 >= s && s1 <= e && e1 >= s && e1 >= e = Ns (fromIntegral $ s1-s) : go s1 e  rgns    -- overlap right
--        | s1 >= s && s1 <= e && e1 >= s            = Ns (fromIntegral $ s1-s) : go s1 e1 rgns    -- overlap right
--        | s1 >= s && s1 <= e            && e1 >= e = error "logic fail"                          -- shouldn't happen
--        | s1 >= s && s1 <= e                       = error "We don't support strandedness"       -- TODO
--        | s1 >= s && s1 >= e && e1 >= s && e1 >= e = Ns (fromIntegral $  e-s): []              -- done, contained in LQ rgn
--        | s1 >= s && s1 >= e            && e1 <= e = error "We don't support strandedness"       -- TODO
--        | s1 >= s && s1 >= e                       = error "logic fail"                         -- shouldn't happen
--
    go s e rgns@((r1,s1,e1):rgns')
        | r1 == r0 && s1 < e = case () of
            _ | s1 <= s && e <= e1 -> eqs  (e-s) : []               -- done, contained in HQ rgn
              | s1 <= s && s <  e1 -> eqs (e1-s) : go e1 e rgns'    -- overlap left
              | s1 > s             -> ns  (s1-s) : go s1 e rgns     -- overlap right
              | otherwise          ->              go s  e rgns'    -- shouldn't happen
    go s e _                        = ns   (e-s) : []               -- done, contained in LQ rgn

    r0  = fromIntegral rs
    ns  = Ns   . fromIntegral
    eqs = Eqs2 . fromIntegral

    -- finds the first entry that could overlap [p0,p0+l0)
    binSearch u v
        | u >= v    = u
        | bigger    = binSearch   u   m
        | otherwise = binSearch (m+1) v
      where
        m         = (u+v) `div` 2      -- ensures m >= u && m < v
        (r1,_,e1) = vee U.! m
        bigger    = r1 == r0 && e1 > fromIntegral p0 || r1 > r0

