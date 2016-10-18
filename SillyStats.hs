module SillyStats (
    main_kayvergence,
    main_patterson,
    main_yaddayadda
                  ) where

import BasePrelude
import Numeric.SpecFunctions            ( incompleteBeta )
import System.Console.GetOpt

import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.Vector                     as V
import qualified Data.Vector.Unboxed             as U

import Merge
import Stretch
import Util

data Config = Config {
    conf_blocksize  :: Int,
    conf_noutgroups :: Int,
    conf_nrefpanel  :: Int,
    conf_filter     :: [Variant] -> [Variant],
    conf_msg        :: String,
    conf_nafricans  :: Int,
    conf_reference  :: FilePath }

defaultConfig :: Config
defaultConfig = Config 5000000 1 1 id ""
                       (error "size of reference panel not known")
                       (error "no reference file specified")

transversions :: [Variant] -> [Variant]
transversions = filter tv
  where
    tv v | (v_ref v == 'C' || v_ref v == 'c') && (v_alt v == 'T' || v_alt v == 't') = False
         | (v_ref v == 'T' || v_ref v == 't') && (v_alt v == 'C' || v_alt v == 'c') = False
         | (v_ref v == 'G' || v_ref v == 'g') && (v_alt v == 'A' || v_alt v == 'a') = False
         | (v_ref v == 'A' || v_ref v == 'a') && (v_alt v == 'G' || v_alt v == 'g') = False
         | otherwise                                                                = True

custom :: String -> [Variant] -> [Variant]
custom (a:b:_) = filter valid
  where
    valid v = toUpper (v_ref v) == a && toUpper (v_alt v) == b

showPValue :: Double -> ShowS
showPValue x | x >= 0.002 = showFFloat (Just 3) x
             | otherwise  = showEFloat (Just 1) x

-- Silly counting statistics.  Commonalities:
--
-- - We count stuff, then compute the ratio of some stuff to other stuff.
-- - We perform block jackknife to estimate standard error.
-- - Blocks are weighted according to the value of the denominator above.
-- - The block size is some constant in genomic coordinates.  It doesn't
--   need to be precise, only "large enough".
-- - Given lots of inputs, we compute more lots of statistics.
--
-- We shall scan a list of variants, apply the counting function to each.
-- Then chop into blocks and sum over blocks; keep this table.  Jackknife
-- is a second pass over the table.  We use 'Double's instead of plain
-- 'Int's, because we want to have fractional counts (instead of saampling
-- one representative count randomly) without confusing the Jackknife
-- downstream.  Also, 52 bits for any actual integer counts is still
-- plenty.

type CountFunction = U.Vector Word8 -> U.Vector (Double,Double)

data SillyStats = SillyStats
    { _silly_count       :: Double       -- ^ pseudo-count of successes
    , _silly_total       :: Double       -- ^ pseudo-count of everything
    , _silly_est         :: Double       -- ^ the point estimate
    , _silly_var         :: Double       -- ^ variance (from block jackknife)
    , _silly_p           :: Double }     -- ^ p-value (one-sided beta-test)

gen_stats :: Int -> CountFunction -> [Variant] -> [SillyStats]
gen_stats blk_size cfn vs0 = final_vals
  where
    -- We generate the count (k), the total (n), the result (k/n),
    -- the variance estimate (\sigma^2, from Jackknifing)
    final_vals = [ SillyStats k n (2*k/n -1) (4*v)
                              (if n == 0 then 0/0 else pval (k/n) v)
                 | i <- [0 .. U.length full_counts - 1]
                 , let (k,n) = full_counts U.! i
                 , let v     = blk_jackknife k n $ V.map (U.! i) blockstats ]

    full_counts = V.foldl1' (U.zipWith biplus) blockstats
    blockstats = V.fromList $ foldBlocks vs0

    foldBlocks [] = []
    foldBlocks (v:vs) = case span (near v) vs of
        (ls,rs) -> ((:) $! foldBlock v ls) (foldBlocks rs)

    near v v' = v_chr v == v_chr v' && v_pos v + blk_size > v_pos v'

    foldBlock = foldl' (\acc -> U.zipWith biplus acc . cfn . v_calls) . cfn . v_calls
    biplus (x,y) (z,w) = (x+z,y+w)

    -- p-value.  Use a Beta distribution with a=b=1/8V - 1/2 to match
    -- the variance.  p-value is the CDF evaluated at r (or 1-r,
    -- whichever is <0.5);
    pval r v |  r >= 1   = 0
             |  r <= 0   = 0
             |  r > 0.5  = incompleteBeta ab ab (1-r)
             | otherwise = incompleteBeta ab ab   r
      where ab = recip (8*v) - 0.5

data Pair a b = !a :!: !b deriving (Read, Show, Eq, Ord, Bounded)

infixl 2 :!:

-- Block Jackknife for estimating ratios.  The bias corrected estimator
-- reduces to the plain estimator, so we skip this.  We estimate
-- variance by applying the "Delete-m Jackknife for Unequal m".
--
-- The variance estimator reduces to
-- \[
-- \hat{\sigma}^2 = \frac{1}{G} \sum_{j=1}^G \frac{1}{n-m_j) \left(
--  \frac{k_j^2}{m_j} - 2 \frac{k_j K}{n} + \frac{K^2 m_j}{n^2} \right)
-- \]
--
-- (Pray to FSM that this is actually correct.)
--
-- Arguments are K (sum of the \$k_j\$), N (sum of the \$m_j\$) and the
-- list of pairs of \$(k_j,m_j)\$; result is \$\sigma^2\$.

blk_jackknife :: Double -> Double -> V.Vector (Double,Double) -> Double
blk_jackknife kk nn = divide . V.foldl' step (0 :!: (0::Int))
  where
    divide (s :!: l) = s / fromIntegral l

    step (s :!: l) (kj,mj) | mj == 0   =               s :!: succ l
                           | otherwise = d / (nn-mj) + s :!: succ l
      where
        d = kj*kj / mj - 2 * kj * kk / nn + kk*kk*mj / (nn*nn)

-- --------------- Kayvergence

-- | Produces labels for 'kayvergence' in the same order that statistics
-- are produces.  Input is number of good genomes and labels for
-- individuals, the reference is added internally.
kaylabels :: Int -> [String] -> [(String,String,String)]
kaylabels ngood labels =
    [ ( ref, smp, outg )
    | (iref,ref):goods <- tails $ take (ngood+1) labels'
    , (ioutg,outg) <- goods
    , (ismp,smp) <- labels'
    , ismp /= iref && ismp /= ioutg ]
  where
    labels' = zip [(0::Int) ..] $ "human" : labels

-- | Counter function that computes all sensible Kayvergence ratios.  It
-- receives the number of "good genomes" appearing at the front of the
-- list of calls.  The reference allele is added at the front, it also
-- counts as good.  For the statistic itself, two good genomes are picked
-- (their order doesn't matter) along with another one (good or bad, can be
-- any one, even further up in the list).
--
-- It then counts the total number of differences between the two good
-- genomes, and the number of those where the bad genome matches the first
-- good genome.  Their ratio is the Kayvergence.
kayvergence :: Int -> CountFunction
kayvergence ngood callz = U.fromListN nstats
    [ kaycounts (at ref) (at smp) (at outg)
    | ref <- [0 .. ngood-1]             -- pick from ref and good genomes, but not the last
    , outg <- [ref+1..ngood]            -- pick another good genome
    , smp <- [0 .. ntot-1]              -- pick sample
    , smp /= ref && smp /= outg ]       -- must be different ]

  where
    at i = case callz' U.! i of w -> (fromIntegral $ w .&. 0x3, fromIntegral $ w `shiftR` 2 .&. 0x3)

    callz' = 1 `U.cons` callz       -- the reference allele
    ntot   = U.length callz'        -- number of individuals

    -- number of ways to pick two good genomes and one ordinary one
    nstats = ngood * (ngood+1) * (ntot-2) `div` 2

-- | Counting for 'kayvergence'.  @aa@ is the number of combinations
-- that count as shared variation; @aa+bb@ is the number of combinations
-- that count as variation.  Their ratio will contribute to our final
-- estimate, but we need to weight it by the total number of
-- combinations, which is $\prod_i u_i+v_i$.
-- kaycounts :: Word8 -> Word8 -> Word8 -> (Double, Double)
{-# INLINE kaycounts #-}
kaycounts :: (Double, Double) -> (Double, Double) -> (Double, Double) -> (Double, Double)
kaycounts (u1,v1) (u2,v2) (u3,v3) = if nn > 0 then ( aa/nn, (aa+bb)/nn ) else ( 0,0 )
    where
        aa = u1 * v2 * v3 + v1 * u2 * u3
        bb = u1 * u2 * v3 + v1 * v2 * u3
        nn = (u1 + v1) * (u2 + v2) * (u3 + v3)

opts_kiv :: [ OptDescr (Config -> IO Config) ]
opts_kiv =
    [ Option "r" ["reference"] (ReqArg set_ref  "FILE") "Read reference from FILE (.fa)"
    , Option "n" ["numgood"]   (ReqArg set_ngood "NUM") "The first NUM inputs are \"good\" genomes (1)"
    , Option "J" ["blocksize"] (ReqArg set_jack  "NUM") "Set blocksize for Jackknife to NUM bases (5M)"
    , Option "t" ["transversions"]   (NoArg set_tvonly) "Restrict to transversion sites" ]
  where
    set_ref    a c = return $ c { conf_reference = a }
    set_ngood  a c = readIO a >>= \n -> return $ c { conf_noutgroups = n }
    set_jack   a c = readNumIO a >>= \n -> return $ c { conf_blocksize = n }
    set_tvonly   c = return $ c { conf_filter = transversions }

main_kayvergence :: [String] -> IO ()
main_kayvergence args = do
    ( hefs, Config{..} ) <- parseOpts True defaultConfig (mk_opts "kayvergence" "[hef-file...]" opts_kiv) args
    ref <- readReference conf_reference
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs

    let labels = kaylabels conf_noutgroups hefs
        stats  = gen_stats conf_blocksize (kayvergence conf_noutgroups)
                 $ conf_filter $ merge_hefs False 0 ref inps

    forM_ (zip labels stats) $ \((rn,sn,cn), SillyStats k n r v _) ->
        putStrLn $ "Kiv( " ++ rn ++ "; " ++ sn ++ "; " ++ cn ++ " ) = "
                ++ showFFloat (Just 0) k "/" ++ showFFloat (Just 0) n " = "
                ++ showFFloat (Just 2) (100 * r) "% ± "
                ++ showFFloat (Just 2) (100 * sqrt v) "%"


-- --------------- D-Stats

pattersonlbls :: Int -> Int -> [String] -> [(String,String,String,String)]
pattersonlbls nout nref labels =
    [ ( smp, outg, refA, refB )
    | outg <- outlabels
    , refA:refs' <- tails reflabels
    , refB <- refs'
    , smp <- smplabels ]
  where
    outlabels = take nout labels
    reflabels = take nref $ drop nout labels
    smplabels = drop (nout+nref) labels

-- D-stats.  We compute D(X,O;G,H) where X is the sample, O is an
-- outgroup, G,H are from a reference panel.  So, in analogy to the
-- Kayvergence, we pick one from the outgroup, two from the ref panel,
-- and one sample.
pattersons :: Int -> Int -> CountFunction
pattersons nout nref callz = U.fromListN nstats
    [ abbacounts (ex smp) (ex outg) (ex refA) (ex refB)
    | outg <- U.toList outcalls
    , refA:refs' <- tails $ U.toList refcalls
    , refB <- refs'
    , smp <- U.toList smpcalls ]
  where
    ex w = (fromIntegral $ w .&. 0x3, fromIntegral $ w `shiftR` 2 .&. 0x3)

    outcalls = U.take nout callz
    refcalls = U.take nref $ U.drop nout callz
    smpcalls = U.drop (nout+nref) callz

    -- number of ways to pick one outgroup, two references, and one sample
    nstats = U.length outcalls * U.length smpcalls *
             U.length refcalls * (U.length refcalls -1) `div` 2

abbacounts :: (Double, Double) -> (Double, Double) -> (Double, Double) -> (Double, Double) -> (Double, Double)
abbacounts (u1,v1) (u2,v2) (u3,v3) (u4,v4) = if nn > 0 then ( baba/nn, (abba+baba)/nn ) else ( 0,0 )
      where
        abba = u1 * v2 * v3 * u4 + v1 * u2 * u3 * v4
        baba = v1 * u2 * v3 * u4 + u1 * v2 * u3 * v4
        nn = (u1 + v1) * (u2 + v2) * (u3 + v3) * (u4 + v4)

opts_dstat :: [ OptDescr (Config -> IO Config) ]
opts_dstat =
    [ Option "r" ["reference"]    (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "n" ["numoutgroups"] (ReqArg set_nout "NUM") "The first NUM inputs are outgroups (1)"
    , Option "k" ["numrefpanel"]  (ReqArg set_nref "NUM") "The next NUM inputs are the reference panel"
    , Option "J" ["blocksize"]    (ReqArg set_jack "NUM") "Set blocksize for Jackknife to NUM bases (5M)"
    , Option "t" ["transversions"]     (NoArg set_tvonly) "Restrict to transversion sites"
    , Option [ ] ["custom"]    (ReqArg set_custom "VARS") "Set custom filter" ]
  where
    set_ref  a c = return $ c { conf_reference = a }
    set_nout a c = readIO a >>= \n -> return $ c { conf_noutgroups = n }
    set_nref a c = readIO a >>= \n -> return $ c { conf_nrefpanel = n }
    set_jack a c = readNumIO a >>= \n -> return $ c { conf_blocksize = n }
    set_tvonly c = return $ c { conf_filter = transversions }
    set_custom a c = return $ c { conf_filter = custom a, conf_msg = a ++ " " }

main_patterson :: [String] -> IO ()
main_patterson args = do
    ( hefs, Config{..} ) <- parseOpts True defaultConfig (mk_opts "dstatistics" "[hef-file...]" opts_dstat) args
    ref <- readReference conf_reference
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs

    let labels = pattersonlbls conf_noutgroups conf_nrefpanel hefs
        stats  = gen_stats conf_blocksize (pattersons conf_noutgroups conf_nrefpanel)
                 $ conf_filter $ merge_hefs False conf_noutgroups ref inps

    forM_ (zip labels stats) $ \((sn,cn,r1,r2), SillyStats k n r v p) ->
        putStrLn $ conf_msg
                ++ "D( " ++ r1 ++ ", " ++ r2 ++ "; " ++ sn ++ ", " ++ cn ++ ") = "
                ++ showFFloat (Just 0) k "/" ++ showFFloat (Just 0) n " = "
                ++ showFFloat (Just 2) (100 * r) "% ± "
                ++ showFFloat (Just 2) (100 * sqrt v) "%, p = "
                ++ showPValue p []

-- --------------- Yadda-yadda

yaddalbls :: Int -> Int -> Int -> [String] -> [(String,String,String,String,String)]
yaddalbls nape nout nnea labels =
    [ ( ape, neaA, neaB, smp, afr )
    | ape <- apelabels
    , afr <- afrlabels
    , neaA:neas <- tails nealabels
    , neaB <- neas
    , smp <- smplabels ]
  where
    apelabels = take nape labels
    afrlabels = take nout $ drop nape labels
    nealabels = take nnea $ drop (nout+nape) labels
    smplabels = drop (nape+nout+nnea) labels

-- | Yadda-yadda-counts: We will have to pick one outgroup (C, chimp),
-- two neandertals (N1, N2), one human (H), one human outgroup (Y,
-- yoruba).  Alternatively, we pick two humans and one neandertal, then
-- compute the patterns as (C, H1, H2, N, Y).  The same patterns are
-- counted.  Parameters are the number of true outgroups ('nape',
-- typically the chimp), number of human outgroups ('nout', typically
-- racially pure africans), number of neanderthals ('nnea').
--
-- (Arguably, D-statistics with the chimp as outgroup, two neanderthals,
-- and one human might detect the second admixture while being
-- insensitive to the first.  The same thing could be done with two
-- humans and one neanderthal.  Maybe this is complete bullshit after
-- all.)

yaddayadda :: Int -> Int -> Int -> CountFunction
yaddayadda nape nout nnea callz = U.fromListN nstats
    [ yaddacounts (ex ape) (ex neaA) (ex neaB) (ex smp) (ex afr)
    | ape <- U.toList apecalls
    , afr <- U.toList afrcalls
    , neaA:neas <- tails $ U.toList neacalls
    , neaB <- neas
    , smp <- U.toList smpcalls ]
  where
    ex w = (fromIntegral $ w .&. 0x3, fromIntegral $ w `shiftR` 2 .&. 0x3)

    apecalls = U.take nape callz
    afrcalls = U.take nout $ U.drop nape callz
    neacalls = U.take nnea $ U.drop (nout+nape) callz
    smpcalls = U.drop (nape+nout+nnea) callz

    -- number of ways to pick one outgroup, one african, one sample,
    -- and two neanderthals
    nstats = U.length apecalls * U.length afrcalls * U.length smpcalls
           * U.length neacalls * (U.length neacalls -1) `div` 2

-- | Yadda-yadda:  We compute patterns of the shape (C, N1, N2, H, Y),
-- then count with a positive sign AADDA, AADDD, ADAAD and with a
-- negative sign AADAD, ADADA, ADADD.  (Need to do it twice, because the
-- human reference can be either state.)

yaddacounts :: (Double, Double) -> (Double, Double) -> (Double, Double) -> (Double, Double) -> (Double, Double) -> (Double, Double)
yaddacounts (u1,v1) (u2,v2) (u3,v3) (u4,v4) (u5,v5) = if nn > 0 then ( aa/nn, (aa+bb)/nn ) else ( 0,0 )
      where
        aadda = u1 * u2 * v3 * v4 * u5  +  v1 * v2 * u3 * u4 * v5
        aaddd = u1 * u2 * v3 * v4 * v5  +  v1 * v2 * u3 * u4 * u5
        adaad = u1 * v2 * u3 * u4 * v5  +  v1 * u2 * v3 * v4 * u5

        aadad = u1 * u2 * v3 * u4 * v5  +  v1 * v2 * u3 * v4 * u5
        adada = u1 * v2 * u3 * v4 * u5  +  v1 * u2 * v3 * u4 * v5
        adadd = u1 * v2 * u3 * v4 * v5  +  v1 * u2 * v3 * u4 * u5

        aa = aadda + aaddd + adaad
        bb = aadad + adada + adadd
        nn = (u1 + v1) * (u2 + v2) * (u3 + v3) * (u4 + v4) * (u5 + v5)

opts_yadda :: [ OptDescr (Config -> IO Config) ]
opts_yadda =
    [ Option "r" ["reference"]      (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "n" ["numoutgroups"]   (ReqArg set_nout "NUM") "The first NUM inputs are outgroups (1)"
    , Option "k" ["numafricans"]    (ReqArg set_nafr "NUM") "The next NUM inputs are africans (1)"
    , Option "l" ["numneandertals"] (ReqArg set_nref "NUM") "The next NUM inputs are neanderthals (2)"
    , Option "J" ["blocksize"]      (ReqArg set_jack "NUM") "Set blocksize for Jackknife to NUM bases (5M)"
    , Option "t" ["transversions"]       (NoArg set_tvonly) "Restrict to transversion sites" ]
  where
    set_ref  a c = return $ c { conf_reference = a }
    set_nout a c = readIO a >>= \n -> return $ c { conf_noutgroups = n }
    set_nref a c = readIO a >>= \n -> return $ c { conf_nrefpanel = n }
    set_nafr a c = readIO a >>= \n -> return $ c { conf_nafricans = n }
    set_jack a c = readNumIO a >>= \n -> return $ c { conf_blocksize = n }
    set_tvonly c = return $ c { conf_filter = transversions }

main_yaddayadda :: [String] -> IO ()
main_yaddayadda args = do
    ( hefs, Config{..} ) <- parseOpts True (defaultConfig { conf_nrefpanel = 2 })
                                      (mk_opts "yaddayadda" "[hef-file...]" opts_yadda) args
    ref <- readReference conf_reference
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs

    let labels = yaddalbls conf_noutgroups conf_nafricans conf_nrefpanel hefs
        stats  = gen_stats conf_blocksize (yaddayadda conf_noutgroups conf_nafricans conf_nrefpanel)
                 $ conf_filter $ merge_hefs False (conf_noutgroups+conf_nafricans) ref inps

    forM_ (zip labels stats) $ \((cn,an,n1,n2,sn), SillyStats k n r v p) ->
        putStrLn $ "Y( " ++ cn ++ "; " ++ n1 ++ ", " ++ n2 ++ "; " ++ sn ++ ", " ++ an ++ ") = "
                ++ showFFloat (Just 0) k "/" ++ showFFloat (Just 0) n " = "
                ++ showFFloat (Just 2) (100 * r) "% ± "
                ++ showFFloat (Just 2) (100 * sqrt v) "%, p = "
                ++ showPValue p ""

