-- | Code to read a BAM file in.  We pileup, then sample a base
-- randomly.  We end with the same format we would get if we ran 'vcflc'
-- on an appropriately generated VCF file.

module Bamin ( main_bam ) where

import Bio.Bam               hiding ( Ns )
import Bio.Bam.Pileup
import Bio.Prelude           hiding ( Ns )
import Data.ByteString.Builder      ( toLazyByteString )
import System.Directory             ( renameFile )
import System.Console.GetOpt
import System.IO
import System.Random                ( StdGen, newStdGen, next )

import qualified Data.ByteString.Lazy           as L
import qualified Data.ByteString.Streaming      as S
import qualified Data.IntMap.Strict             as I
import qualified Data.Vector.Unboxed            as U
import qualified Streaming.Prelude              as Q

import Genome
import Lump
import Util

data Pick = Ostrich | Ignore_T | Ignore_A | Ignore_C | Ignore_G
          | Stranded | UDG | Kay | Mateja | Mateja2
  deriving Show

data ConfBam = ConfBam {
    conf_bam_output    :: FilePath,
    conf_bam_reference :: FilePath,
    conf_pick          :: Pick,
    conf_deaminate     :: Bool,
    conf_min_qual      :: Qual,
    conf_min_mapq      :: Qual,
    conf_snip          :: Int,
    conf_snap          :: Int }
  deriving Show

conf_bam0 :: ConfBam
conf_bam0 = ConfBam (error "no output file specified")
                    (error "no reference file specified")
                    Ostrich False (Q 0) (Q 0) 0 0

opts_bam :: [ OptDescr ( ConfBam -> IO ConfBam ) ]
opts_bam =
    [ Option "o" ["output"]             (ReqArg set_output  "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"]          (ReqArg set_ref     "FILE") "Read reference from FILE (.2bit)"
    , Option "q" ["min-qual"]           (ReqArg set_minqual "QUAL") "Discard bases below quality QUAL"
    , Option "m" ["min-mapq"]           (ReqArg set_minmapq "QUAL") "Discard reads below mapq QUAL"
    , Option [ ] ["deaminate"]          (NoArg       set_deaminate) "Artificially deaminate"
    , Option [ ] ["ignore-t"]           (set_pick         Ignore_T) "Ignore T on forward strand"
    , Option [ ] ["ignore-a"]           (set_pick         Ignore_A) "Ignore A on forward strand"
    , Option [ ] ["ignore-c"]           (set_pick         Ignore_C) "Ignore C on forward strand"
    , Option [ ] ["ignore-g"]           (set_pick         Ignore_G) "Ignore G on forward strand"
    , Option [ ] ["stranded"]           (set_pick         Stranded) "Call only G&A on forward strand"
    , Option [ ] ["udg"]                (set_pick              UDG) "Simulate UDG treatment"
    , Option [ ] ["kay"]                (set_pick              Kay) "Ts become Ns"
    , Option [ ] ["mateja"]             (set_pick           Mateja) "Weird stuff"
    , Option [ ] ["mateja2"]            (set_pick          Mateja2) "Weird stuff"
    , Option [ ] ["snip"]               (ReqArg set_snip     "NUM") "Ignore the first NUM bases in each aread"
    , Option [ ] ["snap"]               (ReqArg set_snap     "NUM") "Ignore the last NUM bases in each aread" ]
  where
    set_pick p = NoArg (\c -> return $ c { conf_pick = p })

    set_output  a c = return $ c { conf_bam_output    = a }
    set_ref     a c = return $ c { conf_bam_reference = a }
    set_deaminate c = return $ c { conf_deaminate     = True }

    set_minqual a c = readIO a >>= \b -> return $ c { conf_min_qual = Q b }
    set_minmapq a c = readIO a >>= \b -> return $ c { conf_min_mapq = Q b }
    set_snip    a c = readIO a >>= \b -> return $ c { conf_snip     =   b }
    set_snap    a c = readIO a >>= \b -> return $ c { conf_snap     =   b }


main_bam :: [String] -> IO ()
main_bam args = do
    ( bams, cfg@ConfBam{..} ) <- parseFileOpts conf_bam0 "bamin" "[bam-file...]" opts_bam args
    ref <- readTwoBit conf_bam_reference
    rnd_gen <- newStdGen

    withFile (conf_bam_output ++ "~") WriteMode        $ \hdl ->
        mergeInputs combineCoordinates bams >=> run    $ \hdr ->
            progressBam "bamin" (meta_refs hdr) 10000 (hPutStr stderr) =$
            filterStream ((>= conf_min_mapq) . b_mapq . unpackBam) =$
            concatMapStream (decompose (DmgToken 0)) =$
            pileup =$
            mapStream (collate_pile cfg) =$
            sample_piles rnd_gen conf_deaminate conf_pick ref =$
            S.hPut hdl (encodePiles ref (meta_refs hdr))
    renameFile (conf_bam_output ++ "~") conf_bam_output

collate_pile :: ConfBam -> Pile -> Var0
collate_pile ConfBam{..} p =
    Var { v_refseq = p_refseq p
        , v_loc    = p_pos p
        , v_call   = U.accum (+) (U.replicate 10 0)
            [ (fromIntegral (unN (db_call b)) +o, 1)
            | (o,b) <- map ((,) 0) (fst $ p_snp_pile p) ++ map ((,) 4) (snd $ p_snp_pile p)
            , db_qual b >= conf_min_qual
            , db_dmg_pos b >= conf_snip || db_dmg_pos b < negate conf_snap ] }


-- Our varcall: position, base, and either the counters for the random
-- sampling or the variant code
data Var a = Var { v_refseq  :: !Refseq
                 , v_loc     :: !Int
                 , v_call    :: a }
    deriving Show

-- This 'Var' stores counts of seen bases in the order "ACGT" on forward
-- strand, "ACGT" on reverse strand, forward N, reverse N.
type Var0 = Var (U.Vector Int)

-- Properly called 'Var' in 2bit encoding.
type Var1 = Var Var2b

sample_piles :: Monad m => StdGen -> Bool -> Pick -> RefSeqs -> Enumeratee [ Var0 ] [ Var1 ] m r
sample_piles g0 deam pick cs0 = eneeCheckIfDone (nextChrom g0 (rss_seqs cs0) (Refseq 0))
  where
    nextChrom _ [    ]  _ = return . liftI
    nextChrom g (c:cs) rs = generic g cs (c ()) rs 0

    generic g cs c !rs !pos k = peekStream >>= \case
        Nothing -> return (liftI k)
        Just var1
            | rs /= v_refseq var1
                -> nextChrom g cs (succ rs) k

            | v_loc var1 < pos -> unexpected "unsorted bam?"

            | Just (N2b rb,c') <- unconsRS $ dropRS (fromIntegral $ v_loc var1 - pos) c
                -> let (vv, g')      = if deam then deaminate (v_call var1) g else (v_call var1, g)
                       cc            = post_collect (N2b rb) pick vv
                       (N2b nc, g'') = sample_from cc g'
                   in headStream >>
                      eneeCheckIfDone (generic g'' cs c' rs (v_loc var1 + 1))
                                      (k (Chunk [var1 { v_call = V2b (xor rb nc) }]))
            | otherwise
                -> headStream >> generic g cs RefEnd rs (v_loc var1 + 1) k

encodePiles :: RefSeqs -> Refs -> S.ByteString (Iteratee [Var1] IO) ()
encodePiles ref tgts = S.mwrap $ do
    map1 <- collect I.empty

    when (I.null map1) . liftIO $ do
        hPutStrLn stderr "Found only unexpected sequences.  Is this the right reference?"
        exitFailure

    return $ do S.fromLazy $ toLazyByteString $ encodeHeader ref
                forM_ [0 .. length (rss_chroms ref) - 1] $
                        \i -> S.fromLazy $ unpackLump $ I.findWithDefault noLump i map1
  where
    collect :: MonadIO m => IntMap PackedLump -> Iteratee [Var1] m (IntMap PackedLump)
    collect m = peekStream >>= maybe (return m) (\v -> scan1 m (v_refseq v) >>= collect)

    scan1 :: MonadIO m => IntMap PackedLump -> Refseq -> Iteratee [Var1] m (IntMap PackedLump)
    scan1 m rs = takeWhileE ((==) rs . v_refseq) >=> lift . run $
                    let rn = sq_name $ getRef tgts rs
                    in case elemIndex rn (rss_chroms ref) of

                        Nothing -> do liftIO $ hPrintf stderr "\nSkipping %s.\n" (unpack rn)
                                      m <$ skipToEof

                        Just  i -> do liftIO $ hPrintf stderr "\nTarget %s becomes index %d.\n" (unpack rn) i
                                      lump Q.:> _ <- encodeLumpToMem $ importPile >> Q.yield Break
                                      liftIO $ hPrintf stderr "\nTarget %s becomes index %d, %d bytes.\n"
                                                    (unpack rn) i (L.length $ unpackLump lump)
                                      return $! I.insert i lump m

importPile :: Monad m => Q.Stream (Q.Of Lump) (Iteratee [Var1] m) ()
importPile = normalizeLump $ generic (Refseq 0) 0
  where
    generic !rs !pos = lift tryHead >>= \case
        Nothing -> Q.yield Break
        Just var1
            -- switch to next chromosome
            | rs /= v_refseq var1 -> do
                forM_ [succ rs .. v_refseq var1] $ \_ -> Q.yield Break
                when (v_loc var1 > 0) $
                    Q.yield $ Ns (fromIntegral $ v_loc var1)
                Q.yield $ enc_var (v_call var1)
                generic (v_refseq var1) (v_loc var1 + 1)

            -- gap, creates Ns
            | v_loc var1 >= pos -> do
                when (v_loc var1 > pos) $
                    Q.yield $ Ns (fromIntegral $ v_loc var1 - pos)
                Q.yield $ enc_var (v_call var1)
                generic rs (v_loc var1 + 1)

            | otherwise -> error $ "Got variant position " ++ show (v_refseq var1, v_loc var1)
                                ++ " when expecting " ++ show pos ++ " or higher."

    enc_var :: Var2b -> Lump
    enc_var (V2b 0) = Eqs1 1     -- ref equals alt, could both be N
    enc_var (V2b 1) = Trans1     -- transition
    enc_var (V2b 2) = Compl1     -- complement
    enc_var (V2b 3) = TCompl1    -- trans-complement
    enc_var (V2b _) = Ns 1       -- ref is N, or alt is N


-- | Let's say we count bases A,C,G,T (not N, it's unneeded), separately
-- for the strands.  So we get 8 numbers as input.  We modify the last
-- two, which will become uracil(!) counts.
-- to apply deamination.
deaminate :: U.Vector Int -> StdGen -> (U.Vector Int, StdGen)
deaminate vv gen0 = ( vv', gen2 )
  where
    [a,c,g,t,a',c',g',t',0,0] = U.toList vv
    vv' = U.fromList [a,c-dc,g,t,a',c',g'-dg,t',dc,dg]

    (dc, gen1) = ifrac 50 c  gen0
    (dg, gen2) = ifrac 50 g' gen1

    -- Integer fraction.  ifrac x f returns either floor(x/f) or
    -- ceil (x/f) such that E(ifrac x f) == x/f.
    ifrac f x g0 = let (p,g1) = next g0
                       (y, r) = x `divMod` f
                   in ( y + fromEnum (p `mod` f < r), g1 )



-- | We receive modified (i.e. deaminated) bases here and only vary
-- the way we count.  No randomness involved here.
post_collect :: Nuc2b -> Pick -> U.Vector Int -> U.Vector Int
post_collect ref pp vv = U.fromListN 5 $ go pp
  where
    -- Most of the time, u acts like t and u' like a'.
    [a,c,g,t,a',c',g',t',u,u'] = U.toList vv

    -- Ostrich selection method:  pick blindly
    go Ostrich = [ a+a'+u', c+c', g+g', t+t'+u, 0 ]

    -- "Deal" with deamination by ignoring possibly broken bases:  T on
    -- forward, A on backward strand.
    go Ignore_T = [ a, c+c', g+g', t', 0 ]

    -- Same as Ignore_T, in case someone wants a control experiment...
    go Ignore_A = [ a'+u', c+c', g+g', t+u, 0 ]
    go Ignore_C = [ a+a'+u', c', g, t+t'+u, 0 ]
    go Ignore_G = [ a+a'+u', c, g', t+t'+u, 0 ]

    -- Call only G and A on the forward strand, C and T on the reverse
    -- strand.  This doesn't need to be combined with simulation code:
    -- Whenever a C could turn into a T, we ignore it either way.
    go Stranded = [ a, c', g, t', 0 ]

    -- Simulated UDG treatment: Uracils vanish without trace.
    go UDG = [ a+a', c+c', g+g', t+t', 0 ]

    -- If we pick a T, we turn it into an N.
    go Kay = [ a, c+c', g+g', t', t+a'+u'+u ]

    -- Mateja's method depends on the reference allele:  if the
    -- reference is C and we see at least one (forward) T, we sample
    -- from reverse-strand reads only.  Same for G/A.
    go Mateja | isC &&  t+u  > 0 = [ a'+u',   c',   g',   t',     0 ]
              | isG && a'+u' > 0 = [ a,       c,    g,    t+u,    0 ]
              | otherwise        = [ a+a'+u', c+c', g+g', t+t'+u, 0 ]

    -- Mateja's second method depends on the reference allele:  if the
    -- reference is C, we ignore the forward Ts and sample from the rest.
    go Mateja2 | isC       = [ a+a'+u', c+c', g+g', t',     0 ]
               | isG       = [ a,       c+c', g+g', t+t'+u, 0 ]
               | otherwise = [ a+a'+u', c+c', g+g', t+t'+u, 0 ]

    isC = ref == N2b 1
    isG = ref == N2b 3


-- | Takes a random sample from prepared counts.
sample_from :: U.Vector Int -> StdGen -> (Nuc2b, StdGen)
sample_from vec _gen | U.length vec < 5
    = error "internal error: expected 5 frequencies"

sample_from vec  gen = do
    let s = U.unsafeIndex vec 0 + U.unsafeIndex vec 1 + U.unsafeIndex vec 2
                + U.unsafeIndex vec 3 + U.unsafeIndex vec 4
    if s == 0
        then (N2b 255, gen)
        else let (ii, gen') = next gen
                 i = ii `mod` s
                 nn = case () of
                        _ | i < U.unsafeIndex vec 0                       -> 2     -- A
                          | i < U.unsafeIndex vec 0 + U.unsafeIndex vec 1 -> 1     -- C
                          | i < U.unsafeIndex vec 0 + U.unsafeIndex vec 1
                              + U.unsafeIndex vec 2                       -> 3     -- G
                          | i < U.unsafeIndex vec 0 + U.unsafeIndex vec 1
                              + U.unsafeIndex vec 2 + U.unsafeIndex vec 3 -> 0     -- T
                          | otherwise                                   -> 255     -- N
             in (N2b nn, gen')

