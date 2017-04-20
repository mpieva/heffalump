-- | Code to read a BAM file in.  We pileup, then sample a base
-- randomly.  We end with the same format we would get if we ran 'vcflc'
-- on an appropriately generated VCF file.

module Bamin where

import BasePrelude
import Data.ByteString.Builder      ( hPutBuilder )
import Data.Fix
import Foreign.C
import Foreign.Ptr                  ( Ptr, plusPtr )
import Foreign.Marshal.Array        ( allocaArray )
import Foreign.Marshal.Alloc        ( allocaBytes, alloca )
import Foreign.Storable
import System.Directory             ( renameFile )
import System.Console.GetOpt
import System.IO
import System.Random                ( StdGen, newStdGen, next )

import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Unboxed.Mutable    as U ( IOVector, new, read, write )

import Lump
import NewRef
import Util

data Pick = Ostrich | Ignore_T | Ignore_A | Ignore_C | Ignore_G
          | Stranded | UDG | Kay | Mateja | Mateja2
  deriving Show

data ConfBam = ConfBam {
    conf_bam_output    :: FilePath,
    conf_bam_reference :: FilePath,
    conf_pick          :: Pick,
    conf_deaminate     :: Bool,
    conf_min_qual      :: Int,
    conf_min_mapq      :: Int,
    conf_ignore_indels :: Bool,
    conf_snip          :: Int,
    conf_snap          :: Int }
  deriving Show

conf_bam0 :: ConfBam
conf_bam0 = ConfBam (error "no output file specified") (error "no reference file specified") Ostrich False 0 0 False 0 0

opts_bam :: [ OptDescr ( ConfBam -> IO ConfBam ) ]
opts_bam =
    [ Option "o" ["output"]             (ReqArg set_output  "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"]          (ReqArg set_ref     "FILE") "Read reference from FILE (.fa)"
    , Option "m" ["min-qual"]           (ReqArg set_minqual "QUAL") "Discard bases below quality QUAL"
    , Option "q" ["min-mapq"]           (ReqArg set_minmapq "QUAL") "Discard reads below mapq QUAL"
    , Option [ ] ["deaminate"]          (NoArg       set_deaminate) "Artificially deaminate"
    , Option [ ] ["ignore-t"]           (set_pick         Ignore_T) "Ignore T on forward strand"
    , Option [ ] ["ignore-a"]           (set_pick         Ignore_A) "Ignore A on forward strand"
    , Option [ ] ["ignore-c"]           (set_pick         Ignore_C) "Ignore C on forward strand"
    , Option [ ] ["ignore-g"]           (set_pick         Ignore_G) "Ignore G on forward strand"
    , Option [ ] ["stranded"]           (set_pick         Stranded) "Call only G&A on forward strand"
    , Option [ ] ["udg"]                (set_pick              UDG) "Simulate UDG treatment"
    , Option [ ] ["kay"]                (set_pick              Kay) "Ts become Ns"
    , Option [ ] ["mateja"]             (set_pick           Mateja) "See man page"
    , Option [ ] ["mateja2"]            (set_pick          Mateja2) "See man page"
    , Option [ ] ["snip"]               (ReqArg set_snip     "NUM") "Ignore the first NUM bases in each aread"
    , Option [ ] ["snap"]               (ReqArg set_snap     "NUM") "Ignore the last NUM bases in each aread"
    , Option "c" ["contiguous"]         (NoArg          set_contig) "Use only reads that align without indels" ]
  where
    set_pick p = NoArg (\c -> return $ c { conf_pick = p })

    set_output  a c = return $ c { conf_bam_output    = a }
    set_ref     a c = return $ c { conf_bam_reference = a }
    set_contig    c = return $ c { conf_ignore_indels = True }
    set_deaminate c = return $ c { conf_deaminate     = True }

    set_minqual a c = readIO a >>= \b -> return $ c { conf_min_qual = b }
    set_minmapq a c = readIO a >>= \b -> return $ c { conf_min_mapq = b }
    set_snip    a c = readIO a >>= \b -> return $ c { conf_snip     = b }
    set_snap    a c = readIO a >>= \b -> return $ c { conf_snap     = b }


main_bam :: [String] -> IO ()
main_bam args = do
    ( bams, cfg@ConfBam{..} ) <- parseOpts True conf_bam0 (mk_opts "bamin" "[bam-file...]" opts_bam) args
    ref <- readTwoBit conf_bam_reference

    rnd_gen <- newStdGen
    withFile (conf_bam_output ++ "~") WriteMode        $ \hdl ->
        hPutBuilder hdl . encode ref . importPile
                . sample_piles rnd_gen conf_deaminate conf_pick ref
                . progressBam conf_bam_output . concat
                =<< mapM (htsPileup cfg) bams
    renameFile (conf_bam_output ++ "~") conf_bam_output


progressBam :: String -> [Var0] -> [Var0]
progressBam fp = go 0 0
  where
    go  _  _ [    ] = []
    go rs po (v:vs)
        | rs /= unRefseq (v_refseq v) || po + 10000000 <= v_loc v
            = unsafePerformIO $ do
                hPutStrLn stderr $ fp ++ "@" ++ show (unRefseq (v_refseq v)) ++ ":" ++ show (v_loc v)
                return $ v : go (unRefseq (v_refseq v)) (v_loc v) vs
        | otherwise =  v : go rs po vs

-- Our varcall: position, base, and either the counters for the random
-- sampling or the variant code
data Var a = Var { v_refseq  :: !Refseq
                 , v_loc     :: !Int
                 , v_call    :: a }
    deriving Show

type Var0 = Var (U.Vector Int)
type Var1 = Var Var2b

sample_piles :: StdGen -> Bool -> Pick -> NewRefSeqs -> [ Var0 ] -> [ Var1 ]
sample_piles g0 deam pick cs0 = nextChrom g0 (nrss_seqs cs0) (Refseq 0)
  where
    nextChrom _ [    ]  _ = const []
    nextChrom g (c:cs) rs = generic g cs c rs 0

    generic _  _ _ !_     _ [          ]     = []
    generic g cs c !rs !pos (var1:mvars)
        | rs /= v_refseq var1
            = nextChrom g cs (succ rs) (var1:mvars)

        | Just (N2b rb,c') <- unconsNRS $ dropNRS (fromIntegral $ v_loc var1 - pos) c
            = let (vv, g')      = if deam then deaminate (v_call var1) g else (v_call var1, g)
                  cc            = post_collect (N2b rb) pick vv
                  (N2b nc, g'') = sample_from cc g'
              in var1 { v_call = V2b (xor rb nc) } :
                 generic g'' cs c' rs (v_loc var1 + 1) mvars

        | otherwise
            = generic g   cs NewRefEnd rs (v_loc var1 + 1) mvars


importPile :: [ Var1 ] -> Fix Lump
importPile = generic (Refseq 0) 0
  where
    generic !_     _ [          ] = Fix (Break (Fix Done))
    generic !rs !pos (var1:mvars)
        -- switch to next chromosome
        | rs /= v_refseq var1 =
            Fix $ Break $ generic (succ rs) 0 (var1:mvars)

        -- gap, creates Ns
        | v_loc var1 >= pos =
            Fix $ Ns (fromIntegral $ v_loc var1 - pos) $
            Fix $ encvar (v_call var1) $
            generic rs (v_loc var1 + 1) mvars

        | otherwise = error $ "Got variant position " ++ show (v_refseq var1, v_loc var1)
                           ++ " when expecting " ++ show pos ++ " or higher."

    encvar (V2b 0) = Eqs1 1     -- ref equals alt, could both be N
    encvar (V2b 1) = Trans1     -- transition
    encvar (V2b 2) = Compl1     -- complement
    encvar (V2b 3) = TCompl1    -- trans-complement
    encvar (V2b _) = Ns 1       -- ref is N, or alt is N



data Buf = Buf { buf_ptr   :: !(Ptr Word8)
               , buf_size  :: !Int
               , buf_flush :: Ptr Word8 -> Int -> IO () }

withBuf :: (Ptr Word8 -> Int -> IO ()) -> (Buf -> IO a) -> IO a
withBuf buf_flush  k = allocaBytes buf_size $ \buf_ptr -> do
        pokeByteOff buf_ptr 0 (0::Word32)
        r <- k Buf{..}
        flush_buffer Buf{..}
        return r
  where
    buf_size = 1000000

flush_buffer :: Buf -> IO ()
flush_buffer Buf{..} = do
    used <- peekByteOff buf_ptr 0
    buf_flush (plusPtr buf_ptr 4) (fromIntegral (used::Word32))
    pokeByteOff buf_ptr 0 (0::Word32)

word8 :: Buf -> Word8 -> IO ()
word8 Buf{..} x = do
    used <- peekByteOff buf_ptr 0
    pokeByteOff buf_ptr (used+4) x
    pokeByteOff buf_ptr 0 (succ used)
    when (buf_size - used < 16) $ flush_buffer Buf{..}


{-# INLINE htsPileup #-}
htsPileup :: ConfBam -> FilePath -> IO [ Var0 ]
htsPileup cfg fp = do
    plp <- throwErrnoIfNull "c_pileup_init" $ withCString fp $
                c_pileup_init (fromIntegral $ conf_min_mapq cfg)
                              (fromIntegral $ fromEnum $ conf_ignore_indels cfg)
    repeatedly (step plp)
  where
    vsize = 2048

    repeatedly k = k >>= maybe (return []) (\a -> (:) a <$> unsafeInterleaveIO (repeatedly k))

    step plp = unsafeInterleaveIO                       $
               alloca                                   $ \ptid ->
               alloca                                   $ \ppos ->
               allocaArray (fromIntegral vsize)         $ \vec ->
               c_pileup_step plp ptid ppos vsize vec  >>= \n_plp ->
               if n_plp < 0 then do
                   c_pileup_done plp
                   return Nothing
               else do
                   rs <- Refseq . fromIntegral <$> peek ptid
                   po <-          fromIntegral <$> peek ppos
                   nuc <- foldVec po vec (fromIntegral n_plp) 0 =<< U.new 10
                   return $ Just $ Var rs po nuc

    foldVec po !p !n !i !acc
        | n /= i = do c <- peekElemOff p i
                      maybe_count cfg acc c
                      foldVec po p n (succ i) acc

        | otherwise = U.freeze acc


-- | Filters for the various cutoffs, then counts the base.  Updates a
-- vector of eight counts:  forward ACGT, reverse ACGT.
maybe_count :: ConfBam -> U.IOVector Int -> CUInt -> IO ()
maybe_count ConfBam{..} counts c = do
    let b_qual =        fromIntegral $ c `shiftR`  0 .&. 0x7F
        b_posn =        fromIntegral $ c `shiftR` 12 .&. 0x3FF
        b_rlen =        fromIntegral $ c `shiftR` 22 .&. 0x3FF
        b_revd =                       c .&. 0x80 /= 0
        -- 0,1,2,3 are A,C,G,T.  4 is something we want to ignore.
        b_call =                 case  c `shiftR`  8 .&. 0xF  of
                              b | b == 1    -> 0
                                | b == 2    -> 1
                                | b == 4    -> 2
                                | b == 8    -> 3
                                | otherwise -> 4
        ix = b_call + if b_revd then 4 else 0

    when ( and [ b_qual >= conf_min_qual
               , b_posn >= conf_snip
               , b_posn <= b_rlen - conf_snap
               , b_call < 4 ] )
         ( U.read counts ix >>= U.write counts ix . succ )


-- | Let's say we count bases A,C,G,T (not N, it's unneeded), separately
-- for the strands.  So we get 8 numbers as input.  We modify the last
-- two, which will become uracil(!) counts.
-- to apply deamination.
deaminate :: U.Vector Int -> StdGen -> (U.Vector Int, StdGen)
deaminate vv gen0 = ( vv', gen2 )
  where
    [a,c,g,t,a',c',g',t',0,0] = U.toList vv
    vv' = U.fromListN 8 [a,c-dc,g,t,a',c',g'-dg,t',dc,dg]

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
    -- reference is C, we ignore the forward Ts and sample from the
    -- rest.
    go Mateja2 | isC       = [ a+a'+u', c+c', g+g', t',     0 ]
               | isG       = [ a,       c+c', g+g', t+t'+u, 0 ]
               | otherwise = [ a+a'+u', c+c', g+g', t+t'+u, 0 ]

    isC = ref == N2b 1
    isG = ref == N2b 3


-- | Takes a random sample from prepared counts.
sample_from :: U.Vector Int -> StdGen -> (Nuc2b, StdGen)
sample_from vec _gen | U.length vec < 5 || U.length vec > 6
    = error "internal error: expected 5 or 6 frequencies"

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

data PlpAux

-- struct plp_aux_t *pileup_init( const char *fn )
foreign import ccall unsafe "mini-pileup.h pileup_init"
    c_pileup_init :: CInt -> CInt -> CString -> IO (Ptr PlpAux)

-- An array of 32 bit words is supplied.  For each base, we store one
-- integer:  bits 0..6:  quality
--           bits 7..7:  reversed?
--           bits 8..11: base
--           bits 12..21: position in the read
--           bits 22..31: length of the read

-- int pileup_step( struct plp_aux_t *data, int *tid, int *pos, int vsize, uint32_t *vec ) ;
foreign import ccall unsafe "mini-pileup.h pileup_step"
    c_pileup_step :: Ptr PlpAux -> Ptr CInt -> Ptr CInt -> CInt -> Ptr CUInt -> IO CInt

-- void pileup_done( struct plp_aux_t *data ) ;
foreign import ccall unsafe "mini-pileup.h pileup_done"
    c_pileup_done :: Ptr PlpAux -> IO ()

newtype Refseq = Refseq { unRefseq :: Int } deriving (Show, Enum, Eq)

