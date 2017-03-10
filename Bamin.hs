-- | Code to read a BAM file in.  We pileup, then sample a base
-- randomly.  We end with the same format we would get if we ran 'vcflc'
-- on an appropriately generated VCF file.

module Bamin where

import BasePrelude
import Data.ByteString.Builder      ( hPutBuilder )
import Data.ByteString.Internal     ( c2w )
import Foreign.C
import Foreign.Ptr                  ( Ptr, plusPtr )
import Foreign.Marshal.Array        ( allocaArray )
import Foreign.Marshal.Alloc        ( allocaBytes, alloca )
import Foreign.Storable
import System.Directory             ( renameFile )
import System.Console.GetOpt
import System.IO
import System.Random                ( StdGen, newStdGen, next )

import qualified Data.ByteString.Unsafe as B
import qualified Data.ByteString.Lazy as L
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as U ( IOVector, new, read, write )

import Stretch
import Util

data Pick = Ostrich | Artificial | Ignore_T | Ignore_A | Ignore_C | Ignore_G
          | ArtificialIgnore | Stranded | UDG | Kay | Mateja | Mateja2
          | ArtificialKay | ArtificialMateja | ArtificialMateja2
  deriving Show

data ConfBam = ConfBam {
    conf_bam_output    :: FilePath,
    conf_bam_reference :: FilePath,
    conf_pick          :: Pick,
    conf_min_qual      :: Int,
    conf_min_mapq      :: Int,
    conf_ignore_indels :: Bool,
    conf_snip          :: Int,
    conf_snap          :: Int }
  deriving Show

conf_bam0 :: ConfBam
conf_bam0 = ConfBam (error "no output file specified") (error "no reference file specified") Ostrich 0 0 False 0 0

opts_bam :: [ OptDescr ( ConfBam -> IO ConfBam ) ]
opts_bam =
    [ Option "o" ["output"]             (ReqArg set_output  "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"]          (ReqArg set_ref     "FILE") "Read reference from FILE (.fa)"
    , Option "m" ["min-qual"]           (ReqArg set_minqual "QUAL") "Discard bases below quality QUAL"
    , Option "q" ["min-mapq"]           (ReqArg set_minmapq "QUAL") "Discard reads below mapq QUAL"
    , Option [ ] ["deaminate"]          (set_pick       Artificial) "Artificially deaminate"
    , Option [ ] ["ignore-t"]           (set_pick         Ignore_T) "Ignore T on forward strand"
    , Option [ ] ["ignore-a"]           (set_pick         Ignore_A) "Ignore A on forward strand"
    , Option [ ] ["ignore-c"]           (set_pick         Ignore_C) "Ignore C on forward strand"
    , Option [ ] ["ignore-g"]           (set_pick         Ignore_G) "Ignore G on forward strand"
    , Option [ ] ["deaminate-ignore-t"] (set_pick ArtificialIgnore) "Artificially deaminate, then ignore T"
    , Option [ ] ["stranded"]           (set_pick         Stranded) "Call only G&A on forward strand"
    , Option [ ] ["udg"]                (set_pick              UDG) "Simulate UDG treatment"
    , Option [ ] ["kay"]                (set_pick              Kay) "Ts become Ns"
    , Option [ ] ["deaminate-kay"]      (set_pick    ArtificialKay) "Artificially deaminate, then --kay"
    , Option [ ] ["mateja"]             (set_pick           Mateja) "See man page"
    , Option [ ] ["mateja2"]            (set_pick          Mateja2) "See man page"
    , Option [ ] ["deaminate-mateja"]   (set_pick ArtificialMateja) "Artificially deaminate, then --mateja"
    , Option [ ] ["deaminate-mateja2"] (set_pick ArtificialMateja2) "Artificially deaminate, then --mateja2"
    , Option [ ] ["snip"]               (ReqArg set_snip     "NUM") "Ignore the first NUM bases in each aread"
    , Option [ ] ["snap"]               (ReqArg set_snap     "NUM") "Ignore the last NUM bases in each aread"
    , Option "c" ["contiguous"]         (NoArg          set_contig) "Use only reads that align without indels" ]
  where
    set_pick p = NoArg (\c -> return $ c { conf_pick = p })

    set_output  a c = return $ c { conf_bam_output    = a }
    set_ref     a c = return $ c { conf_bam_reference = a }
    set_contig    c = return $ c { conf_ignore_indels = True }

    set_minqual a c = readIO a >>= \b -> return $ c { conf_min_qual = b }
    set_minmapq a c = readIO a >>= \b -> return $ c { conf_min_mapq = b }
    set_snip    a c = readIO a >>= \b -> return $ c { conf_snip     = b }
    set_snap    a c = readIO a >>= \b -> return $ c { conf_snap     = b }


main_bam :: [String] -> IO ()
main_bam args = do
    ( bams, cfg@ConfBam{..} ) <- parseOpts True conf_bam0 (mk_opts "bamin" "[bam-file...]" opts_bam) args
    (_,ref) <- readReference conf_bam_reference

    rnd_gen <- newStdGen
    withFile (conf_bam_output ++ "~") WriteMode        $ \hdl ->
        hPutBuilder hdl . encode_hap . importPile ref
                . sample_piles rnd_gen conf_pick ref
                . progress conf_bam_output . concat
                =<< mapM (htsPileup cfg) bams
    renameFile (conf_bam_output ++ "~") conf_bam_output

progress :: String -> [Var0] -> [Var0]
progress fp = go 0 0
  where
    go  _  _ [    ] = []
    go rs po (v:vs)
        | rs /= unRefseq (v_refseq v) || po + 10000000 <= v_pos v
            = unsafePerformIO $ do
                hPutStrLn stderr $ fp ++ "@" ++ show (unRefseq (v_refseq v)) ++ ":" ++ show (v_pos v)
                return $ v : go (unRefseq (v_refseq v)) (v_pos v) vs
        | otherwise =  v : go rs po vs

-- Our varcall: position, base, and counters for the rnd. sampling
data Var a = Var { v_refseq  :: !Refseq
                 , v_pos     :: !Int
                 , v_call    :: a }
    deriving Show

type Var0 = Var (U.Vector Int)
type Var1 = Var NucCode

sample_piles :: StdGen -> Pick -> Reference -> [ Var0 ] -> [ Var1 ]
sample_piles g0 pick (Reference cs0) = nextChrom g0 cs0 (Refseq 0)
  where
    nextChrom _ [    ]  _ = const []
    nextChrom g (c:cs) rs = generic g cs c rs 0

    generic _  _ _ !_     _ [          ]     = []
    generic g cs c !rs !pos (var1:mvars)
        | rs /= v_refseq var1 = nextChrom g cs (succ rs) (var1:mvars)
        | L.null c' =                        generic g   cs c' rs (v_pos var1) mvars
        | otherwise = var1 { v_call = nc } : generic g'' cs c' rs (v_pos var1) mvars
      where
        c' = L.drop (fromIntegral $ v_pos var1 - pos) c
        (cc, g') = post_collect (L.head c') pick (v_call var1) g
        (nc, g'') = sample_from cc g'



importPile :: Reference -> [ Var1 ] -> Stretch
importPile (Reference cs0) = nextChrom cs0 (Refseq 0)
  where
    nextChrom [    ]  _ = const $ Break Done
    nextChrom (c:cs) rs = generic cs c rs 0

    generic  _ _ !_     _ [    ]     = Break Done
    generic cs c !rs !pos (var1:mvars)
        | rs /= v_refseq var1 = Break $ nextChrom cs (succ rs) (var1:mvars)

        -- long gap, creates Ns
        | v_pos var1 >= pos + 2 = let l  = (v_pos var1 - pos) `div` 2
                                      c' = L.drop (fromIntegral $ 2*l) c
                                  in Ns (fromIntegral l) $ generic cs c' rs (pos + 2*l) (var1:mvars)

        -- small gap, have to emit codes
        | v_pos var1 == pos + 1 = Chrs (NucCode 0) (v_call var1) $ generic cs (L.drop 2 c) rs (pos+2) mvars

        -- positions must match now
        | v_pos var1 == pos = case mvars of
                -- two variant calls next to each other
                -- if both are reference, we can try and build a stretch
                var2:vars | v_refseq var2 == rs && v_pos var2 == pos+1 ->
                    if isVar c var1 || isVar (L.drop 1 c) var2
                      then Chrs (v_call var1) (v_call var2) $ generic cs (L.drop 2 c) rs (pos+2) vars
                      else matches cs (L.drop 2 c) rs 1 (pos+2) vars

                -- one followed by gap, will become an N
                _ -> do Chrs (v_call var1) (NucCode 0) $ generic cs (L.drop 2 c) rs (pos+2) mvars

        | otherwise = error $ "Got variant position " ++ show (v_refseq var1, v_pos var1)
                           ++ " when expecting " ++ show pos ++ " or higher."

    -- To extend a stretch of matches, we need
    -- two non-vars at the next two positions
    matches cs c !rs !num !pos (var1:var2:vars)
        | v_refseq var1 == rs && v_pos var1 == pos &&
          v_refseq var2 == rs && v_pos var2 == pos+1 &&
          not (isVar c var1) && not (isVar (L.drop 1 c) var2)
            = matches cs (L.drop 2 c) rs (num+1) (pos+2) vars

    -- anything else, we dump the stretch out and pass the buck
    matches cs c rs num pos vars = Eqs1 num $ generic cs c rs pos vars

    isVar c v | L.null  c = False
              | low (L.head c) == 110 = False
              | otherwise = low (tr1 (v_call v)) /= low (L.head c)

    tr1 (NucCode w) = B.unsafeIndex iupac_chars . fromIntegral $ w



data Buf = Buf { buf_ptr :: Ptr Word8
               -- , buf_used :: Int -- will cheat here :(
               , buf_size :: Int
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
                   nuc <- foldVec po vec (fromIntegral n_plp) 0 =<< U.new 8
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

    {- dam <- (==) 0 . (`mod` 50) <$> getStdRandom next

    let nc = case conf_pick of
            _ | b_qual < conf_min_qual      -> 5
              | b_posn < conf_snip          -> 5
              | b_posn > b_rlen - conf_snap -> 5

            -- Ostrich selection method:  pick blindly
            Ostrich -> b_call

            -- Ostrich selection method, also deaminates artificially:  2%
            -- of the time, a C becomes a T instead.  Likewise, 2% of the
            -- time, a G in a reversed-alignment becomes an A.
            Artificial | b_call == 1 && not b_revd && dam -> 3      -- C, !rev
                       | b_call == 2 &&     b_revd && dam -> 0      -- G, rev
                       | otherwise                        -> b_call

            -- "Deal" with deamination by ignoring possibly broken bases:  T
            -- on forward, A on backward strand.
            Ignore_T | b_call == 0 &&     b_revd -> 5         -- A, reverse
                     | b_call == 3 && not b_revd -> 5         -- T, forward
                     | otherwise                 -> b_call

            -- Same as Ignore_T, in case someone wants a control experiment...
            Ignore_A | b_call == 3 &&     b_revd -> 5
                     | b_call == 0 && not b_revd -> 5
                     | otherwise                 -> b_call
            Ignore_C | b_call == 2 &&     b_revd -> 5
                     | b_call == 1 && not b_revd -> 5
                     | otherwise                 -> b_call
            Ignore_G | b_call == 1 &&     b_revd -> 5
                     | b_call == 2 && not b_revd -> 5
                     | otherwise                 -> b_call

            -- Introduce artificial deamination, then "deal" with it by
            -- ignoring damaged bases.
            ArtificialIgnore | b_call == 0 -> if     b_revd        then 5 else 0   -- A
                             | b_call == 1 -> if not b_revd && dam then 5 else 1   -- C
                             | b_call == 2 -> if     b_revd && dam then 5 else 2   -- G
                             | b_call == 3 -> if not b_revd        then 5 else 3   -- T
                             | otherwise   ->                           5

            -- Call only G and A on the forward strand, C and T on the
            -- reverse strand.  This doesn't need to be combined with
            -- simulation code: Whenever a C could turn into a T, we ignore
            -- it either way.
            Stranded | b_call == 0 -> if b_revd then 5 else 0             -- A
                     | b_call == 1 -> if b_revd then 1 else 5             -- C
                     | b_call == 2 -> if b_revd then 5 else 2             -- G
                     | b_call == 3 -> if b_revd then 3 else 5             -- T
                     | otherwise   -> 5

            -- Simulated UDG treatment: 2% of Cs on forward strand are
            -- lost, and so are 2% of Gs on reverse strand.
            UDG | b_call == 1 && not b_revd && dam -> 5      -- C, !rev
                | b_call == 2 &&     b_revd && dam -> 5      -- G, rev
                | otherwise                        -> b_call

            -- If we pick a T, we turn it into an N.
            Kay | b_call == 0 &&     b_revd -> 4         -- A, reverse
                | b_call == 3 && not b_revd -> 4         -- T, forward
                | otherwise                 -> b_call

            ArtificialKay | b_call == 0 -> if     b_revd        then 4 else 0   -- A
                          | b_call == 1 -> if not b_revd && dam then 4 else 1   -- C
                          | b_call == 2 -> if     b_revd && dam then 4 else 2   -- G
                          | b_call == 3 -> if not b_revd        then 4 else 3   -- T
                          | otherwise   ->                           5 -}


    -- return $ counts U.// [(nc, succ $ counts U.! nc)]

-- | Let's say we count bases A,C,G,T (not N, it's unneeded), separately
-- for the strands.  So we get 8 numbers as input.  Depending on the
-- sampling scheme, we now modify the counts into five (this time *with*
-- N) effective counts, then sample.
--
-- XXX  Argh, this doesn't work:  We don't have the reference allele.
-- So we need to count, put the eight counts into the Var0 record, zip
-- with the reference and *then* postprocess and sample.
post_collect :: Word8 -> Pick -> U.Vector Int -> StdGen -> (U.Vector Int, StdGen)
post_collect ref pp vv gen0 = ( U.fromListN 5 $ go pp, gen2 )
  where
    [a,c,g,t,a',c',g',t'] = U.toList vv

    (dc, gen1) = ifrac 50 c gen0
    (dg, gen2) = ifrac 50 g gen1

    -- Integer fraction.  ifrac x f returns either floor(x/f) or ceil
    -- (x/f) such that E(ifrac x f) == x/f.
    ifrac f x g0 = let (p,g1) = next g0
                       (y, r) = x `divMod` f
                   in ( y + fromEnum (p `mod` f < r), g1 )

    -- Ostrich selection method:  pick blindly
    go Ostrich = [ a+a', c+c', g+g', t+t', 0 ]

    -- Ostrich selection method, also deaminates artificially:  2%
    -- of the time, a C becomes a T instead.  Likewise, 2% of the
    -- time, a G in a reversed-alignment becomes an A.
    go Artificial = [ a+a'+dg, c+c'-dc, g+g'-dg, t+t'+dc, 0 ]

    -- "Deal" with deamination by ignoring possibly broken bases:  T on
    -- forward, A on backward strand.
    go Ignore_T = [ a, c+c', g+g', t', 0 ]

    -- Same as Ignore_T, in case someone wants a control experiment...
    go Ignore_A = [ a', c+c', g+g', t, 0 ]
    go Ignore_C = [ a+a', c', g, t+t', 0 ]
    go Ignore_G = [ a+a', c, g', t+t', 0 ]

    -- Introduce artificial deamination, then "deal" with it by ignoring
    -- damaged bases.
    go ArtificialIgnore = [ a, c+c'-dc, g+g'-dg, t', 0 ]

    -- Call only G and A on the forward strand, C and T on the reverse
    -- strand.  This doesn't need to be combined with simulation code:
    -- Whenever a C could turn into a T, we ignore it either way.
    go Stranded = [ a, c', g, t', 0 ]

    -- Simulated UDG treatment: 2% of Cs on forward strand are lost, and
    -- so are 2% of Gs on reverse strand.
    go UDG = [ a+a', c+c'-dc, g+g'-dg, t+t', 0 ]

    -- If we pick a T, we turn it into an N.
    go Kay           = [ a, c+c',    g+g',    t', t+a'       ]
    go ArtificialKay = [ a, c+c'-dc, g+g'-dg, t', t+a'+dc+dg ]

    -- Mateja's method depends on the reference allele:  if the
    -- reference is C and we see at least one (forward) T, we sample
    -- from reverse-strand reads only.  Same for G/A.
    go Mateja | isC && t  > 0 = [ a',     c',   g',   t', 0 ]
              | isG && a' > 0 = [   a,  c,    g,    t,    0 ]
              | otherwise     = [ a+a', c+c', g+g', t+t', 0 ]

    go ArtificialMateja | isC && t +dc > 0 = [   a',   c',   g',   t', 0 ]
                        | isG && a'+dg > 0 = [ a,    c,    g,    t,    0 ]
                        | otherwise        = [ a+a', c+c', g+g', t+t', 0 ]

    -- Mateja's second method depends on the reference allele:  if the
    -- reference is C, we ignore the forward Ts and sample from the
    -- rest.
    go Mateja2 | isC       = [ a+a', c+c', g+g',   t', 0 ]
               | isG       = [ a   , c+c', g+g', t+t', 0 ]
               | otherwise = [ a+a', c+c', g+g', t+t', 0 ]

    go ArtificialMateja2 | isC       = [ a+a'   , c+c'-dc, g+g'   ,   t'   , 0 ]
                         | isG       = [ a      , c+c'-dc, g+g'-dg, t+t'+dc, 0 ]
                         | otherwise = [ a+a'+dg, c+c'-dc, g+g'-dg, t+t'+dc, 0 ]

    isC = ref == c2w 'C' || ref == c2w 'c'
    isG = ref == c2w 'G' || ref == c2w 'g'

-- | Takes a random sample from prepared counts.
sample_from :: U.Vector Int -> StdGen -> (NucCode, StdGen)
sample_from vec gen | U.length vec == 5 || U.length vec == 6 = do
    let s = U.unsafeIndex vec 0 + U.unsafeIndex vec 1 + U.unsafeIndex vec 2
                + U.unsafeIndex vec 3 + U.unsafeIndex vec 4
    if s == 0
        then (NucCode 0, gen)
        else let (ii, gen') = next gen
                 i = ii `mod` s
                 nn = case () of
                        _ | i < U.unsafeIndex vec 0                        -> 11
                          | i < U.unsafeIndex vec 0 + U.unsafeIndex vec 1  -> 12
                          | i < U.unsafeIndex vec 0 + U.unsafeIndex vec 1
                              + U.unsafeIndex vec 2                        -> 13
                          | i < U.unsafeIndex vec 0 + U.unsafeIndex vec 1
                              + U.unsafeIndex vec 2 + U.unsafeIndex vec 3  -> 14
                          | otherwise                                      ->  0
             in (NucCode nn, gen')

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

