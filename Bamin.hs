-- | Code to read a BAM file in.  We pileup, then sample a base
-- randomly.  We end with the same format we would get if we ran 'vcflc'
-- on an appropriately generated VCF file.

module Bamin where

import BasePrelude
import Data.ByteString.Builder      ( hPutBuilder )
import Foreign.C
import Foreign.Ptr                  ( Ptr, plusPtr )
import Foreign.Marshal.Array        ( allocaArray )
import Foreign.Marshal.Alloc        ( allocaBytes, alloca )
import Foreign.Storable
import System.Directory             ( renameFile )
import System.Console.GetOpt
import System.IO
import System.Random                ( getStdRandom, next )

import qualified Data.ByteString.Unsafe as B
import qualified Data.ByteString.Lazy as L
import qualified Data.Vector.Unboxed as U

import Stretch
import Util

data Pick = Ostrich | Artificial | Ignore | ArtificialIgnore | Stranded | UDG
  deriving Show

data ConfBam = ConfBam {
    conf_bam_output :: FilePath,
    conf_bam_reference :: FilePath,
    conf_pick :: Pick }
  deriving Show

conf_bam0 :: ConfBam
conf_bam0 = ConfBam (error "no output file specified") (error "no reference file specified") Ostrich

opts_bam :: [ OptDescr ( ConfBam -> IO ConfBam ) ]
opts_bam =
    [ Option "o" ["output"]      (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"]      (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option [ ] ["deaminate"]        (NoArg set_deaminate) "Artificially deaminate"
    , Option [ ] ["ignore-t"]           (NoArg  set_ignore) "Ignore T on forward strand"
    , Option [ ] ["deaminate-ignore-t"] (NoArg set_deamign) "Artificially deaminate, then ignore T"
    , Option [ ] ["stranded"]          (NoArg set_stranded) "Call only G&A on forward strand"
    , Option [ ] ["udg"]                    (NoArg set_udg) "Simulate UDG treatment" ]
  where
    set_output  a c = return $ c { conf_bam_output    = a }
    set_ref     a c = return $ c { conf_bam_reference = a }
    set_deaminate c = return $ c { conf_pick = Artificial }
    set_ignore    c = return $ c { conf_pick = Ignore }
    set_deamign   c = return $ c { conf_pick = ArtificialIgnore }
    set_stranded  c = return $ c { conf_pick = Stranded }
    set_udg       c = return $ c { conf_pick = UDG }

main_bam :: [String] -> IO ()
main_bam args = do
    ( bams, ConfBam{..} ) <- parseOpts True conf_bam0 (mk_opts "bamin" "[bam-file...]" opts_bam) args
    ref <- readReference conf_bam_reference

    withFile (conf_bam_output ++ "~") WriteMode        $ \hdl ->
        hPutBuilder hdl . encode_hap . importPile ref
                . progress conf_bam_output . concat
                =<< mapM (htsPileup conf_pick) bams
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

-- Our varcall: position, base, and counter for the rnd. sampling
data Var0 = Var0 { v_refseq :: !Refseq
                 , v_pos    :: !Int
                 , v_call   :: !NucCode }
    deriving Show

-- Hrm.  Need reference sequence.
importPile :: Reference -> [ Var0 ] -> Stretch
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
htsPileup :: Pick -> FilePath -> IO [ Var0 ]
htsPileup pick fp = do
    plp <- throwErrnoIfNull "c_pileup_init" $ withCString fp $ c_pileup_init
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
                   nuc <- foldVec po vec (fromIntegral n_plp) 0 (U.replicate 5 0)
                   return $ Just $ Var0 rs po nuc

    foldVec po !p !n !i !acc
        | n == i = sample_from acc
        | otherwise = do c <- peekElemOff p i
                         maybe_count pick acc c >>= foldVec po p n (succ i)


maybe_count :: Pick -> U.Vector Int -> CUInt -> IO (U.Vector Int)
maybe_count pick counts c = do
    let -- b_qual =                       c `shiftR`  0 .&. 0x7F
        -- b_posn =        fromIntegral $ c `shiftR` 12 .&. 0x3FF
        -- b_rlen =        fromIntegral $ c `shiftR` 22 .&. 0x3FF
        b_revd =                       c .&. 0x80 /= 0
        b_call =                 case  c `shiftR`  8 .&. 0xF  of
                              b | b == 1    -> 0
                                | b == 2    -> 1
                                | b == 4    -> 2
                                | b == 8    -> 3
                                | otherwise -> 4

    dam <- (==) 0 . (`mod` 50) <$> getStdRandom next

    let nc = case pick of
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
            Ignore | b_call == 0 &&     b_revd -> 4         -- A, reverse
                   | b_call == 3 && not b_revd -> 4         -- T, forward
                   | otherwise                 -> b_call

            -- Introduce artificial deamination, then "deal" with it by
            -- ignoring damaged bases.
            ArtificialIgnore | b_call == 0 -> if     b_revd        then 4 else 0   -- A
                             | b_call == 1 -> if not b_revd && dam then 4 else 1   -- C
                             | b_call == 2 -> if     b_revd && dam then 4 else 2   -- G
                             | b_call == 3 -> if not b_revd        then 4 else 3   -- T
                             | otherwise   ->                           4

            -- Call only G and A on the forward strand, C and T on the
            -- reverse strand.  This doesn't need to be combined with
            -- simulation code: Whenever a C could turn into a T, we ignore
            -- it either way.
            Stranded | b_call == 0 -> if b_revd then 4 else 0             -- A
                     | b_call == 1 -> if b_revd then 1 else 4             -- C
                     | b_call == 2 -> if b_revd then 4 else 2             -- G
                     | b_call == 3 -> if b_revd then 3 else 4             -- T
                     | otherwise   -> 4

            -- Simulated UDG treatment: 2% of Cs on forward strand are
            -- lost, and so are 2% of Gs on reverse strand.
            UDG | b_call == 1 && not b_revd && dam -> 4      -- C, !rev
                | b_call == 2 &&     b_revd && dam -> 4      -- G, rev
                | otherwise                        -> b_call

    return $ counts U.// [(nc, succ $ counts U.! nc)]

-- | Takes a random sample from prepared counts.
sample_from :: U.Vector Int -> IO NucCode
sample_from vec | U.length vec == 5 = do
    let s = U.unsafeIndex vec 0 + U.unsafeIndex vec 1 + U.unsafeIndex vec 2 + U.unsafeIndex vec 3
    if s == 0
        then return (NucCode 0)
        else do
            i <- (`mod` s) <$> getStdRandom next
            return . NucCode $ case () of
                _ | i < U.unsafeIndex vec 0                                             -> 11
                  | i < U.unsafeIndex vec 0 + U.unsafeIndex vec 1                       -> 12
                  | i < U.unsafeIndex vec 0 + U.unsafeIndex vec 1 + U.unsafeIndex vec 2 -> 13
                  | otherwise                                                           -> 14


data PlpAux

-- struct plp_aux_t *pileup_init( const char *fn )
foreign import ccall unsafe "mini-pileup.h pileup_init"
    c_pileup_init :: CString -> IO (Ptr PlpAux)

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

