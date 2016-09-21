-- | Code to read a BAM file in.  We pileup, then sample a base
-- randomly.  We end with the same format we would get if we ran 'vcflc'
-- on an appropriately generated VCF file.

module Bamin where

import Control.Monad                ( when )
import Data.Bits
import Data.ByteString.Builder      ( hPutBuilder )
import Data.Word                    ( Word8, Word32 )
import Foreign.C
import Foreign.Ptr                  ( Ptr, plusPtr )
import Foreign.Marshal.Array        ( allocaArray )
import Foreign.Marshal.Alloc        ( allocaBytes, alloca )
import Foreign.Storable
import System.Console.GetOpt
import System.IO
import System.IO.Unsafe             ( unsafeInterleaveIO )
import System.Random                ( randomRIO )

import qualified Data.ByteString.Unsafe as B
import qualified Data.ByteString.Lazy as L

import Stretch
import Util

data ConfBam = ConfBam {
    conf_bam_output :: FilePath,
    conf_bam_reference :: FilePath,
    conf_pick :: FilePath -> IO [ Var0 ] }
  -- deriving Show

conf_bam0 :: ConfBam
conf_bam0 = ConfBam (error "no output file specified") (error "no reference file specified") ostrich

opts_bam :: [ OptDescr ( ConfBam -> IO ConfBam ) ]
opts_bam =
    [ Option "o" ["output"]      (ReqArg set_output "FILE") "Write output to FILE (.hef)"
    , Option "r" ["reference"]      (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option [ ] ["deaminate"]        (NoArg set_deaminate) "Artificially deaminate"
    , Option [ ] ["ignore-t"]           (NoArg  set_ignore) "Ignore T on forward strand"
    , Option [ ] ["deaminate-ignore-t"] (NoArg set_deamign) "Artificially deaminate, then ignore T"
    , Option [ ] ["stranded"]          (NoArg set_stranded) "Call only G&A on forward strand" ]
  where
    set_output  a c = return $ c { conf_bam_output    = a }
    set_ref     a c = return $ c { conf_bam_reference = a }
    set_deaminate c = return $ c { conf_pick = artificial_ostrich }
    set_ignore    c = return $ c { conf_pick = ignore_t }
    set_deamign   c = return $ c { conf_pick = ignore_artificial_t }
    set_stranded  c = return $ c { conf_pick = stranded }

main_bam :: [String] -> IO ()
main_bam args = do
    ( bams, ConfBam{..} ) <- parseOpts True conf_bam0 (mk_opts "bamin" "[bam-file...]" opts_bam) args
    ref <- readReference conf_bam_reference

    withFile conf_bam_output WriteMode                          $ \hdl ->
        hPutBuilder hdl . encode_hap . importPile ref . concat
                =<< mapM conf_pick bams

-- Our varcall: position, base, and counter for the rnd. sampling
data Var0 = Var0 { v_refseq :: !Refseq
                 , v_pos    :: !Int
                 , v_call   :: !NucCode
                 , v_count  :: !Int }
    deriving Show

zeroVar :: Refseq -> Int -> Var0
zeroVar rs po = Var0 rs po (NucCode 0) 0

data SomeBase = SomeBase { b_call :: !Nucleotides
                         , b_qual :: !Qual
                         , b_posn :: !Int
                         , b_rlen :: !Int
                         , b_revd :: !Bool }
    deriving Show

-- Sample randomly with a transformation function.  The transformation
-- function can try and turn a base into something else, or throw it
-- away by turning it into an N (code 0).
{-# INLINE pickWith #-}
pickWith :: (SomeBase -> IO NucCode) -> (Var0 -> SomeBase -> IO Var0)
pickWith f v0@(Var0 rs po n0 i) sb = do
    nc <- f sb
    case nc of
        NucCode 0 -> return v0
        _         -> do p <- randomRIO (0,i)
                        return $! Var0 rs po (if p == 0 then nc else n0) (succ i)

-- | Ostrich selection method:  pick blindly, but only actual bases.
ostrich :: FilePath -> IO [ Var0 ]
ostrich = htsPileup (pickWith pick) zeroVar
  where
    pick SomeBase{ b_call = b } | b == nucsA = return $ NucCode 11
                                | b == nucsC = return $ NucCode 12
                                | b == nucsG = return $ NucCode 13
                                | b == nucsT = return $ NucCode 14
                                | otherwise  = return $ NucCode  0

-- | Ostrich selection method:  pick blindly, but only actual bases.
-- Also deaminates artificially:  2% of the time, a C becomes a T
-- instead.  Likewise, 2% of the time, a G in a reversed-alignment
-- becomes an A.
artificial_ostrich :: FilePath -> IO [ Var0 ]
artificial_ostrich = htsPileup (pickWith pick) zeroVar
  where
    pick SomeBase{ b_call = b, b_revd = r }
        | b == nucsA          = return $ NucCode 11
        | b == nucsC &&     r = return $ NucCode 12
        | b == nucsG && not r = return $ NucCode 13
        | b == nucsT          = return $ NucCode 14
        | b == nucsC && not r = do p <- randomRIO (0,50::Int)
                                   return . NucCode $ if p == 0 then 14 else 12
        | b == nucsG &&     r = do p <- randomRIO (0,50::Int)
                                   return . NucCode $ if p == 0 then 11 else 13
        | otherwise           = return $ NucCode  0

-- | "Deal" with deamination by ignoring possibly broken bases:  T on
-- forward, C on backward strand.
ignore_t :: FilePath -> IO [ Var0 ]
ignore_t = htsPileup (pickWith pick) zeroVar
  where
    pick SomeBase{ b_call = b, b_revd = r }
        | b == nucsA = return . NucCode $ if   r   then 0 else 11
        | b == nucsC = return $ NucCode 12
        | b == nucsG = return $ NucCode 13
        | b == nucsT = return . NucCode $ if not r then 0 else 14
        | otherwise  = return $ NucCode  0

-- | Introduce artificial deamination, then "deal" with it by ignoring
-- damaged bases.
ignore_artificial_t :: FilePath -> IO [ Var0 ]
ignore_artificial_t = htsPileup (pickWith pick) zeroVar
  where
    pick SomeBase{ b_call = b, b_revd = False }
        | b == nucsA = return $ NucCode 11
        | b == nucsG = return $ NucCode 13
        | b == nucsC = do p <- randomRIO (0,50::Int)
                          return . NucCode $ if p == 0 then 0 else 12
        | otherwise  = return $ NucCode 0

    pick SomeBase{ b_call = b, b_revd = True }
        | b == nucsC = return $ NucCode 12
        | b == nucsT = return $ NucCode 14
        | b == nucsG = do p <- randomRIO (0,50::Int)
                          return . NucCode $ if p == 0 then 0 else 13
        | otherwise  = return $ NucCode  0

-- | Call only G and A on the forward strand, C and T on the reverse
-- strand.  This doesn't need to be combined with simulation code:
-- Whenever a C could turn into a T, we ignore it either way.
stranded :: FilePath -> IO [ Var0 ]
stranded = htsPileup (pickWith pick) zeroVar
  where
    pick SomeBase{ b_call = b, b_revd = False }
        | b == nucsA = return $ NucCode 11
        | b == nucsG = return $ NucCode 13
        | otherwise  = return $ NucCode 0

    pick  SomeBase{ b_call = b, b_revd = True }
        | b == nucsC = return $ NucCode 12
        | b == nucsT = return $ NucCode 14
        | otherwise  = return $ NucCode 0


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

    {- # INLINE getTwo # -}
    -- getTwo k = tryHead >>= \a -> tryHead >>= \b -> k a b

    isVar c v | L.null  c = False
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
htsPileup :: (a -> SomeBase -> IO a) -> (Refseq -> Int -> a) -> FilePath -> IO [ a ]
htsPileup cons nil fp = do
    plp <- throwErrnoIfNull "c_pileup_init" $ withCString fp $ c_pileup_init
    loop (step plp)
  where
    vsize = 256

    loop k = k >>= maybe (return []) (\a -> (:) a <$> unsafeInterleaveIO (loop k))

    step plp = unsafeInterleaveIO                       $
               alloca                                   $ \ptid ->
               alloca                                   $ \ppos ->
               allocaArray (fromIntegral vsize)         $ \vec ->
               c_pileup_step plp ptid ppos vsize vec  >>= \n_plp ->
               if n_plp < 0 then do
                   c_pileup_done plp
                   return Nothing
               else do
                   acc0 <- nil <$> (Refseq . fromIntegral <$> peek ptid)
                               <*> (         fromIntegral <$> peek ppos)
                   Just <$> foldVec vec (fromIntegral n_plp) 0 acc0

    foldVec !p !n !i !acc
        | n == i = return acc
        | otherwise = do c <- peekElemOff p i
                         let somebase = SomeBase
                                { b_qual = Q    . fromIntegral $ c `shiftR`  0 .&. 0x7F
                                , b_revd =                       c .&. 0x80 /= 0
                                , b_call = Nucs . fromIntegral $ c `shiftR`  8 .&. 0xF
                                , b_posn =        fromIntegral $ c `shiftR` 12 .&. 0x3FF
                                , b_rlen =        fromIntegral $ c `shiftR` 22 .&. 0x3FF }
                         cons acc somebase >>= foldVec p n (succ i)


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

newtype Qual = Q { unQ :: Word8 } deriving Show

newtype Nucleotides = Nucs { unNucs :: Word8 } deriving (Show, Eq)

nucsA, nucsC, nucsG, nucsT :: Nucleotides
nucsA = Nucs 1
nucsC = Nucs 2
nucsG = Nucs 4
nucsT = Nucs 8
