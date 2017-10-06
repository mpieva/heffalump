module BcfScan ( readBcf, decodeBcf ) where

-- Minimalistic BCF reader.  We only care about the genotype of the
-- first individual.  Should allow for a handful of shortcuts...
--
-- This assumes that GT is the first field encodes for each sample.  If
-- that assumption is violated, you probably get an incomprehensible and
-- misleading error message.  Doing this right just isn't worth the
-- increase in code complexity.

import Bio.Prelude
import Util                     ( decomp )
import VcfScan                  ( RawVariant(..), hashChrom )

import Foreign.C.Types          ( CChar )
import Foreign.Ptr              ( Ptr, plusPtr )
import Foreign.Storable         ( peekByteOff, pokeByteOff )
import Streaming
import System.IO                ( withFile, IOMode(..) )

import qualified Data.ByteString                 as B
import qualified Data.ByteString.Char8           as C
import qualified Data.ByteString.Internal        as B
import qualified Data.ByteString.Streaming       as S
import qualified Data.ByteString.Unsafe          as B
import qualified Data.Vector.Unboxed             as V
import qualified Streaming.Prelude               as Q

readBcf :: FilePath -> (Stream (Of RawVariant) IO () -> IO r) -> IO r
readBcf fp k = withFile fp ReadMode $ k . decodeBcf . decomp . S.fromHandle

-- Skip over header, for now we won't parse it.  This fails if GT is not
-- the first individual field that's encoded, but that should always be
-- the case.
decodeBcf :: MonadIO m => S.ByteString m r -> Stream (Of RawVariant) m r
decodeBcf = lift . S.toStrict . S.splitAt 9 >=> parse_hdr
  where
    parse_hdr (hdr :> rest)
        | "BCF\2" == B.take 4 hdr
            = do let l_text = slow_word32 5 hdr
                 txt :> body <- lift . S.toStrict $ S.splitAt l_text rest
                 lift (S.nextChunk body) >>= \case
                    Left r -> pure r
                    Right (s,ss) -> getvars (mkSymTab txt) ss s
        | otherwise = error "not a BCFv2 file"

    slow_word32 o s =     fromIntegral (B.index s $ o+0)
                      .|. fromIntegral (B.index s $ o+1) `shiftL`  8
                      .|. fromIntegral (B.index s $ o+2) `shiftL` 16
                      .|. fromIntegral (B.index s $ o+3) `shiftL` 24


getvars :: MonadIO m => V.Vector Int -> S.ByteString m r -> B.ByteString -> Stream (Of RawVariant) m r
getvars !tab strs !str = go 0
  where
    go !off
        | B.length str < 32 + off = lift (S.nextChunk strs) >>= \case
                                        Right (s,ss) -> getvars tab ss (B.append (B.drop off str) s)
                                        Left r | B.length str == off -> pure r
                                               | otherwise           -> error "Short record."
        | otherwise = do
            (!l_shared, !l_indiv) <- liftIO . B.unsafeUseAsCString str $ \p -> do
                                            (,) <$> peek32 p off <*> peek32 p (off+4)
            let !l_tot = l_shared + l_indiv + 8

            if B.length str < fromIntegral l_tot + off
                -- This looks kinda unsafe, but I think we should never
                -- hit the end here on properly formatted files.
                then do (hd,tl) <- either (error "Short record.") id <$> lift (S.nextChunk strs)
                        getvars tab tl (B.append (B.drop off str) hd)

                else do !v1 <- liftIO $ B.unsafeUseAsCString str $ \p0 -> do
                                      let !p = plusPtr p0 off
                                      !refid    <- fromIntegral        <$> peek32 p  8
                                      !rv_pos   <- fromIntegral . succ <$> peek32 p 12
                                      !n_allls  <- (`shiftR` 16)       <$> peek32 p 24
                                      !rv_vars  <- get_als n_allls (plusPtr p 32)
                                      !rv_gt    <- get_gts         (plusPtr p (fromIntegral l_shared + 8))
                                      return $! RawVariant{ rv_chrom = tab V.! refid, .. }
                        v1 `Q.cons` go (fromIntegral l_tot + off)


-- skip over variant ID, then get alleles
get_als :: Word32 -> Ptr CChar -> IO B.ByteString
get_als n !p = do !k1 <- peek8 p 0
                  case k1 of
                        0xF7 -> do !k2 <- peek8 p 1
                                   case k2 of
                                        0x01 -> peek8  p 2 >>= kont . (+3) . fromIntegral  -- should be plenty
                                        0x02 -> peek16 p 2 >>= kont . (+4) . fromIntegral  -- but isn't  :-(
                                        0x03 -> peek32 p 2 >>= kont . (+6) . fromIntegral  -- not even close  :,-(
                                        x -> error $ "Huh? " ++ show x
                        tp | tp .&. 0xF == 7 -> kont (1 + fromIntegral (tp `shiftR` 4))
                        _                    -> error "string expected"
  where
    kont !sk = let !p' = plusPtr p sk in get_als' NoFrags 0 n p'


data Frags = Frag !(Ptr CChar) !Int Frags | NoFrags

get_als' :: Frags -> Int -> Word32 -> Ptr CChar -> IO B.ByteString
get_als' !acc !l 0 !_ = B.createUptoN l $ cpfrags acc l
  where
    cpfrags (Frag ps ln fs) !o !p = do forM_ [0..ln-1] $ \i -> do
                                            x <- peekByteOff ps i
                                            pokeByteOff p (o-ln-1+i) (x::Word8)
                                       pokeByteOff p (o-1) (fromIntegral (ord ',') :: Word8)
                                       cpfrags fs (o-ln-1) p
    cpfrags  NoFrags         _  _ = return (l-1)

get_als' !acc !l n !p = do !k1 <- peek8 p 0
                           case k1 of
                                    0xF7 -> do !k2 <- peek8 p 1
                                               case k2 of
                                                    0x01 -> peek8  p 2 >>= kont 3 . fromIntegral -- should be plenty
                                                    0x02 -> peek16 p 2 >>= kont 4 . fromIntegral -- but isn't  :-(
                                                    0x03 -> peek32 p 2 >>= kont 6 . fromIntegral -- not even close  :,-(
                                                    x -> error $ "Huh? " ++ show x
                                    tp | tp .&. 0xF == 7 -> kont 1 (fromIntegral $ tp `shiftR` 4)
                                    _                    -> error "string expected"
  where
    kont !sk !ln = get_als' (Frag (plusPtr p sk) ln acc) (l+ln+1) (n-1) (plusPtr p (sk+ln))


peek8 :: Ptr a -> Int -> IO Word8
peek8 = peekByteOff

peek16 :: Ptr a -> Int -> IO Word16
peek16 = peekByteOff

peek32 :: Ptr a -> Int -> IO Word32
peek32 = peekByteOff

get_gts :: Ptr CChar -> IO Word16
get_gts p = do !k1 <- peek8 p 0
               let !ks = case k1 of 1 -> 2; 2 -> 3; 3 -> 5; _ -> error "WTF?"  -- key size, value ignored
               !tp_byte <- peek8 p ks
                -- we support haploid and diploid, and Word8 and Word16
               case tp_byte of
                    0x11 -> (.|.) 0xFF00 . fromIntegral <$> peek8  p (ks+1)
                    0x12 -> (.|.) 0xFF00 . fromIntegral <$> peek16 p (ks+1)

                    0x21 -> do !x <- fromIntegral <$> peek8 p (ks+1)
                               !y <- fromIntegral <$> peek8 p (ks+2)
                               return $! y `shiftL` 8 .|. x

                    0x22 -> do !x <- peek16 p (ks+1)
                               !y <- peek16 p (ks+3)
                               return $! y `shiftL` 8 .|. x

                    _    -> error $ "only haploid or diploid calls are supported " ++ showHex tp_byte []

mkSymTab :: C.ByteString -> V.Vector Int
mkSymTab = V.fromList . mapMaybe parse . C.lines
  where
    parse l = case C.splitAt (C.length key) l of
        (u,v) | u == key  -> Just $! hashChrom (C.takeWhile (/=',') v)
              | otherwise -> Nothing
    key = "##contig=<ID="

