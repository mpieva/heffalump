module Lump
    ( Lump(..)
    , mergeLumps
    , mergeLumpsDense
    , PackedLump(..)
    , noLump
    , diff
    , diff2
    , getRefPath
    , encodeLumpToMem
    , decodeMany
    , encode
    , encodeHeader
    , decode
    , debugLump
    , make_hap
    , normalizeLump
    , encTwoVars
    , Frag(..)
    , patch
    ) where

import Bio.Prelude                   hiding ( Ns )
import Data.ByteString.Builder
import Data.ByteString.Internal             ( c2w, w2c )
import Data.ByteString.Unsafe               ( unsafePackMallocCStringLen )
import Data.Hashable                        ( hash )
import Foreign.Marshal.Alloc                ( mallocBytes )
import Streaming
import Streaming.Prelude                    ( mapOf )
import System.Directory                     ( doesFileExist )
import System.IO                            ( withFile, IOMode(..) )

import qualified Data.ByteString            as B
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString.Streaming  as S
import qualified Data.HashSet               as H
import qualified Data.Vector                as V
import qualified Data.Vector.Unboxed        as U
import qualified Streaming.Prelude          as Q
import qualified Stretch                    as S ( Stretch(..) )

import Genome
import Stretch ( Stretch, decode_dip, decode_hap, NucCode(..) )
import Util    ( decomp, unexpected )

-- | Improved diff representation and encoding.  Goals:
--
-- * support one or two known alleles
-- * do not force SNPs into pairs
-- * allow short indels
-- * allow many operations /without access to the reference/
--
-- Outline:
--
-- We store SNP bases as patch to the reference.  Instead of a base, we
-- store the change needed:  none, transition, complement,
-- trans-complement.  Deletions just have a length, insertions have a
-- length and a sequence of bases.  Complex alleles (such as two different
-- alleles) would need to be broken up into separate variants---haven't
-- thought about the details.
--
-- The code we use:
--
--      0b00000000: break
--
--      0b00000001: transition
--      0b00000010: complement
--      0b00000011: trans-complement
--
--      0b00000101: ref+trans
--      0b00000110: trans+trans
--      0b00000111: ref+compl
--      0b00001000: trans+compl
--      0b00001001: compl+compl
--      0b00001010: ref+tcompl
--      0b00001011: trans+tcompl
--      0b00001100: compl+tcompl
--      0b00001101: tcompl+tcompl
--
--      0b0001yxxx: monoallelic indel
--      0b0010yxxx: heterozygous indel
--      0b0011yxxx: diallelic indel
--          if xxx == 0, one length byte follows.
--          if y == 1, it's an insertion and xxx bases follow (2 bits per base,
--                     ATGC, padded to full bytes)
--
--      0b01xxxxxx: short uncalled stretch
--      0b10xxxxxx: short monoallelic matches
--      0b11xxxxxx: short diallelic matches
--          if xxxxxx == 111111, one length byte follows
--          if xxxxxx == 111110, two length bytes follow
--          if xxxxxx == 111101, three length bytes follow
--          if xxxxxx == 111100, four length bytes follow
--
-- We are left with three "reserved" codes (0x04, 0x0E, 0x0F) and three
-- nonsensical ones (0x40, 0x80, 0xC0).
--
-- Note the bases: bit 0 codes for transition, bit 1 for complement.
-- We might sometimes want to output variants without having a reference
-- available.  In that case, we code as follows:  I (Identical), O
-- (transitiOn), P (comPlement), X (trans-complement)

data Lump
    = Ns !Int                          -- ^ uncalled stretch

    | Eqs1 !Int                        -- ^ stretch of monoallelic matches
    | Trans1                           -- ^ monoallelic transversion
    | Compl1                           -- ^ monoallelic complement
    | TCompl1                          -- ^ monoallelic trans-complement

    | Eqs2 !Int                        -- ^ stretch of diallelic matches
    | RefTrans                         -- ^ nine more diallelic SNP codes (*groan*)
    | Trans2
    | RefCompl
    | TransCompl
    | Compl2
    | RefTCompl
    | TransTCompl
    | ComplTCompl
    | TCompl2

    | Del1 !Int                        -- ^ deletion with other allele unknown
    | Del2 !Int                        -- ^ homozygous deletion
    | DelH !Int                        -- ^ deletion, the other allele is the reference

    | Ins1 !Seq1                       -- ^ insertion with other allele unknown
    | Ins2 !Seq1                       -- ^ homozygous insertion
    | InsH !Seq1                       -- ^ insertion, the other allele is the reference

    | Break                            -- ^ break marker (end-of-chromosome)
  deriving Show

type Seq1 = U.Vector Word8

debugLump :: Stream (Of Lump) IO r -> IO r
debugLump = go (0::Int) (0::Int)
  where
    go c i = Q.next >=> \case
      Left r -> return r
      Right (x,l) -> do
            putStrLn $ shows (c,i) "\t" ++ show x
            case x of
                Break  -> go (c+1) 0 l
                Ns   n -> go c (i+n) l
                Eqs1 n -> go c (i+n) l
                Eqs2 n -> go c (i+n) l

                Del1 _ -> go c i l
                Del2 _ -> go c i l
                DelH _ -> go c i l

                Ins1 _ -> go c i l
                Ins2 _ -> go c i l
                InsH _ -> go c i l

                _      -> go c (i+1) l


normalizeLump :: Monad m => Stream (Of Lump) m r -> Stream (Of Lump) m r
normalizeLump = Q.unfoldr (Q.next >=> either (return . Left) go)
  where
    wrapE = either pure (uncurry Q.cons)

    go (Ns 0,l) = Q.next l
    go (Ns n,l) = Q.next l >>= \case
        Right (Ns n',l') -> go (Ns (n+n'),l')
        l'               -> return $ Right (Ns n, wrapE l')

    go (Eqs1 0,l) = Q.next l
    go (Eqs1 n,l) = Q.next l >>= \case
        Right (Eqs1 n',l') -> go (Eqs1 (n+n'),l')
        l'                 -> return $ Right (Eqs1 n, wrapE l')

    go (Eqs2 0,l) = Q.next l
    go (Eqs2 n,l) = Q.next l >>= \case
        Right (Eqs2 n',l') -> go (Eqs2 (n+n'),l')
        l'                 -> return $ Right (Eqs2 n, wrapE l')

    go (Del1 0,l) = Q.next l
    go (Del1 n,l) = Q.next l >>= \case
        Right (Del1 n',l') -> go (Del1 (n+n'),l')
        l'                 -> return $ Right (Del1 n, wrapE l')

    go (Del2 0,l) = Q.next l
    go (Del2 n,l) = Q.next l >>= \case
        Right (Del2 n',l') -> go (Del2 (n+n'),l')
        l'                 -> return $ Right (Del2 n, wrapE l')

    go (DelH 0,l) = Q.next l
    go (DelH n,l) = Q.next l >>= \case
        Right (DelH n',l') -> go (DelH (n+n'),l')
        l'                 -> return $ Right (DelH n, wrapE l')

    go (Ins1 s,l) | U.null s = Q.next l
    go (Ins1 s,l) = Q.next l >>= \case
        Right (Ins1 s',l') -> go (Ins1 (s<>s'),l')
        l'                 -> return $ Right (Ins1 s, wrapE l')

    go (Ins2 s,l) | U.null s = Q.next l
    go (Ins2 s,l) = Q.next l >>= \case
        Right (Ins2 s',l') -> go (Ins2 (s<>s'),l')
        l'                 -> return $ Right (Ins2 s, wrapE l')

    go (InsH s,l) | U.null s = Q.next l
    go (InsH s,l) = Q.next l >>= \case
        Right (InsH s',l') -> go (InsH (s<>s'),l')
        l'                 -> return $ Right (InsH s, wrapE l')

    go lump = return $ Right lump

-- | Lowers all calls to haploid.
-- (Not sure this is really a good idea anymore...)
make_hap :: Lump -> Lump
make_hap  Trans2  = Trans1
make_hap  Compl2  = Compl1
make_hap TCompl2  = TCompl1
make_hap (Eqs2 n) = Eqs1 n
make_hap (Del2 n) = Del1 n
make_hap (Ins2 s) = Ins1 s
make_hap       x  = x

newtype PackedLump = PackLump { unpackLump :: L.ByteString }

-- | The encoding of an empty chromosome:  just the 'Break'
{-# INLINE noLump #-}
noLump :: PackedLump
noLump = PackLump (L.singleton '\0')

{-# INLINE encodeLumpToMem #-}
encodeLumpToMem :: MonadIO m => Stream (Of Lump) m r -> m (Of PackedLump r)
encodeLumpToMem = fmap (mapOf PackLump) . S.toLazy . encodeLump . normalizeLump

-- | Encodes a 'Lump' to a compact binary representation.  Some
-- functionality is duplicated in 'encodeGenome' and 'noLump'; if the
-- encoding is ever changed, it needs to be tracked over there.
{-# INLINE encodeLump #-}
encodeLump :: MonadIO m => Stream (Of Lump) m r -> S.ByteString m r
encodeLump = go
  where
    go s = S.mwrap $ do p <- liftIO $ mallocBytes 0x8000
                        (o,r) <- fill p 0x8000 s 0
                        c <- liftIO $ unsafePackMallocCStringLen (castPtr p,0x8000)
                        return $ do S.chunk (B.take o c)
                                    either pure go r

    -- An insert is never longer than 255, so it needs no more than
    -- ceil(255/4)+2 bytes, hence the magic 66 below.  Once we have 66
    -- bytes, the next code is guaranteed to fit.
    fill p l s o
        | l - o < 66 = return (o, Right s)      -- need more space
        | otherwise  = Q.next s >>= either (return . (,) o . Left) (\case
            (Ns n,s')                         -> stretchOf 0x40 n >>= fill p l s'
            (Eqs1 n,s')                       -> stretchOf 0x80 n >>= fill p l s'
            (Trans1,s')                       -> wrd8 0x01 >>= fill p l s'
            (Compl1,s')                       -> wrd8 0x02 >>= fill p l s'
            (TCompl1,s')                      -> wrd8 0x03 >>= fill p l s'

            (Eqs2 n,s')                       -> stretchOf 0xC0 n >>= fill p l s'
            (RefTrans,s')                     -> wrd8 0x05 >>= fill p l s'
            (Trans2,s')                       -> wrd8 0x06 >>= fill p l s'
            (RefCompl,s')                     -> wrd8 0x07 >>= fill p l s'
            (TransCompl,s')                   -> wrd8 0x08 >>= fill p l s'
            (Compl2,s')                       -> wrd8 0x09 >>= fill p l s'
            (RefTCompl,s')                    -> wrd8 0x0A >>= fill p l s'
            (TransTCompl,s')                  -> wrd8 0x0B >>= fill p l s'
            (ComplTCompl,s')                  -> wrd8 0x0C >>= fill p l s'
            (TCompl2,s')                      -> wrd8 0x0D >>= fill p l s'

            (Del1 n,s')                       -> indelOf 0x10 n >>= fill p l s'
            (DelH n,s')                       -> indelOf 0x20 n >>= fill p l s'
            (Del2 n,s')                       -> indelOf 0x30 n >>= fill p l s'

            (Ins1 sq,s')                      -> indelOf 0x18 (U.length sq) >>= seqOf sq >>= fill p l s'
            (InsH sq,s')                      -> indelOf 0x28 (U.length sq) >>= seqOf sq >>= fill p l s'
            (Ins2 sq,s')                      -> indelOf 0x38 (U.length sq) >>= seqOf sq >>= fill p l s'

            (Break,s')                        -> wrd8 0x00 >>= fill p l s')
      where
        stretchOf k n
            | n < 0x3C      = liftIO $ do pokeElemOff p o (k .|. fromIntegral n)
                                          return $ o+1
            | n < 0x100     = liftIO $ do pokeElemOff p (o+0) (k .|. 0x3F)
                                          pokeElemOff p (o+1) (fromIntegral  n)
                                          return $ o+2
            | n < 0x10000   = liftIO $ do pokeElemOff p (o+0) (k .|. 0x3E)
                                          pokeElemOff p (o+1) (fromIntegral (n             .&. 0xff))
                                          pokeElemOff p (o+2) (fromIntegral (n `shiftR`  8 .&. 0xff))
                                          return $ o+3
            | n < 0x1000000 = liftIO $ do pokeElemOff p (o+0) (k .|. 0x3D)
                                          pokeElemOff p (o+1) (fromIntegral (n             .&. 0xff))
                                          pokeElemOff p (o+2) (fromIntegral (n `shiftR`  8 .&. 0xff))
                                          pokeElemOff p (o+3) (fromIntegral (n `shiftR` 16 .&. 0xff))
                                          return $ o+4
            | otherwise     = liftIO $ do pokeElemOff p (o+0) (k .|. 0x3C)
                                          pokeElemOff p (o+1) (fromIntegral (n             .&. 0xff))
                                          pokeElemOff p (o+2) (fromIntegral (n `shiftR`  8 .&. 0xff))
                                          pokeElemOff p (o+3) (fromIntegral (n `shiftR` 16 .&. 0xff))
                                          pokeElemOff p (o+4) (fromIntegral (n `shiftR` 24 .&. 0xff))
                                          return $ o+5

        indelOf k n
            | n == 0        = error "empty indel"
            | n < 8         = liftIO $ pokeElemOff p o (k .|. fromIntegral n) >> return (o+1)
            | n < 0x100     = liftIO $ pokeElemOff p o k >> pokeElemOff p (o+1) (fromIntegral n) >> return (o+2)
            | otherwise     = error $ "long indel: " ++ show n

        seqOf sq oo
            | U.length sq == 0 = return oo
            | U.length sq == 1 = liftIO $ pokeElemOff p oo (U.unsafeIndex sq 0) >> return (oo+1)
            | U.length sq == 2 = liftIO $ pokeElemOff p oo (U.unsafeIndex sq 0 .|.
                                                            U.unsafeIndex sq 1 `shiftL` 2) >> return (oo+1)
            | U.length sq == 3 = liftIO $ pokeElemOff p oo (U.unsafeIndex sq 0 .|.
                                                            U.unsafeIndex sq 1 `shiftL` 2 .|.
                                                            U.unsafeIndex sq 2 `shiftL` 4) >> return (oo+1)
            | otherwise        = liftIO (pokeElemOff p oo (U.unsafeIndex sq 0 .|.
                                                           U.unsafeIndex sq 1 `shiftL` 2 .|.
                                                           U.unsafeIndex sq 2 `shiftL` 4 .|.
                                                           U.unsafeIndex sq 3 `shiftL` 6))
                                 >> seqOf (U.unsafeDrop 4 sq) (oo+1)

        wrd8 w = liftIO $ pokeElemOff p o w >> return (o+1)

decodeLump :: Monad m => S.ByteString m r -> Stream (Of Lump) m r
decodeLump = Q.unfoldr decode1
  where
    decode1 = S.nextByte >=> \case
        Left r                                 -> return $ Left r
        Right (c,s1)
            | c == 0x00                        -> return $ Right (Break,s1)
            | c == 0x01                        -> return $ Right (Trans1,s1)
            | c == 0x02                        -> return $ Right (Compl1,s1)
            | c == 0x03                        -> return $ Right (TCompl1,s1)
            | c == 0x04                        -> unexpected (show c)
            | c == 0x05                        -> return $ Right (RefTrans,s1)
            | c == 0x06                        -> return $ Right (Trans2,s1)
            | c == 0x07                        -> return $ Right (RefCompl,s1)
            | c == 0x08                        -> return $ Right (TransCompl,s1)
            | c == 0x09                        -> return $ Right (Compl2,s1)
            | c == 0x0A                        -> return $ Right (RefTCompl,s1)
            | c == 0x0B                        -> return $ Right (TransTCompl,s1)
            | c == 0x0C                        -> return $ Right (ComplTCompl,s1)
            | c == 0x0D                        -> return $ Right (TCompl2,s1)
            | c == 0x0E                        -> unexpected (show c)
            | c == 0x0F                        -> unexpected (show c)

            | c .&. 0xF8 == 0x10               -> del_of Del1 c s1
            | c .&. 0xF8 == 0x20               -> del_of DelH c s1
            | c .&. 0xF8 == 0x30               -> del_of Del2 c s1

            | c .&. 0xF8 == 0x18               -> ins_of Ins1 c s1
            | c .&. 0xF8 == 0x28               -> ins_of InsH c s1
            | c .&. 0xF8 == 0x38               -> ins_of Ins2 c s1

            | c .&. 0xC0 == 0x40               -> stretch_of Ns c s1
            | c .&. 0xC0 == 0x80               -> stretch_of Eqs1 c s1
            | c .&. 0xC0 == 0xC0               -> stretch_of Eqs2 c s1

            | otherwise                        -> error $ "Impossibru! " ++ show c

    stretch_of cons c s1
        | c .&. 0x3F == 0x3F  = do Right (w0,s2) <- S.nextByte s1
                                   return $ Right (cons (fromIntegral w0),s2)
        | c .&. 0x3F == 0x3E  = do Right (w0,s2) <- S.nextByte s1
                                   Right (w1,s3) <- S.nextByte s2
                                   return $ Right (cons (fromIntegral w0 .|.
                                                         fromIntegral w1 `shiftL`  8),s3)
        | c .&. 0x3F == 0x3D  = do Right (w0,s2) <- S.nextByte s1
                                   Right (w1,s3) <- S.nextByte s2
                                   Right (w2,s4) <- S.nextByte s3
                                   return $ Right (cons (fromIntegral w0 .|.
                                                         fromIntegral w1 `shiftL`  8 .|.
                                                         fromIntegral w2 `shiftL` 16),s4)
        | c .&. 0x3F == 0x3C  = do Right (w0,s2) <- S.nextByte s1
                                   Right (w1,s3) <- S.nextByte s2
                                   Right (w2,s4) <- S.nextByte s3
                                   Right (w3,s5) <- S.nextByte s4
                                   return $ Right (cons (fromIntegral w0 .|.
                                                         fromIntegral w1 `shiftL`  8 .|.
                                                         fromIntegral w2 `shiftL` 16 .|.
                                                         fromIntegral w3 `shiftL` 24),s5)
        | otherwise           = return $ Right (cons (fromIntegral (c .&. 0x3F)),s1)

    del_of cons c s1
        | c .&. 0x07 == 0  = do Right (hd,tl) <- S.nextByte s1
                                return $ Right (cons (fromIntegral hd),tl)
        | otherwise        = return $ Right (cons (fromIntegral c .&. 0x07),s1)

    ins_of cons c s1
        | c .&. 0x07 == 0  = do Right (hd,tl) <- S.nextByte s1
                                seq_of cons U.empty hd tl
        | otherwise        = seq_of cons U.empty (c .&. 0x07) s1



    seq_of cons acc 0 s1 = return $ Right (cons acc,s1)

    seq_of cons acc 1 s1 = do
        Right (hd,tl) <- S.nextByte s1
        let acc' = acc `U.snoc` (hd `shiftR` 0 .&. 0x3)
        return $ Right (cons acc',tl)

    seq_of cons acc 2 s1 = do
        Right (hd,tl) <- S.nextByte s1
        let acc' = acc `U.snoc` (hd `shiftR` 0 .&. 0x3)
                       `U.snoc` (hd `shiftR` 2 .&. 0x3)
        return $ Right (cons acc',tl)

    seq_of cons acc 3 s1 = do
        Right (hd,tl) <- S.nextByte s1
        let acc' = acc `U.snoc` (hd `shiftR` 0 .&. 0x3)
                       `U.snoc` (hd `shiftR` 2 .&. 0x3)
                       `U.snoc` (hd `shiftR` 4 .&. 0x3)
        return $ Right (cons acc',tl)

    seq_of cons acc n s1 = do
        Right (hd,tl) <- S.nextByte s1
        let acc' = acc `U.snoc` (hd `shiftR` 0 .&. 0x3)
                       `U.snoc` (hd `shiftR` 2 .&. 0x3)
                       `U.snoc` (hd `shiftR` 4 .&. 0x3)
                       `U.snoc` (hd `shiftR` 6 .&. 0x3)
        seq_of cons acc' (n-4) tl


-- | Merging without access to a reference sequence.  This code doesn't
-- believe in Indels and skips over them.
--
-- 'noutgroups':  number of outgroups

mergeLumps :: Monad m => Int -> V.Vector (Stream (Of Lump) m ()) -> Stream (Of [Variant]) m ()
mergeLumps !noutgroups =
    Q.filter (not . null) .
    -- it's only a variant if at least one alt called
    Q.map (filter (U.any has_alt . U.drop noutgroups . v_calls)) .
    mergeLumpsWith (V.minimum . V.map skiplen . V.drop noutgroups)
  where
    has_alt c = c .&. 0xC /= 0

    skiplen (Just (Ns   n)) = n
    skiplen (Just (Eqs1 n)) = n
    skiplen (Just (Eqs2 n)) = n
    skiplen (Just  Break  ) = maxBound
    skiplen  Nothing        = maxBound
    skiplen  _              = 0

-- Like mergeLumps, but produces at least one SNP record for each site.
-- Since all variants are produced (and then some), outgroups need no
-- special treatment.
mergeLumpsDense :: Monad m => V.Vector (Stream (Of Lump) m ()) -> Stream (Of [Variant]) m ()
mergeLumpsDense =
    -- If there is no actual variant in here, turn the first into
    -- pseudo-var; else keep them all
    Q.map (\vs ->
        let vs' = filter (U.any has_alt . v_calls) vs
        in if null vs' then [ (head vs) { v_alt = V2b 255 } ] else vs' ) .
    Q.filter (not . null) .
    mergeLumpsWith (V.minimum . V.map skiplen)
  where
    has_alt c = c .&. 0xC /= 0

    -- Don't skip over equal stretches.  This will cause the code in
    -- mergeLumpsWith to attempt to create four variants at each
    -- position, unless all samples have no calls.  (This neatly skips
    -- over N stretches in the reference.)
    skiplen (Just (Ns   n)) = n
    skiplen (Just (Eqs1 _)) = 0
    skiplen (Just (Eqs2 _)) = 0
    skiplen (Just  Break  ) = maxBound
    skiplen  Nothing        = maxBound
    skiplen  _              = 0


-- | Merges multiple 'Lump's and generates a stream of 'Variant's.  This
-- is efficient by skipping over long invariant stretches.  The
-- 'skipLen' argument is expected to compute how far ahead the process
-- could jump, but it is only allowed to look at the first code in the
-- 'Lump'y stream.
mergeLumpsWith :: Monad m => (V.Vector (Maybe Lump) -> Int) -> V.Vector (Stream (Of Lump) m r) -> Stream (Of [Variant]) m ()
mergeLumpsWith skipLen = go 0 0
    -- Merging stretches.  We define a 'Variant' as anything that is
    -- different from the reference.  Therefore, 'Eqs' ('Eqs1') and 'Ns'
    -- never create a 'Variant' and we can skip forwards.  Else, if at
    -- least one sample has some sort of call, we need to collect the
    -- alleles, produce up to four variants, then output and skip
    -- forward by one base.
  where
    go :: Monad m => Int -> Int -> V.Vector (Stream (Of Lump) m r) -> Stream (Of [Variant]) m ()
    go !ix !pos = lift . V.mapM Q.next >=> \case
             -- all samples 'Done', no more 'Variant's to produce
        smps | V.all isLeft  smps -> pure ()

             -- all samples 'Break' or are done, next chromosome
             | V.all isBreak smps -> go (succ ix) 0 (V.map skipBreak smps)

             -- We ignore outgroups in computing the longest stretch.
             -- That way, outgroups don't get to define variants, but
             -- participate in variants found in another way.
             | otherwise ->
                -- pure operation, doesn't consume the streams
                let msmps = V.map (either (const Nothing) (Just . fst)) smps
                in case skipLen msmps of
                    -- no stretch, have to look for vars
                    0 -> mkVar ix pos msmps `Q.cons`
                         go ix (succ pos) (unS $ V.mapM (S . skipStretchE 1) smps)

                    -- a stretch of l non-variants can be skipped over
                    l -> go ix (pos + fromIntegral l) (unS $ V.mapM (S . skipStretchE l) smps)

    isBreak (Left         _ ) = True
    isBreak (Right (Break,_)) = True
    isBreak               _   = False

    skipBreak (Left         x ) =     pure x
    skipBreak (Right (Break,s)) =          s
    skipBreak (Right     (a,s)) = Q.cons a s

    -- Skip over exactly l sites.
    skipStretchS :: Monad m => Int -> Stream (Of Lump) m r -> Stream (Of Lump) m r
    skipStretchS !l = lift . Q.next >=> skipStretchE l

    skipStretchE :: Monad m => Int -> Either r (Lump, Stream (Of Lump) m r) -> Stream (Of Lump) m r
    skipStretchE !_ (Left  r) = pure r
    skipStretchE !l (Right s) = skipStretch l s

    skipStretch :: Monad m => Int -> (Lump, Stream (Of Lump) m r) -> Stream (Of Lump) m r
    skipStretch !l _ | l < 0  = unexpected ""

    skipStretch !l (Del1 _,s) = skipStretchS l s
    skipStretch !l (Del2 _,s) = skipStretchS l s
    skipStretch !l (DelH _,s) = skipStretchS l s

    skipStretch !l (Ins1 _,s) = skipStretchS l s
    skipStretch !l (Ins2 _,s) = skipStretchS l s
    skipStretch !l (InsH _,s) = skipStretchS l s

    skipStretch !_ (Break, s) = Break `Q.cons` s

    skipStretch !l (Ns   !n,s) | l == n    = s
                               | l <  n    = Ns (n-l) `Q.cons` s
                               | otherwise = skipStretchS (l-n) s
    skipStretch !l (Eqs1 !n,s) | l == n    = s
                               | l <  n    = Eqs1 (n-l) `Q.cons` s
                               | otherwise = skipStretchS (l-n) s
    skipStretch !l (Eqs2 !n,s) | l == n    = s
                               | l <  n    = Eqs2 (n-l) `Q.cons` s
                               | otherwise = skipStretchS (l-n) s

    skipStretch  0      (a,s) = Q.cons a s
    skipStretch !l      (_,s) = skipStretchS (l-1) s

    -- Find all the variants, anchored on the reference allele, and
    -- split them.  Misfitting alleles are not counted.
    mkVar :: Int -> Int -> V.Vector (Maybe Lump) -> [Variant]
    mkVar ix pos ss = [ Variant ix pos ref_indet (V2b alt) calls
                      | (alt, ct) <- zip [1..3] [ct_trans, ct_compl, ct_tcompl]
                      , let calls = V.convert $ V.map (maybe 0 ct) ss ]

    ref_indet = N2b 255

    -- Variant codes:  #ref + 4 * #alt
    ct_trans :: Lump -> Word8
    ct_trans (Eqs1      _) = 1
    ct_trans  Trans1       = 4

    ct_trans (Eqs2      _) = 2
    ct_trans  RefTrans     = 5
    ct_trans  Trans2       = 8
    ct_trans  _            = 0

    ct_compl :: Lump -> Word8
    ct_compl (Eqs1      _) = 1
    ct_compl  Compl1       = 4

    ct_compl (Eqs2      _) = 2
    ct_compl  RefCompl     = 5
    ct_compl  Compl2       = 8
    ct_compl  _            = 0

    ct_tcompl :: Lump -> Word8
    ct_tcompl (Eqs1      _) = 1
    ct_tcompl  TCompl1      = 4

    ct_tcompl (Eqs2      _) = 2
    ct_tcompl  RefTCompl    = 5
    ct_tcompl  TCompl2      = 8
    ct_tcompl  _            = 0


-- | This gunk is need to make a map over a 'Vector' strict.  Looks
-- ugly, but is perfectly servicable.
newtype S a = S { unS :: a }

instance Functor S where
    {-# INLINE fmap #-}
    fmap f (S a) = S (f a)

instance Applicative S where
    {-# INLINE pure #-}
    pure = S
    {-# INLINE (<*>) #-}
    S f <*> S a = S (f a)

instance Monad S where
    {-# INLINE return #-}
    return = S
    {-# INLINE (>>=) #-}
    S !a >>= k = k a


-- | Converts old style 'Stretch' to new style 'Lump'.  Unfortunately,
-- this needs a reference, but at least we can use a new style
-- 'RefSeq', too.
stretchToLump :: Monad m => RefSeqs -> Stream (Of Stretch) m r -> Stream (Of Lump) m r
stretchToLump nrs = normalizeLump . go1 (rss_seqs nrs)
  where
    go1 [     ] l = lift $ Q.effects l
    go1 (r0:rs) l = go (r0 ()) l
      where
        go r = lift . Q.next >=> \case
            Right (S.Chrs a b,k) -> call a (call b (`go` k)) r
            Right (S.Ns     c,k) -> Ns   (fromIntegral $ c+c) `Q.cons` go (dropRS (fromIntegral $ c+c) r) k
            Right (S.Eqs    c,k) -> Eqs2 (fromIntegral $ c+c) `Q.cons` go (dropRS (fromIntegral $ c+c) r) k
            Right (S.Eqs1   c,k) -> Eqs1 (fromIntegral $ c+c) `Q.cons` go (dropRS (fromIntegral $ c+c) r) k
            Right (S.Break   ,k) -> Break                     `Q.cons` go1 rs k
            Left z               -> pure z

    call :: Monad m => NucCode -> (RefSeq -> Stream (Of Lump) m r) -> RefSeq -> Stream (Of Lump) m r
    call (NucCode b) k r = case unconsRS r of
        Nothing     -> lift $ Q.effects (k r)
        Just (a,r') -> encTwoNucs a (nc2vs U.! fromIntegral b) `Q.cons` k r'

    nc2vs :: U.Vector Word8
    !nc2vs = U.fromList [ 255,  5, 10, 15, 0,        -- NACGT
                            6,  7,  1, 11, 2, 3,     -- MRWSYK
                           17, 18, 19, 16, 255 ]     -- acgtN


-- | Main decoder.  Switches behavior based on header.  "HEF\0",
-- "HEF\1" are the old format and go through conversion;  For
-- the old format, we need a reference, and we have to trust it's a
-- compatible one.  "HEF\2" is the new format.  If a reference is given,
-- we check if we can use it.
--
-- We try to support legacy files.  These could start with anything, but
-- in practice start with either (Ns 5000) or (Ns 5001).  Anything else
-- is rejected as unknown junk.
decode :: Monad m => Either String RefSeqs -> S.ByteString m r -> Stream (Of Lump) m r
decode (Left err) str = do hd :> tl <- lift . S.toStrict $ S.splitAt 4 str
                           if hd == "HEF\3"
                               then decodeLump . S.drop 5 $ S.dropWhile (/= 0) tl
                               else error err

decode (Right  r) str = do hd :> tl <- lift . S.toStrict $ S.splitAt 4 str
                           case B.take 3 hd :> B.drop 3 hd of
                               "HEF" :> "\0"     -> stretchToLump r $ decode_dip tl
                               "HEF" :> "\1"     -> stretchToLump r $ decode_hap tl
                               "HEF" :> "\3"     -> go tl
                               "\176\113\2" :> _ -> stretchToLump r $ decode_dip $ S.chunk hd >> tl
                               "\177\113\2" :> _ -> stretchToLump r $ decode_dip $ S.chunk hd >> tl
                               _                 -> error "File format not recognized?"
  where
    go s = do hs :> s' <- lift . S.toStrict . S.splitAt 4 . S.drop 1 . S.dropWhile (/= 0) $ s
              let hv = B.foldr (\b a -> fromIntegral b .|. shiftL a 8) 0 hs :: Word32
              if hv == fromIntegral (hash (rss_chroms r, rss_lengths r))
                  then decodeLump s'
                  else error "Incompatible reference genome."


getRefPath :: Monad m => S.ByteString m r -> m (Maybe FilePath, S.ByteString m r)
getRefPath str = do
    key :> rest <- S.toStrict $ S.splitAt 4 str
    if "HEF\3" == key
      then do fp :> rest' <- S.toStrict $ S.break (== 0) rest
              return ( Just $ unpack fp, S.chunk key >> S.chunk fp >> rest' )
      else return ( Nothing, S.chunk key >> rest )

-- I don't like to infect everything with 'ResourceT', and apparently
-- the only way around it is to use 'withFile', so we might as well do
-- it right in here.
decodeMany :: Maybe FilePath -> [FilePath]
           -> ( Either String RefSeqs -> V.Vector (Stream (Of Lump) IO ()) -> IO r ) -> IO r
decodeMany mrs fs kk =
    withFiles fs ReadMode $ \hdls -> do
        (rps,raws) <- unzip <$> mapM (getRefPath . decomp . S.fromHandle) hdls
        rs <- case mrs of
                Just fp -> Right <$> liftIO (readTwoBit fp)
                Nothing -> do
                    fps <- filterM (liftIO . doesFileExist) $ catMaybes rps
                    case fps of
                        fp : _ -> Right <$> liftIO (readTwoBit fp)
                        [    ] -> return . Left $ "No reference found.  Looked for it at "
                                    ++ intercalate ", " (nubHash $ catMaybes rps) ++ "."
        kk rs (V.fromList $ map (decode rs) raws)
  where
    withFiles [      ] _iom k = k []
    withFiles (fp:fps)  iom k =
        withFile fp iom $ \hdl ->
            withFiles fps iom $
                k . (:) hdl

    nubHash = H.toList . H.fromList


-- | Encode a 'Lump' and enough information about the 'RefSeqs' to be
-- (1) able to find it again and (2) to make sure we got the right one
-- when operating on a 'Lump'.
encode :: MonadIO m => RefSeqs -> Stream (Of Lump) m r -> S.ByteString m r
encode r s = do S.toStreamingByteString (encodeHeader r)
                encodeLump s

encodeHeader :: RefSeqs -> Builder
encodeHeader r = byteString "HEF\3" <>
                 byteString (rss_path r) <> word8 0 <>
                 int32LE (fromIntegral (hash (rss_chroms r, rss_lengths r)))

-- | Diff between a reference and a sample in string format.  The sample
-- is treated as diploid.
--
-- XXX  We should deal with gaps, now that we can actually encode them
{-# INLINE gendiff #-}
gendiff :: Monad m => (a -> RefSeqView a) -> a -> S.ByteString m r -> Stream (Of Lump) m r
gendiff view = generic
  where
    isN        c = c == c2w 'N' || c == c2w 'n' || c == c2w '-'
    isNN (N2b a) = a > 3

    eq (N2b a) b = b == c2w 'Q' || b == c2w 'q' ||
                   "TCAG" `B.index` fromIntegral a == b ||
                   "tcag" `B.index` fromIntegral a == b

    generic !ref !smp = case view ref of
            NilRef                     -> lift $ S.effects smp
            l :== ref'                 -> ns l ref' (S.drop (fromIntegral l) smp)
            u :^  ref' | isNN u        -> ns 1 ref' (S.drop               1  smp)
            u :^  ref' ->
              lift (S.nextByte smp) >>= \case
                Left r                 -> pure r
                Right (x, smp')
                    | x .&. 0x7F <= 32 -> generic ref  smp'
                    | isN  x           -> ns    1 ref' smp'
                    | eq u x           -> eqs_2 1 ref' smp'
                    | otherwise        -> encVar u x `Q.cons` generic ref' smp'

    eqs_2 !n !ref !smp = case view ref of
            NilRef                     -> Eqs2 n `Q.cons` lift (S.effects smp)
            l :== ref'                 -> Eqs2 n `Q.cons` ns l ref' (S.drop (fromIntegral l) smp)
            u :^  ref' | isNN u        -> Eqs2 n `Q.cons` ns 1 ref' (S.drop               1  smp)
            u :^  ref'                 -> lift (S.nextByte smp) >>= \case
                Left r                 -> Q.cons (Eqs2 n) (pure r)
                Right (x, smp')
                    | x .&. 0x7F <= 32 -> eqs_2       n  ref  smp'
                    | eq u x           -> eqs_2 (succ n) ref' smp'
                    | isN  x           -> Q.cons (Eqs2 n) $ ns 1 ref' smp'
                    | otherwise        -> Q.cons (Eqs2 n) $ Q.cons (encVar u x) $ generic ref' smp'

    ns !n !ref !smp = case view ref of
            NilRef                     -> Ns n `Q.cons` lift (S.effects smp)
            l :== ref'                 -> ns (l+n) ref' (S.drop (fromIntegral l) smp)
            u :^  ref' | isNN u        -> ns (1+n) ref' (S.drop               1  smp)
            u :^  ref' ->
              lift (S.nextByte smp) >>= \case
                Left r                 -> Ns n `Q.cons` pure r
                Right (x, smp')
                    | x .&. 0x7F <= 32 -> ns       n  ref  smp'
                    | isN  x           -> ns (succ n) ref' smp'
                    | eq u x           -> Ns n `Q.cons` eqs_2 1 ref' smp'
                    | otherwise        -> Q.cons (Ns n) $ Q.cons (encVar u x) $ generic ref' smp'


-- two alleles in bits 0,1 and 2,3
fromAmbCode :: Word8 -> Word8
fromAmbCode c | c == c2w 't' =  0
              | c == c2w 'c' =  5
              | c == c2w 'a' = 10
              | c == c2w 'g' = 15

              | c == c2w 's' =  7
              | c == c2w 'w' =  2
              | c == c2w 'm' =  6
              | c == c2w 'k' =  3
              | c == c2w 'r' = 11
              | c == c2w 'y' =  1
              | otherwise    = unexpected ("fromAmbCode " ++ show c)

encVar :: Nuc2b -> Word8 -> Lump
encVar r c = encTwoNucs r $ fromAmbCode (c .|. 32)

encTwoNucs :: Nuc2b -> Word8 -> Lump
encTwoNucs (N2b r) ns = encTwoVars $ xor ns (r .|. shiftL r 2)

encTwoVars :: Word8 -> Lump
encTwoVars vs = fromMaybe (Ns 1) $ vv V.!? fromIntegral vs
  where
    !vv = V.fromList [ Eqs2 1,    RefTrans,    RefCompl,    RefTCompl
                     , RefTrans,  Trans2,      TransCompl,  TransTCompl
                     , RefCompl,  TransCompl,  Compl2,      ComplTCompl
                     , RefTCompl, TransTCompl, ComplTCompl, TCompl2

                     , Eqs1 1, Trans1, Compl1, TCompl1
                     , Eqs1 1, Trans1, Compl1, TCompl1
                     , Eqs1 1, Trans1, Compl1, TCompl1
                     , Eqs1 1, Trans1, Compl1, TCompl1 ]

diff2 :: Monad m => (() -> RefSeq) -> S.ByteString m r -> Stream (Of Lump) m r
diff2 = gendiff viewRS . ($ ())

diff :: Monad m => B.ByteString -> B.ByteString -> Stream (Of Lump) m ()
diff s1 = gendiff viewBS s1 . S.fromStrict
  where
    viewBS = maybe NilRef (\(a,b) -> fromCode (a .|. 32) :^ b) . B.uncons

    fromCode c | c == c2w 't' = N2b 0
               | c == c2w 'c' = N2b 1
               | c == c2w 'a' = N2b 2
               | c == c2w 'g' = N2b 3
               | otherwise =  N2b 255


data Frag a = Short !Char a | Long !L.ByteString a
  deriving Functor

patch :: Monad m => RefSeq -> Stream (Of Lump) m r -> Stream Frag m (Stream (Of Lump) m r)
patch ref = case unconsRS ref of
    Just (N2b hd,tl) -> lift . Q.next >=> \case
        Left  r ->          yields $ Long (unpackRS ref) (pure r)
        Right l -> case l of
            (Break   ,s) -> yields $ Long (unpackRS ref) s

            (Eqs2   n,s) -> wrap $ Long          (unpackRS $ takeRS n ref) $ patch (dropRS n ref) s
            (Eqs1   n,s) -> wrap $ Long          (unpackRS $ takeRS n ref) $ patch (dropRS n ref) s
            (Ns     n,s) -> wrap $ Long (L.replicate (fromIntegral n) 'N') $ patch (dropRS n ref) s

            (Del1   n,s) -> wrap $ Long (L.replicate (fromIntegral n) '-') $ patch (dropRS n ref) s
            (Del2   n,s) -> wrap $ Long (L.replicate (fromIntegral n) '-') $ patch (dropRS n ref) s
            (DelH   n,s) -> wrap $ Long          (unpackRS $ takeRS n ref) $ patch (dropRS n ref) s

            (Ins1   _,s) -> patch ref s
            (Ins2   _,s) -> patch ref s
            (InsH   _,s) -> patch ref s

            (Trans1     ,s) -> step "CTGA" s
            (Trans2     ,s) -> step "CTGA" s
            (Compl1     ,s) -> step "AGTC" s
            (Compl2     ,s) -> step "AGTC" s
            (TCompl1    ,s) -> step "GACT" s
            (TCompl2    ,s) -> step "GACT" s

            (RefTrans   ,s) -> step "YYRR" s
            (RefCompl   ,s) -> step "WSWS" s
            (TransCompl ,s) -> step "MKKM" s
            (RefTCompl  ,s) -> step "KMMK" s
            (TransTCompl,s) -> step "SWSW" s
            (ComplTCompl,s) -> step "RRYY" s
          where
            step cs s = wrap $ Short (if hd > 3 then 'N' else w2c $ B.index cs (fromIntegral hd)) (patch tl s)

    Nothing      -> clear
      where
        clear :: Monad m => Stream (Of Lump) m r -> Stream Frag m (Stream (Of Lump) m r)
        clear = lift . Q.next >=> \case
            Left          r  -> pure (pure r)
            Right (Break, a) -> pure a
            Right (_,     a) -> clear a

