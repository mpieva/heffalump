{-# LANGUAGE OverloadedLists, LambdaCase #-}
module Lump where

import Bio.Prelude                   hiding ( Ns )
import Data.ByteString.Builder
import Data.ByteString.Internal             ( c2w )
import Data.Hashable                        ( hash )
import Foreign.Marshal.Alloc                ( mallocBytes )
import Streaming
import Streaming.Prelude                    ( mapOf )
import System.Directory                     ( doesFileExist )

import qualified Data.ByteString            as B
import qualified Data.ByteString.Char8      as BC
import qualified Data.ByteString.Lazy       as L
import qualified Data.ByteString.Lazy.Char8 as C
import qualified Data.ByteString.Streaming  as S
import qualified Data.ByteString.Unsafe     as B
import qualified Data.Vector                as V
import qualified Data.Vector.Unboxed        as U

import Stretch ( Stretch, decode_dip, decode_hap, NucCode(..) )
import qualified Stretch  as S ( Stretch(..) )
import NewRef
import Util ( decomp )

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

data Lump a
    = Ns !Int a                         -- ^ uncalled stretch

    | Eqs1 !Int a                       -- ^ stretch of monoallelic matches
    | Trans1 a                          -- ^ monoallelic transversion
    | Compl1 a                          -- ^ monoallelic complement
    | TCompl1 a                         -- ^ monoallelic trans-complement

    | Eqs2 !Int a                       -- ^ stretch of diallelic matches
    | RefTrans a                        -- ^ nine more diallelic SNP codes (*groan*)
    | Trans2 a
    | RefCompl a
    | TransCompl a
    | Compl2 a
    | RefTCompl a
    | TransTCompl a
    | ComplTCompl a
    | TCompl2 a

    | Del1 !Int a                       -- ^ deletion with other allele unknown
    | Del2 !Int a                       -- ^ homozygous deletion
    | DelH !Int a                       -- ^ deletion, the other allele is the reference

    | Ins1 !Seq1 a                      -- ^ insertion with other allele unknown
    | Ins2 !Seq1 a                      -- ^ homozygous insertion
    | InsH !Seq1 a                      -- ^ insertion, the other allele is the reference

    | Break a                           -- ^ break marker (end-of-chromosome)
    | Done                              -- ^ end of data stream
  deriving (Functor, Show)

type Seq1 = U.Vector Word8

debugLump :: Fix Lump -> IO ()
debugLump = debugLump' 0 0

debugLump' :: Int -> Int -> Fix Lump -> IO ()
debugLump' c i fl = case unFix fl of
    Done          -> return ()
    Break       l -> do putStrLn $ shows (c,i) "\tBreak" ;           debugLump' (c+1) 0 l
    Ns        n l -> do putStrLn $ shows (c,i) "\tNs   " ++ show n ; debugLump' c (i+n) l
    Eqs1      n l -> do putStrLn $ shows (c,i) "\tEqs1 " ++ show n ; debugLump' c (i+n) l
    Eqs2      n l -> do putStrLn $ shows (c,i) "\tEqs2 " ++ show n ; debugLump' c (i+n) l

    Trans1      l -> do putStrLn $ shows (c,i) "\tTrans1 " ;      debugLump' c (i+1) l
    Compl1      l -> do putStrLn $ shows (c,i) "\tCompl1 " ;      debugLump' c (i+1) l
    TCompl1     l -> do putStrLn $ shows (c,i) "\tTCompl1 " ;     debugLump' c (i+1) l

    RefTrans    l -> do putStrLn $ shows (c,i) "\tRefTrans " ;    debugLump' c (i+1) l
    Trans2      l -> do putStrLn $ shows (c,i) "\tTrans2 " ;      debugLump' c (i+1) l
    RefCompl    l -> do putStrLn $ shows (c,i) "\tRefCompl " ;    debugLump' c (i+1) l
    TransCompl  l -> do putStrLn $ shows (c,i) "\tTransCompl " ;  debugLump' c (i+1) l
    Compl2      l -> do putStrLn $ shows (c,i) "\tCompl2 " ;      debugLump' c (i+1) l
    RefTCompl   l -> do putStrLn $ shows (c,i) "\tRefTCompl " ;   debugLump' c (i+1) l
    TransTCompl l -> do putStrLn $ shows (c,i) "\tTransTCompl " ; debugLump' c (i+1) l
    ComplTCompl l -> do putStrLn $ shows (c,i) "\tComplTCompl " ; debugLump' c (i+1) l
    TCompl2     l -> do putStrLn $ shows (c,i) "\tTCompl2 " ;     debugLump' c (i+1) l

    Del1      n l -> do putStrLn $ shows (c,i) "\tDel1 " ++ shows n " " ; debugLump' c i l
    Del2      n l -> do putStrLn $ shows (c,i) "\tDel2 " ++ shows n " " ; debugLump' c i l
    DelH      n l -> do putStrLn $ shows (c,i) "\tDelH " ++ shows n " " ; debugLump' c i l

    Ins1      s l -> do putStrLn $ shows (c,i) "\tIns1 " ++ shows s " " ; debugLump' c i l
    Ins2      s l -> do putStrLn $ shows (c,i) "\tIns2 " ++ shows s " " ; debugLump' c i l
    InsH      s l -> do putStrLn $ shows (c,i) "\tInsH " ++ shows s " " ; debugLump' c i l


{-# DEPRECATED normalizeLump "use normalizeLump'" #-}
normalizeLump :: Fix Lump -> Fix Lump
normalizeLump = ana (go . unFix)
  where
    go (Ns 0 l) = unFix l
    go (Ns n l) = case unFix l of
        Ns n' l' -> go $ Ns (n+n') l'
        _        -> Ns n l

    go (Eqs1 0 l) = unFix l
    go (Eqs1 n l) = case unFix l of
        Eqs1 n' l' -> go $ Eqs1 (n+n') l'
        _          -> Eqs1 n l

    go (Eqs2 0 l) = unFix l
    go (Eqs2 n l) = case unFix l of
        Eqs2 n' l' -> go $ Eqs2 (n+n') l'
        _          -> Eqs2 n l

    go (Del1 0 l) = unFix l
    go (Del1 n l) = case unFix l of
        Del1 n' l' -> go $ Del1 (n+n') l'
        _          -> Del1 n l

    go (Del2 0 l) = unFix l
    go (Del2 n l) = case unFix l of
        Del2 n' l' -> go $ Del2 (n+n') l'
        _          -> Del2 n l

    go (DelH 0 l) = unFix l
    go (DelH n l) = case unFix l of
        DelH n' l' -> go $ DelH (n+n') l'
        _          -> DelH n l

    go (Ins1 n l) | U.null n = unFix l
    go (Ins1 n l) = case unFix l of
        Ins1 n' l' -> go $ Ins1 (n U.++ n') l'
        _          -> Ins1 n l

    go (Ins2 n l) | U.null n = unFix l
    go (Ins2 n l) = case unFix l of
        Ins2 n' l' -> go $ Ins2 (n U.++ n') l'
        _          -> Ins2 n l

    go (InsH n l) | U.null n = unFix l
    go (InsH n l) = case unFix l of
        InsH n' l' -> go $ InsH (n U.++ n') l'
        _          -> InsH n l

    go lump = lump

normalizeLump' :: Monad m => Stream Lump m r -> Stream Lump m r
normalizeLump' = unfold (inspect >=> either (return . Left) go)
  where
    wrapE = either pure wrap

    go (Ns 0 l) = inspect l
    go (Ns n l) = inspect l >>= \case
        Right (Ns n' l') -> go $ Ns (n+n') l'
        l'               -> return $ Right (Ns n (wrapE l'))

    go (Eqs1 0 l) = inspect l
    go (Eqs1 n l) = inspect l >>= \case
        Right (Eqs1 n' l') -> go $ Eqs1 (n+n') l'
        l'                 -> return $ Right (Eqs1 n (wrapE l'))

    go (Eqs2 0 l) = inspect l
    go (Eqs2 n l) = inspect l >>= \case
        Right (Eqs2 n' l') -> go $ Eqs2 (n+n') l'
        l'                 -> return $ Right (Eqs2 n (wrapE l'))

    go (Del1 0 l) = inspect l
    go (Del1 n l) = inspect l >>= \case
        Right (Del1 n' l') -> go $ Del1 (n+n') l'
        l'                 -> return $ Right (Del1 n (wrapE l'))

    go (Del2 0 l) = inspect l
    go (Del2 n l) = inspect l >>= \case
        Right (Del2 n' l') -> go $ Del2 (n+n') l'
        l'                 -> return $ Right (Del2 n (wrapE l'))

    go (DelH 0 l) = inspect l
    go (DelH n l) = inspect l >>= \case
        Right (DelH n' l') -> go $ DelH (n+n') l'
        l'                 -> return $ Right (DelH n (wrapE l'))

    go (Ins1 s l) | U.null s = inspect l
    go (Ins1 s l) = inspect l >>= \case
        Right (Ins1 s' l') -> go $ Ins1 (s<>s') l'
        l'                 -> return $ Right (Ins1 s (wrapE l'))

    go (Ins2 s l) | U.null s = inspect l
    go (Ins2 s l) = inspect l >>= \case
        Right (Ins2 s' l') -> go $ Ins2 (s<>s') l'
        l'                 -> return $ Right (Ins2 s (wrapE l'))

    go (InsH s l) | U.null s = inspect l
    go (InsH s l) = inspect l >>= \case
        Right (InsH s' l') -> go $ InsH (s<>s') l'
        l'                 -> return $ Right (InsH s (wrapE l'))

    go lump = return $ Right lump

-- | Lowers all calls to haploid.  (Not sure
-- if this is really a good idea anymore...)
make_hap :: Monad m => Stream Lump m r -> Stream Lump m r
make_hap = maps step
  where
    step (Ns    n a) = Ns    n a
    step (Eqs1  n a) = Eqs1  n a
    step (Eqs2  n a) = Eqs1  n a
    step (Trans1  a) = Trans1  a
    step (Trans2  a) = Trans1  a
    step (Compl1  a) = Compl1  a
    step (Compl2  a) = Compl1  a
    step (TCompl1 a) = TCompl1 a
    step (TCompl2 a) = TCompl1 a

    step (RefTrans    a) = RefTrans    a
    step (RefCompl    a) = RefCompl    a
    step (TransCompl  a) = TransCompl  a
    step (RefTCompl   a) = RefTCompl   a
    step (TransTCompl a) = TransTCompl a
    step (ComplTCompl a) = ComplTCompl a

    step (Del1 n a) = Del1 n a
    step (Del2 n a) = Del1 n a
    step (DelH n a) = DelH n a

    step (Ins1 s a) = Ins1 s a
    step (Ins2 s a) = Ins1 s a
    step (InsH s a) = InsH s a

    step (Break  a) = Break  a


newtype PackedLump = PackLump { unpackLump :: L.ByteString }

-- | The encoding of an empty chromosome:  just the 'Break'
noLump :: PackedLump
noLump = PackLump (L.singleton 0)

{-# INLINE encodeLumpToMem #-}
encodeLumpToMem :: MonadIO m => Stream Lump m r -> m (Of PackedLump r)
encodeLumpToMem = fmap (mapOf PackLump) . S.toLazy . encodeLump . normalizeLump'

{-# INLINE encodeLump #-}
encodeLump :: MonadIO m => Stream Lump m r -> S.ByteString m r
encodeLump = go
  where
    go s = S.mwrap $ do p <- liftIO $ mallocBytes 0x8000
                        (o,r) <- fill p 0x8000 s 0
                        c <- liftIO $ B.unsafePackMallocCStringLen (castPtr p,0x8000)
                        return $ do S.chunk (B.take o c)
                                    either pure go r

    -- An insert is never longer than 255, so it needs no more than
    -- ceil(255/4)+2 bytes, hence the magic 66 below.  Once we have 66
    -- bytes, the next code is guaranteed to fit.
    fill p l s o
        | l - o < 66 = return (o, Right s)      -- need more space
        | otherwise  = inspect s >>= either (return . (,) o . Left) (\case
            Ns n s'                          -> stretchOf 0x40 n >>= fill p l s'
            Eqs1 n s'                        -> stretchOf 0x80 n >>= fill p l s'
            Trans1 s'                        -> wrd8 0x01 >>= fill p l s'
            Compl1 s'                        -> wrd8 0x02 >>= fill p l s'
            TCompl1 s'                       -> wrd8 0x03 >>= fill p l s'

            Eqs2 n s'                        -> stretchOf 0xC0 n >>= fill p l s'
            RefTrans s'                      -> wrd8 0x05 >>= fill p l s'
            Trans2 s'                        -> wrd8 0x06 >>= fill p l s'
            RefCompl s'                      -> wrd8 0x07 >>= fill p l s'
            TransCompl s'                    -> wrd8 0x08 >>= fill p l s'
            Compl2 s'                        -> wrd8 0x09 >>= fill p l s'
            RefTCompl s'                     -> wrd8 0x0A >>= fill p l s'
            TransTCompl s'                   -> wrd8 0x0B >>= fill p l s'
            ComplTCompl s'                   -> wrd8 0x0C >>= fill p l s'
            TCompl2 s'                       -> wrd8 0x0D >>= fill p l s'

            Del1 n s'                        -> indelOf 0x10 n >>= fill p l s'
            DelH n s'                        -> indelOf 0x20 n >>= fill p l s'
            Del2 n s'                        -> indelOf 0x30 n >>= fill p l s'

            Ins1 sq s'                       -> indelOf 0x18 (U.length sq) >>= seqOf sq >>= fill p l s'
            InsH sq s'                       -> indelOf 0x28 (U.length sq) >>= seqOf sq >>= fill p l s'
            Ins2 sq s'                       -> indelOf 0x38 (U.length sq) >>= seqOf sq >>= fill p l s'

            Break s'                         -> wrd8 0x00 >>= fill p l s')
    -- encode1  Done                              = mempty
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

decodeLump :: L.ByteString -> Fix Lump
decodeLump = ana decode1
  where
    decode1 s = case L.uncons s of
        Nothing                                -> Done
        Just (c,s1)
            | c == 0x00                        -> Break s1
            | c == 0x01                        -> Trans1 s1
            | c == 0x02                        -> Compl1 s1
            | c == 0x03                        -> TCompl1 s1
            | c == 0x04                        -> error $ "unexpected " ++ show c
            | c == 0x05                        -> RefTrans s1
            | c == 0x06                        -> Trans2 s1
            | c == 0x07                        -> RefCompl s1
            | c == 0x08                        -> TransCompl s1
            | c == 0x09                        -> Compl2 s1
            | c == 0x0A                        -> RefTCompl s1
            | c == 0x0B                        -> TransTCompl s1
            | c == 0x0C                        -> ComplTCompl s1
            | c == 0x0D                        -> TCompl2 s1
            | c == 0x0E                        -> error $ "unexpected " ++ show c
            | c == 0x0F                        -> error $ "unexpected " ++ show c

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
        | c .&. 0x3F == 0x3F  = cons (fromIntegral (L.index s1 0)) (L.drop 1 s1)
        | c .&. 0x3F == 0x3E  = cons (fromIntegral (L.index s1 0) .|.
                                      fromIntegral (L.index s1 1) `shiftL`  8) (L.drop 2 s1)
        | c .&. 0x3F == 0x3D  = cons (fromIntegral (L.index s1 0) .|.
                                      fromIntegral (L.index s1 1) `shiftL`  8 .|.
                                      fromIntegral (L.index s1 2) `shiftL` 16) (L.drop 3 s1)
        | c .&. 0x3F == 0x3C  = cons (fromIntegral (L.index s1 0) .|.
                                      fromIntegral (L.index s1 1) `shiftL`  8 .|.
                                      fromIntegral (L.index s1 2) `shiftL` 16 .|.
                                      fromIntegral (L.index s1 3) `shiftL` 24) (L.drop 4 s1)
        | otherwise           = cons (fromIntegral (c .&. 0x3F)) s1

    del_of cons c s1
        | c .&. 0x07 == 0  = cons (fromIntegral (L.head s1)) (L.tail s1)
        | otherwise        = cons (fromIntegral c .&. 0x07) s1

    ins_of cons c s1
        | c .&. 0x07 == 0  = seq_of cons U.empty (L.head s1) (L.tail s1)
        | otherwise        = seq_of cons U.empty (c .&. 0x07) s1



    seq_of cons acc 0 s1 = cons acc s1

    seq_of cons acc 1 s1 = cons acc' (L.tail s1)
      where
        acc' = acc `U.snoc` (L.head s1 `shiftR` 0 .&. 0x3)

    seq_of cons acc 2 s1 = cons acc' (L.tail s1)
      where
        acc' = acc U.++ [ L.head s1 `shiftR` 0 .&. 0x3
                        , L.head s1 `shiftR` 2 .&. 0x3 ]

    seq_of cons acc 3 s1 = cons acc' (L.tail s1)
      where
        acc' = acc U.++ [ L.head s1 `shiftR` 0 .&. 0x3
                        , L.head s1 `shiftR` 2 .&. 0x3
                        , L.head s1 `shiftR` 4 .&. 0x3 ]

    seq_of cons acc n s1 = seq_of cons acc' (n-4) (L.tail s1)
      where
        acc' = acc U.++ [ L.head s1 `shiftR` 0 .&. 0x3
                        , L.head s1 `shiftR` 2 .&. 0x3
                        , L.head s1 `shiftR` 4 .&. 0x3
                        , L.head s1 `shiftR` 6 .&. 0x3 ]


-- | Merging without access to a reference sequence.  This code doesn't
-- believe in Indels and skips over them.
--
-- 'noutgroups':  number of outgroups

mergeLumps :: Int -> V.Vector (Fix Lump) -> [[Variant]]
mergeLumps !noutgroups = filter (not . null) . go 0 0
    -- Merging stretches.  We define a 'Variant' as anything that is
    -- different from the reference.  Therefore, 'Eqs' ('Eqs1') and 'Ns'
    -- never create a 'Variant' and we can skip forwards.  A 'Done' is
    -- an error.  Else, if at least one sample has some sort of call, we
    -- need to collect the alleles, produce up to four variants, then
    -- output and skip forward by one base.
  where
    go :: Int -> Int -> V.Vector (Fix Lump) -> [[Variant]]
    go !ix !pos !smps
        -- all samples 'Done', no more 'Variant's to produce
        | V.all (isDone . unFix) smps = []

        -- all samples 'Break' or are done, next chromosome
        | V.all (isBreak . unFix) smps = go (succ ix) 0 (V.map (skipBreak . unFix) smps)

        -- We ignore outgroups in computing the longest stretch.  That
        -- way, outgroups don't get to define variants, but participate
        -- in variants found in another way.
        | otherwise = case longestStretch $ V.drop noutgroups smps of

            -- no stretch, have to look for vars
            0 -> mkVar ix    pos                                            smps $
                 go ix (succ pos) (unS $ V.mapM (S . skipStretch 1 . unFix) smps)

            -- a stretch of l non-variants can be skipped over
            l -> go ix (pos + fromIntegral l) (unS $ V.mapM (S . skipStretch l . unFix) smps)


    isDone Done = True
    isDone    _ = False

    isBreak  Done     = True
    isBreak (Break _) = True
    isBreak        _  = False

    skipBreak (Break l) =     l
    skipBreak        l  = Fix l

    -- Skip over exactly l sites.
    skipStretch :: Int -> Lump (Fix Lump) -> Fix Lump
    skipStretch !l _ | l < 0  = error "WTF?!"

    skipStretch !l (Del1 _ s) = skipStretch l (unFix s)
    skipStretch !l (Del2 _ s) = skipStretch l (unFix s)
    skipStretch !l (DelH _ s) = skipStretch l (unFix s)

    skipStretch !l (Ins1 _ s) = skipStretch l (unFix s)
    skipStretch !l (Ins2 _ s) = skipStretch l (unFix s)
    skipStretch !l (InsH _ s) = skipStretch l (unFix s)

    skipStretch _ (Break    s) = Fix $ Break s
    skipStretch _  Done        = Fix $ Done

    skipStretch !l (Ns   !n s) | l == n = s
                               | l <  n = Fix $ Ns (n-l) s
                               | otherwise = skipStretch (l-n) (unFix s)
    skipStretch !l (Eqs1 !n s) | l == n = s
                               | l <  n = Fix $ Eqs1 (n-l) s
                               | otherwise = skipStretch (l-n) (unFix s)
    skipStretch !l (Eqs2 !n s) | l == n = s
                               | l <  n = Fix $ Eqs2 (n-l) s
                               | otherwise = skipStretch (l-n) (unFix s)

    skipStretch  0          s  = Fix s
    skipStretch !l (Trans1  s) = skipStretch (l-1) (unFix s)
    skipStretch !l (Compl1  s) = skipStretch (l-1) (unFix s)
    skipStretch !l (TCompl1 s) = skipStretch (l-1) (unFix s)

    skipStretch !l (RefTrans    s) = skipStretch (l-1) (unFix s)
    skipStretch !l (Trans2      s) = skipStretch (l-1) (unFix s)
    skipStretch !l (RefCompl    s) = skipStretch (l-1) (unFix s)
    skipStretch !l (TransCompl  s) = skipStretch (l-1) (unFix s)
    skipStretch !l (Compl2      s) = skipStretch (l-1) (unFix s)
    skipStretch !l (RefTCompl   s) = skipStretch (l-1) (unFix s)
    skipStretch !l (TransTCompl s) = skipStretch (l-1) (unFix s)
    skipStretch !l (ComplTCompl s) = skipStretch (l-1) (unFix s)
    skipStretch !l (TCompl2     s) = skipStretch (l-1) (unFix s)


    longestStretch :: V.Vector (Fix Lump) -> Int
    longestStretch = V.minimum . V.map (get_stretch_len . unFix)
      where
        get_stretch_len (Ns   n _) = n
        get_stretch_len (Eqs1 n _) = n
        get_stretch_len (Eqs2 n _) = n
        get_stretch_len (Break  _) = maxBound
        get_stretch_len  Done      = maxBound
        get_stretch_len  _         = 0


    -- Find all the variants, anchored on the reference allele, and
    -- split them.  Misfitting alleles are not counted.
    mkVar :: Int -> Int -> V.Vector (Fix Lump) -> [[Variant]] -> [[Variant]]
    mkVar ix pos ss = (:)
            [ Variant ix pos ref_indet (V2b alt) calls
            | (alt, ct) <- zip [1..3] [ct_trans, ct_compl, ct_tcompl]
            , let calls = V.convert $ V.map (ct . unFix) ss
            -- it's only a variant if at least one alt called
            , U.any (\c -> c .&. 0xC /= 0) (U.drop noutgroups calls) ]

    ref_indet = N2b 255

    -- Variant codes:  #ref + 4 * #alt
    ct_trans :: Lump a -> Word8
    ct_trans (Eqs1      _ _) = 1
    ct_trans (Trans1      _) = 4

    ct_trans (Eqs2      _ _) = 2
    ct_trans (RefTrans    _) = 5
    ct_trans (Trans2      _) = 8
    ct_trans              _  = 0

    ct_compl :: Lump a -> Word8
    ct_compl (Eqs1      _ _) = 1
    ct_compl (Compl1      _) = 4

    ct_compl (Eqs2      _ _) = 2
    ct_compl (RefCompl    _) = 5
    ct_compl (Compl2      _) = 8
    ct_compl              _  = 0

    ct_tcompl :: Lump a -> Word8
    ct_tcompl (Eqs1      _ _) = 1
    ct_tcompl (TCompl1     _) = 4

    ct_tcompl (Eqs2      _ _) = 2
    ct_tcompl (RefTCompl   _) = 5
    ct_tcompl (TCompl2     _) = 8
    ct_tcompl              _  = 0


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
-- 'NewRefSeq', too.
stretchToLump :: NewRefSeqs -> Stretch -> Fix Lump
stretchToLump nrs = normalizeLump . go1 (nrss_seqs nrs)
  where
    go1 [     ] _ = Fix Done
    go1 (r0:rs) l = go (r0 ()) l
      where
        go r (S.Chrs a b k) = call a (call b (flip go k)) r
        go r (S.Ns     c k) = Fix $ Ns   (fromIntegral $ c+c) $ go (dropNRS (fromIntegral $ c+c) r) k
        go r (S.Eqs    c k) = Fix $ Eqs2 (fromIntegral $ c+c) $ go (dropNRS (fromIntegral $ c+c) r) k
        go r (S.Eqs1   c k) = Fix $ Eqs1 (fromIntegral $ c+c) $ go (dropNRS (fromIntegral $ c+c) r) k
        go _ (S.Break    k) = Fix $ Break                     $ go1 rs k
        go _  S.Done        = Fix   Done

    call :: NucCode -> (NewRefSeq -> Fix Lump) -> NewRefSeq -> Fix Lump
    call (NucCode b) k r = case unconsNRS r of
        Nothing     -> Fix Done
        Just (a,r') -> Fix . encTwoNucs a (nc2vs U.! fromIntegral b) $ k r'

    nc2vs :: U.Vector Word8
    !nc2vs = [ 255,  5, 10, 15, 0,        -- NACGT
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
decode :: Either String NewRefSeqs -> L.ByteString -> Fix Lump
decode (Left err)  str | "HEF\3" `L.isPrefixOf` str = decodeLump . L.drop 5 $ L.dropWhile (/= 0) str
                       | otherwise                  = error err

decode (Right  r) str | "HEF\0"  `L.isPrefixOf` str = stretchToLump r $ decode_dip (L.drop 4 str)
                      | "HEF\1"  `L.isPrefixOf` str = stretchToLump r $ decode_hap (L.drop 4 str)
                      | "HEF\3"  `L.isPrefixOf` str = go (L.drop 4 str)
                      | 176:113:2:_ <- L.unpack str = stretchToLump r $ decode_dip str
                      | 177:113:2:_ <- L.unpack str = stretchToLump r $ decode_dip str
                      | otherwise                   = error "File format not recognized?"
  where
    go s = let (hs,s') = L.splitAt 4 . L.drop 1 . L.dropWhile (/= 0) $ s
               hv = L.foldr (\b a -> fromIntegral b .|. shiftL a 8) 0 hs :: Word32
           in if hv == fromIntegral (hash (nrss_chroms r, nrss_lengths r))
              then decodeLump s'
              else error "Incompatible reference genome."


getRefPath :: L.ByteString -> Maybe FilePath
getRefPath str | "HEF\3" `L.isPrefixOf` str = Just . C.unpack . L.takeWhile (/= 0) $ L.drop 4 str
               | otherwise                  = Nothing

decodeMany :: Maybe FilePath -> [FilePath] -> IO ( Either String NewRefSeqs, V.Vector (Fix Lump) )
decodeMany mrs fs = do
    raws <- mapM (fmap decomp . L.readFile) fs
    rs <- case mrs of
            Just fp -> Right <$> readTwoBit fp
            Nothing -> do
                fps <- filterM doesFileExist $ mapMaybe getRefPath raws
                case fps of
                    fp : _ -> Right <$> readTwoBit fp
                    [    ] -> return . Left $ "No reference found.  Looked for it at "
                                ++ intercalate ", " (mapMaybe getRefPath raws) ++ "."
    return ( rs, V.fromList $ map (decode rs) raws )



-- | Encode a 'Lump' and enough information about the 'NewRefSeqs' to be
-- (1) able to find it again and (2) to make sure we got the right one
-- when operating on a 'Lump'.
encode :: MonadIO m => NewRefSeqs -> Stream Lump m r -> S.ByteString m r
encode r s = do S.toStreamingByteString (encodeHeader r)
                encodeLump s

encodeHeader :: NewRefSeqs -> Builder
encodeHeader r = byteString "HEF\3" <>
                 byteString (nrss_path r) <> word8 0 <>
                 int32LE (fromIntegral (hash (nrss_chroms r, nrss_lengths r)))

-- | Diff between a reference and a sample in string format.  The sample
-- is treated as diploid.
gendiff :: (a -> RefSeqView a) -> a -> L.ByteString -> Fix Lump -> Fix Lump
gendiff view = generic
  where
    isN        c = c == c2w 'N' || c == c2w 'n' || c == c2w '-'
    eq (N2b a) b = b == c2w 'Q' || b == c2w 'q' ||
                   "TCAG" `B.index` fromIntegral a == b ||
                   "tcag" `B.index` fromIntegral a == b

    generic !ref !smp =
        case (L.uncons smp, view ref) of
            (Nothing, _)    -> id
            (_, NilRef)     -> id
            (_, l :== ref') -> Fix . Ns l . generic ref' (L.drop (fromIntegral l) smp)
            (Just (x, smp'), u :^ ref')
                | x .&. 0x7F <= 32             -> generic ref  smp'
                | isN  x    -> Fix . Ns       1 . generic ref' smp'
                | eq u x    -> Fix . Eqs2     1 . generic ref' smp'
                | otherwise -> Fix . encVar u x . generic ref' smp'


-- | Generic diff between two aligned strings using "streaming" code.
--
-- XXX  We should deal with gaps, now that we can actually encode them
{-# INLINE gendiff' #-}
gendiff' :: Monad m => (a -> RefSeqView a) -> a -> S.ByteString m r -> Stream Lump m r
gendiff' view = generic
  where
    isN        c = c == c2w 'N' || c == c2w 'n' || c == c2w '-'
    isNN (N2b a) = a > 3

    eq (N2b a) b = b == c2w 'Q' || b == c2w 'q' ||
                   "TCAG" `B.index` fromIntegral a == b ||
                   "tcag" `B.index` fromIntegral a == b

    generic !ref !smp = effect $ case view ref of
            NilRef              -> pure <$> S.effects smp
            l :== ref'          -> ns l ref' (S.drop (fromIntegral l) smp)
            u :^  ref' | isNN u -> ns 1 ref' (S.drop               1  smp)
            u :^  ref' -> S.nextByte smp >>= \case
                Left r                 -> return $ pure r
                Right (x, smp')
                    | x .&. 0x7F <= 32 -> return $ generic ref  smp'
                    | isN  x           -> ns    1 ref' smp'
                    | eq u x           -> eqs_2 1 ref' smp'
                    | otherwise        -> return $ yields (encVar u x smp') >>= generic ref'

    eqs_2 !n !ref !smp = case view ref of
            NilRef              -> yields . Eqs2 n <$> S.effects smp
            l :== ref'          -> return $ yields (Eqs2 n (S.drop (fromIntegral l) smp)) >>= effect . ns l ref'
            u :^  ref' | isNN u -> return $ yields (Eqs2 n (S.drop               1  smp)) >>= effect . ns 1 ref'
            u :^  ref'                 -> S.nextByte smp >>= \case
                Left r                 -> return $ yields (Eqs2 n r)
                Right (x, smp')
                    | x .&. 0x7F <= 32 -> eqs_2       n  ref  smp'
                    | eq u x           -> eqs_2 (succ n) ref' smp'
                    | isN  x           -> return $ yields (Eqs2 n smp') >>= effect . ns 1 ref'
                    | otherwise        -> return $ yields (Eqs2 n smp') >>= yields . encVar u x >>= generic ref'

    ns !n !ref !smp = case view ref of
            NilRef              -> yields . Ns n <$> S.effects smp
            l :== ref'          -> ns (l+n) ref' (S.drop (fromIntegral l) smp)
            u :^  ref' | isNN u -> ns (1+n) ref' (S.drop               1  smp)
            u :^  ref' -> S.nextByte smp >>= \case
                Left r                 -> return $ yields (Ns n r)
                Right (x, smp')
                    | x .&. 0x7F <= 32 -> ns       n  ref  smp'
                    | isN  x           -> ns (succ n) ref' smp'
                    | eq u x           -> return $ yields (Ns n smp') >>= effect . eqs_2 1 ref'
                    | otherwise        -> return $ yields (Ns n smp') >>= yields . encVar u x >>= generic ref'


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
              | otherwise    = error $ "[fromAmbCode] What?! " ++ show c

encVar :: Nuc2b -> Word8 -> a -> Lump a
encVar r c = encTwoNucs r $ fromAmbCode (c .|. 32)

encTwoNucs :: Nuc2b -> Word8 -> a -> Lump a
encTwoNucs (N2b r) ns = encTwoVars $ xor ns (r .|. shiftL r 2)

encTwoVars :: Word8 -> a -> Lump a
encTwoVars vs = fromMaybe (Ns 1) $ vv V.!? fromIntegral vs
  where
    !vv = [ Eqs2 1,    RefTrans,    RefCompl,    RefTCompl
          , RefTrans,  Trans2,      TransCompl,  TransTCompl
          , RefCompl,  TransCompl,  Compl2,      ComplTCompl
          , RefTCompl, TransTCompl, ComplTCompl, TCompl2

          , Eqs1 1, Trans1, Compl1, TCompl1
          , Eqs1 1, Trans1, Compl1, TCompl1
          , Eqs1 1, Trans1, Compl1, TCompl1
          , Eqs1 1, Trans1, Compl1, TCompl1 ]


diff2 :: (() -> NewRefSeq) -> L.ByteString -> Fix Lump -> Fix Lump
diff2 = gendiff viewNRS . ($ ())

diff2' :: Monad m => (() -> NewRefSeq) -> S.ByteString m r -> Stream Lump m r
diff2' = gendiff' viewNRS . ($ ())

diff' :: Monad m => L.ByteString -> S.ByteString m r -> Stream Lump m r
diff' = gendiff' viewLBS
  where
    viewLBS = maybe NilRef (\(a,b) -> fromCode (a .|. 32) :^ b) . L.uncons

    fromCode c | c == c2w 't' = N2b 0
               | c == c2w 'c' = N2b 1
               | c == c2w 'a' = N2b 2
               | c == c2w 'g' = N2b 3
               | otherwise =  N2b 255

diff :: L.ByteString -> L.ByteString -> Fix Lump -> Fix Lump
diff = gendiff viewLBS
  where
    viewLBS = maybe NilRef (\(a,b) -> fromCode (a .|. 32) :^ b) . L.uncons

    fromCode c | c == c2w 't' = N2b 0
               | c == c2w 'c' = N2b 1
               | c == c2w 'a' = N2b 2
               | c == c2w 'g' = N2b 3
               | otherwise =  N2b 255


data Frag = Short !Char Frag | Long !L.ByteString Frag | Term (Fix Lump)

patch :: NewRefSeq -> Fix Lump -> Frag
patch ref (Fix l) = case unconsNRS ref of
    Just (N2b hd,tl) -> case l of
        Break    s -> Long (unpackNRS ref) $ Term s
        Done       -> Long (unpackNRS ref) $ Term $ Fix Done

        Eqs2   n s -> Long        (unpackNRS $ takeNRS n ref) $ patch (dropNRS n ref) s
        Eqs1   n s -> Long        (unpackNRS $ takeNRS n ref) $ patch (dropNRS n ref) s
        Ns     n s -> Long (C.replicate (fromIntegral n) 'N') $ patch (dropNRS n ref) s

        Del1   n s -> Long (C.replicate (fromIntegral n) '-') $ patch (dropNRS n ref) s
        Del2   n s -> Long (C.replicate (fromIntegral n) '-') $ patch (dropNRS n ref) s
        DelH   n s -> Long        (unpackNRS $ takeNRS n ref) $ patch (dropNRS n ref) s

        Ins1   _ s -> patch ref s
        Ins2   _ s -> patch ref s
        InsH   _ s -> patch ref s

        Trans1      s -> step "CTGA" s
        Trans2      s -> step "CTGA" s
        Compl1      s -> step "AGTC" s
        Compl2      s -> step "AGTC" s
        TCompl1     s -> step "GACT" s
        TCompl2     s -> step "GACT" s

        RefTrans    s -> step "YYRR" s
        RefCompl    s -> step "WSWS" s
        TransCompl  s -> step "MKKM" s
        RefTCompl   s -> step "KMMK" s
        TransTCompl s -> step "SWSW" s
        ComplTCompl s -> step "RRYY" s
      where
        step cs s = Short (if hd > 3 then 'N' else BC.index cs (fromIntegral hd)) (patch tl s)

    Nothing      -> clear (Fix l)
      where
        clear :: Fix Lump -> Frag
        clear (Fix m) = case m of
            Ns        _ a -> clear a
            Eqs1      _ a -> clear a
            Eqs2      _ a -> clear a

            Trans1      a -> clear a
            Compl1      a -> clear a
            TCompl1     a -> clear a

            RefTrans    a -> clear a
            Trans2      a -> clear a
            RefCompl    a -> clear a
            TransCompl  a -> clear a
            Compl2      a -> clear a
            RefTCompl   a -> clear a
            TransTCompl a -> clear a
            ComplTCompl a -> clear a
            TCompl2     a -> clear a

            Del1      _ a -> clear a
            Del2      _ a -> clear a
            DelH      _ a -> clear a

            Ins1      _ a -> clear a
            Ins2      _ a -> clear a
            InsH      _ a -> clear a

            Break       a -> Term a
            Done          -> Term (Fix Done)


-- XXX  Stolen from data-fix.  If we migrate completely to streaming,
-- this isn't needed anymore.
newtype Fix f = Fix { unFix :: f (Fix f) }

ana :: Functor f => (a -> f a) -> (a -> Fix f)
ana f = Fix . fmap (ana f) . f
