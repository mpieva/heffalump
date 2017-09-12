{-# LANGUAGE LambdaCase #-}
module Util where

import Bio.Prelude
import Streaming
import System.Console.GetOpt
import System.IO                        ( hPutStrLn, stderr )

import qualified Data.ByteString                 as B
import qualified Data.ByteString.Lazy            as L
import qualified Data.ByteString.Lazy.Internal   as L ( ByteString(..) )
import qualified Data.ByteString.Streaming       as SB ( nextByte, cons' )
import qualified Data.ByteString.Streaming.Char8 as S
import qualified Codec.Compression.Zlib.Internal as Z

decomp :: L.ByteString -> L.ByteString
decomp s0 = case L.uncons s0 of
    Just (31, s') -> case L.uncons s' of
        Just (139,_) -> decomp_loop s0
        _            -> s0 -- no decompression needed
    _                -> s0 -- no decompression needed
  where
    decomp_loop :: L.ByteString -> L.ByteString
    decomp_loop s = case L.uncons s of
        Just (31, s') -> case L.uncons s' of
            Just (139,_) -> Z.foldDecompressStreamWithInput L.Chunk decomp_loop throw
                            (Z.decompressST Z.gzipOrZlibFormat Z.defaultDecompressParams) s
            _            -> L.empty  -- ignores trailing garbage
        _                -> L.empty  -- ignores trailing garbage

-- | Checks if the input is GZip at all, returns it unchanged if it
-- isn't.  Else runs gunzip and discards the remaining stream.  (Some
-- HETFA files appear to have junk at the end...)
decomp' :: MonadIO m => S.ByteString m r -> S.ByteString m r
decomp' s0 = S.mwrap $ SB.nextByte s0 >>= \case
    Right (31, s') -> SB.nextByte s' >>= \case
        Right (139,s'') -> return . S.drained $ gunzip (SB.cons' 31 (SB.cons' 139 s''))
        Right ( c, s'') -> return $ SB.cons' 31 (SB.cons' c s'')
        Left     r      -> return $ SB.cons' 31 (pure r)
    Right ( c, s') -> return $ SB.cons' c s'
    Left     r     -> return $ pure r

gunzip :: MonadIO m => S.ByteString m r -> S.ByteString m (S.ByteString m r)
gunzip = go $ Z.decompressIO Z.gzipOrZlibFormat Z.defaultDecompressParams
  where
    -- get next chunk, make sure it is empty iff the input ended
    go (Z.DecompressInputRequired next) inp =
        lift (S.nextChunk inp) >>= \case
            Left r          -> liftIO (next B.empty) >>= flip go (pure r)
            Right (ck,inp')
                | B.null ck -> go (Z.DecompressInputRequired next) inp'
                | otherwise -> liftIO (next ck) >>= flip go inp'

    go (Z.DecompressOutputAvailable outchunk next) inp =
        S.chunk outchunk >> liftIO next >>= flip go inp

    go (Z.DecompressStreamEnd inchunk) inp =
        pure (S.chunk inchunk >> inp)

    go (Z.DecompressStreamError derr) _inp =
        liftIO (throwM derr)


mk_opts :: String -> String -> [OptDescr (b -> IO b)] -> [OptDescr (b -> IO b)]
mk_opts cmd moreopts ods = ods'
  where
    ods' = ods ++ [ Option "h?" ["help","usage"] (NoArg usage) "Display this usage information" ]

    usage _ = do pn <- getProgName
                 hPutStrLn stderr $ usageInfo
                    ("Usage: " ++ pn ++ ' ': cmd ++ " [options...] " ++ moreopts) ods
                 exitSuccess

parseOpts :: Bool -> a -> [OptDescr (a -> IO a)] -> [String] -> IO ([String], a)
parseOpts fileargs def ods args = do
    let (opts, files, errs) = getOpt Permute ods args
    unless (null errs) $ do mapM_ (hPutStrLn stderr) errs >> exitFailure
    unless (fileargs || null files) $ do
        mapM_ (hPutStrLn stderr . (++) "unexpected argument " . show) files >> exitFailure
    (,) files <$> foldl (>>=) (return def) opts


data FastaSeq m r = FastaSeq !B.ByteString (S.ByteString m r) deriving Functor

readSampleFa :: MonadResource m => FilePath -> Stream (FastaSeq m) m ()
readSampleFa = parseFasta . decomp' . S.readFile

-- | Parsing a FastA file results in pairs of sequence names and
-- sequences.  The sequences are still text with their case preserved.
-- This assumes that the character '>' does not show up anywhere except
-- to introduce a header and that downstream can deal with whitespace
-- (especially line breaks).  Annoying as the use of
-- streaming-bytestring is here, it completely prevents otherwise
-- unmanageable memory leaks.

parseFasta :: Monad m => S.ByteString m r -> Stream (FastaSeq m) m r
parseFasta = go . ignoreBody . S.lines
  where
    go :: Monad m => Stream (S.ByteString m) m r -> Stream (FastaSeq m) m r
    go s = effect (inspect s >>= go')

    go' :: Monad m => Either r (S.ByteString m (Stream (S.ByteString m) m r)) -> m (Stream (FastaSeq m) m r)
    go' (Left    r) = return $ pure r
    go' (Right hdr) = do
                h :> ls <- S.toStrict hdr
                let name = B.drop 1 $ B.takeWhile (not . is_space_w8) h
                return $ yields (FastaSeq name (S.concat (takeBody ls))) >>= go

    is_space_w8 :: Word8 -> Bool
    is_space_w8 x = (x .&. 0x7F) <= 32

    -- Take lines as long as they don't start with ">"
    takeBody :: Monad m => Stream (S.ByteString m) m r -> Stream (S.ByteString m) m (Stream (S.ByteString m) m r)
    takeBody s = effect $ inspect s >>= \case
        -- empty input stream:  return an empty stream that returns an
        -- empty stream.
        Left r -> return $ pure (pure r)     -- ?!
        -- inspect first line...
        Right line -> S.uncons line >>= \case
            -- empty line:  return it and continue
            Left      s'      ->  return $ yields (pure s') >>= takeBody

            -- header line:  return empty stream that returns an empty
            -- stream that starts with the restored line.
            Right ('>',line') -> return $ pure (wrap (S.cons '>' line'))

            -- body line:  return a stream that begins with it, recurse
            -- on its tail
            Right ( c ,line') -> return $ yields (S.cons c line') >>= takeBody

    ignoreBody :: Monad m => Stream (S.ByteString m) m r -> Stream (S.ByteString m) m r
    ignoreBody s = effect $ inspect s >>= \case
        -- empty input stream:  return an empty stream
        Left r -> return $ pure r
        -- inspect first line...
        Right line -> S.uncons line >>= \case
            -- empty line:  continue
            Left      s'      ->  return $ ignoreBody s'

            -- header line:
            -- return a stream that starts with the restored line.
            Right ('>',line') -> return $ wrap (S.cons '>' line')

            -- body line:  continue
            Right ( _ ,line') -> return $ ignoreBody $ effect (S.effects line')


readNumIO :: String -> IO Int
readNumIO s = case reads s of
    [(n,[ ])] -> return n
    [(n,"k")] -> return $ n * 1000
    [(n,"M")] -> return $ n * 1000000
    [(n,"G")] -> return $ n * 1000000000
    _         -> fail $ "unable to parse: " ++ show s


