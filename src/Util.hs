module Util
    ( decomp
    , gunzip
    , gzip
    , mk_opts
    , parseOpts
    , parseFileOpts
    , parseFasta
    , readSampleFa
    , readNumIO
    , FastaSeq(..)
    , unexpected
    ) where

import Bio.Prelude
import Streaming
import System.Console.GetOpt
import System.IO                        ( hPutStrLn, stderr )

import qualified Data.ByteString                 as B
import qualified Data.ByteString.Streaming       as Q ( nextByte, cons' )
import qualified Data.ByteString.Streaming.Char8 as S
import qualified Codec.Compression.Zlib.Internal as Z

-- | Checks if the input is GZip at all, returns it unchanged if it
-- isn't.  Else runs gunzip.
decomp :: MonadIO m => S.ByteString m r -> S.ByteString m r
decomp = decompWith return

-- | Checks if the input is GZip at all, and runs gunzip if it is.  If
-- it isn't, it runs 'k' on the input.
decompWith :: MonadIO m
           => (S.ByteString m r -> m (S.ByteString m r))
           -> S.ByteString m r -> S.ByteString m r
decompWith k s0 = S.mwrap $ Q.nextByte s0 >>= \case
    Right (31, s') -> Q.nextByte s' >>= \case
        Right (139,s'') -> return $ gunzip (Q.cons' 31 (Q.cons' 139 s''))
        Right ( c, s'') -> return $ Q.cons' 31 (Q.cons' c s'')
        Left     r      -> return $ Q.cons' 31 (pure r)
    Right ( c, s') -> k $ Q.cons' c s'
    Left     r     -> k $ pure r

-- | Decompresses a gzip stream.  If the leftovers look like another
-- gzip stream, it recurses (some files look like gzip and contain
-- bgzip).  Otherwise, the leftover are discarded (some HETFA files
-- appear to have junk at the end).
gunzip :: MonadIO m => S.ByteString m r -> S.ByteString m r
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
        -- decompress leftovers if possible, else discard
        decompWith (fmap pure . S.effects) (S.chunk inchunk >> inp)

    go (Z.DecompressStreamError derr) _inp =
        liftIO (throwM derr)

gzip :: MonadIO m => S.ByteString m r -> S.ByteString m r
gzip = go $ Z.compressIO Z.gzipFormat Z.defaultCompressParams
  where
    -- get next chunk, make sure it is empty iff the input ended
    go (Z.CompressInputRequired next) inp =
        lift (S.nextChunk inp) >>= \case
            Left r          -> liftIO (next B.empty) >>= flip go (pure r)
            Right (ck,inp')
                | B.null ck -> go (Z.CompressInputRequired next) inp'
                | otherwise -> liftIO (next ck) >>= flip go inp'

    go (Z.CompressOutputAvailable outchunk next) inp =
        S.chunk outchunk >> liftIO next >>= flip go inp

    go Z.CompressStreamEnd inp = lift (S.effects inp)


mk_opts :: String -> String -> [OptDescr (b -> IO b)] -> [OptDescr (b -> IO b)]
mk_opts cmd moreopts ods = ods'
  where
    ods' = ods ++ [ Option "h?" ["help","usage"] (NoArg usage) "Display this usage information" ]

    usage _ = do pn <- getProgName
                 hPutStrLn stderr $ usageInfo
                    ("Usage: " ++ pn ++ ' ': cmd ++ " [options...] " ++ moreopts) ods
                 exitSuccess

parseOpts :: a -> [OptDescr (a -> IO a)] -> [String] -> IO a
parseOpts def ods args = do
    let (opts, files, errs) = getOpt Permute ods args
    unless (null errs) $ do mapM_ (hPutStrLn stderr) errs >> exitFailure
    unless (null files) $ do
        mapM_ (hPutStrLn stderr . (++) "unexpected argument " . show) files >> exitFailure
    foldl (>>=) (return def) opts

parseFileOpts :: a -> [OptDescr (a -> IO a)] -> [String] -> IO ([String], a)
parseFileOpts def ods args = do
    let (opts, files, errs) = getOpt Permute ods args
    unless (null errs) $ do mapM_ (hPutStrLn stderr) errs >> exitFailure
    (,) files <$> foldl (>>=) (return def) opts


data FastaSeq m r = FastaSeq !B.ByteString (S.ByteString m r) deriving Functor

readSampleFa :: MonadResource m => FilePath -> Stream (FastaSeq m) m ()
readSampleFa = parseFasta . decomp . S.readFile

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
    go = lift . inspect >=> go'

    go' :: Monad m => Either r (S.ByteString m (Stream (S.ByteString m) m r)) -> Stream (FastaSeq m) m r
    go' (Left    r) = pure r
    go' (Right hdr) = do
                h :> ls <- lift $ S.toStrict hdr
                let name = B.drop 1 $ B.takeWhile (not . is_space_w8) h
                yields (FastaSeq name (S.concat (takeBody ls))) >>= go

    is_space_w8 :: Word8 -> Bool
    is_space_w8 x = (x .&. 0x7F) <= 32

    -- Take lines as long as they don't start with ">"
    takeBody :: Monad m => Stream (S.ByteString m) m r -> Stream (S.ByteString m) m (Stream (S.ByteString m) m r)
    takeBody = lift . inspect >=> \case
        -- empty input stream:  return an empty stream that returns an
        -- empty stream.
        Left r -> pure (pure r)     -- ?!
        -- inspect first line...
        Right line -> lift (S.uncons line) >>= \case
            -- empty line:  return it and continue
            Left      s'      -> yields (pure s') >>= takeBody

            -- header line:  return empty stream that returns an empty
            -- stream that starts with the restored line.
            Right ('>',line') -> pure (wrap (S.cons '>' line'))

            -- body line:  return a stream that begins with it, recurse
            -- on its tail
            Right ( c ,line') -> yields (S.cons c line') >>= takeBody

    ignoreBody :: Monad m => Stream (S.ByteString m) m r -> Stream (S.ByteString m) m r
    ignoreBody = lift . inspect >=> \case
        -- empty input stream:  return an empty stream
        Left r -> pure r
        -- inspect first line...
        Right line -> lift (S.uncons line) >>= \case
            -- empty line:  continue
            Left      s'      -> ignoreBody s'

            -- header line:
            -- return a stream that starts with the restored line.
            Right ('>',line') -> wrap (S.cons '>' line')

            -- body line:  continue
            Right ( _ ,line') -> ignoreBody $ effect (S.effects line')

unexpected :: String -> a
unexpected msg = error $ "Ph'nglui mglw'nafh Cthulhu R'lyeh wgah'nagl fhtagn."
        ++ if null msg then "" else "\n (" ++ msg ++ ")"

readNumIO :: String -> IO Int
readNumIO s = case reads s of
    [(n,[ ])] -> return n
    [(n,"k")] -> return $ n * 1000
    [(n,"M")] -> return $ n * 1000000
    [(n,"G")] -> return $ n * 1000000000
    _         -> fail $ "unable to parse: " ++ show s


