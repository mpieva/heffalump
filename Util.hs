module Util where

import BasePrelude
import System.Console.GetOpt
import System.IO                        ( hPutStrLn, stderr )

import qualified Data.ByteString                 as S
import qualified Data.ByteString.Lazy            as B
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.ByteString.Lazy.Internal   as L ( ByteString(..) )
import qualified Codec.Compression.Zlib.Internal as Z

decomp :: L.ByteString -> L.ByteString
decomp s = case B.uncons s of
    Just (31, s') -> case B.uncons s' of
        Just (139,_) -> Z.foldDecompressStreamWithInput L.Chunk decomp throw
                        (Z.decompressST Z.gzipOrZlibFormat Z.defaultDecompressParams) s
        _            -> s
    _                -> s

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


newtype Reference = Reference [L.ByteString]

-- Reference is FastA format, treated as haploid.
readReference :: FilePath -> IO ([S.ByteString], Reference)
readReference fp = (id *** Reference) . unzip . parseFasta chroms . L.lines . decomp <$> L.readFile fp

-- Samples in FastA format are treated as diploid.
readSampleFa :: FilePath -> IO [L.ByteString]
readSampleFa fp = map snd . parseFasta chroms . L.lines . decomp <$> L.readFile fp

parseFasta :: [L.ByteString] -> [L.ByteString] -> [(S.ByteString, L.ByteString)]
parseFasta [    ]  _ = []
parseFasta (c:cs) ls =
    case dropWhile (not . isHeader) ls of
        [   ] -> error $ "expected chromosome " ++ show c ++ " not found"
        h:ls' -> case span isBody ls' of
            (mine,rest) -> ( S.concat . concatMap L.toChunks . take 1 . L.words $ L.tail h, L.concat mine )
                           : parseFasta cs rest
  where
    isBody s = L.null s || L.head s /= '>'
    isHeader s = case L.words s of
        nm:_ -> nm == '>' `L.cons` c
        [  ] -> False

-- Cheap 'toLower' function, best applied to FastA sequences only.
{-# INLINE low #-}
low :: Word8 -> Word8
low !x = x .|. 32

-- Cheap 'toUpper' function, best applied to FastA sequences only.
{-# INLINE up #-}
up :: Word8 -> Word8
up !x = x .&. complement 32

-- | Our expected chromosomes.
chroms :: [L.ByteString]
chroms = L.words "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
-- chroms = L.words "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"

readNumIO :: String -> IO Int
readNumIO s = case reads s of
    [(n,[ ])] -> return n
    [(n,"k")] -> return $ n * 1000
    [(n,"M")] -> return $ n * 1000000
    [(n,"G")] -> return $ n * 1000000000
    _         -> fail $ "unable to parse: " ++ show s


