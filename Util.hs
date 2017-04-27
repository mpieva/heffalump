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
decomp s0 = case B.uncons s0 of
    Just (31, s') -> case B.uncons s' of
        Just (139,_) -> decomp_loop s0
        _            -> s0
    _                -> s0
  where
    decomp_loop :: L.ByteString -> L.ByteString
    decomp_loop s = case B.uncons s of
        Just (31, s') -> case B.uncons s' of
            Just (139,_) -> Z.foldDecompressStreamWithInput L.Chunk decomp_loop throw
                            (Z.decompressST Z.gzipOrZlibFormat Z.defaultDecompressParams) s
            _            -> L.empty  -- ignores trailing garbage
        _                -> L.empty  -- ignores trailing garbage

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


-- Samples in FastA format are treated as diploid.
{-# DEPRECATED readSampleFa "No good." #-}
readSampleFa :: FilePath -> IO [( S.ByteString, L.ByteString )]
readSampleFa fp = parseFasta . L.lines . decomp <$> L.readFile fp

-- | Parsing a FastA file results in pairs of sequence names and
-- sequences.  The sequences are still text with their case preserved.
{-# DEPRECATED parseFasta "No good." #-}
parseFasta :: [L.ByteString] -> [(S.ByteString, L.ByteString)]
parseFasta  = go . dropWhile (not . L.isPrefixOf ">")
  where
    go [    ] = [ ]
    go (h:ls) = ( name, L.concat body ) : parseFasta rest
      where
        (body,rest) = break (L.isPrefixOf ">") ls
        name = S.concat . concatMap L.toChunks . take 1 . L.words $ L.drop 1 h


-- Cheap 'toLower' function, best applied to FastA sequences only.
{-# INLINE low #-}
{-# DEPRECATED low "Do we need this?" #-}
low :: Word8 -> Word8
low !x = x .|. 32

-- Cheap 'toUpper' function, best applied to FastA sequences only.
{-# INLINE up #-}
{-# DEPRECATED up "Do we need this?" #-}
up :: Word8 -> Word8
up !x = x .&. complement 32

-- | Our expected chromosomes.
{-# DEPRECATED chroms "Use proper reference." #-}
chroms :: [L.ByteString]
chroms = L.words "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

readNumIO :: String -> IO Int
readNumIO s = case reads s of
    [(n,[ ])] -> return n
    [(n,"k")] -> return $ n * 1000
    [(n,"M")] -> return $ n * 1000000
    [(n,"G")] -> return $ n * 1000000000
    _         -> fail $ "unable to parse: " ++ show s


