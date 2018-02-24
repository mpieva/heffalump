module Treemix( main_treemix ) where

import Bio.Prelude
import Control.Monad.Trans.Writer       ( Writer, execWriter, tell )
import Data.ByteString.Builder          ( Builder, intDec, char7, byteString )
import Data.ByteString.Streaming        ( concatBuilders, toHandle, toStreamingByteString )
import Streaming
import System.Console.GetOpt
import System.FilePath                  ( takeBaseName )
import System.IO                        ( stdout, withFile, IOMode(..) )
import Text.Regex.Posix

import qualified Data.ByteString.Char8           as B
import qualified Data.HashMap.Strict             as H
import qualified Data.Vector.Storable            as U
import qualified Streaming.Prelude               as Q

import Bed
import Genome
import Lump
import Util

data ConfTmx = ConfTmx {
    conf_noutgroups :: Int,
    conf_indivs     :: Maybe FilePath,
    conf_reference  :: Maybe FilePath,
    conf_regions    :: Maybe FilePath,
    conf_transv     :: Bool,
    conf_chroms     :: Maybe Regex,
    conf_output     :: (Handle -> IO ()) -> IO () }

defaultConfTmx :: ConfTmx
defaultConfTmx = ConfTmx 0 Nothing Nothing Nothing False Nothing ($ stdout)

opts_treemix :: [ OptDescr (ConfTmx -> IO ConfTmx) ]
opts_treemix =
    [ Option "o" ["output"]        (ReqArg set_output "FILE") "Write output to FILE (.tmx.gz)"
    , Option "r" ["reference"]     (ReqArg set_ref    "FILE") "Read reference from FILE (.2bit)"
    , Option "i" ["individuals"]   (ReqArg set_indiv  "FILE") "Read individuals from FILE (.ind)"
    , Option "n" ["numoutgroups"]  (ReqArg set_nout    "NUM") "The first NUM individuals are outgroups (0)"
    , Option "t" ["transversions"] (NoArg  set_transv       ) "Restrict to transversion variants"
    , Option "c" ["chromosomes"]   (ReqArg set_chrs "REGEX" ) "Analyze subset of chromosomes"
    , Option "a" ["autosomes"]     (NoArg  set_autosomes    ) "Analyze only autosomes"
    , Option "x" ["x-chromosome"]  (NoArg  set_xchrom       ) "Analyze only X chromsome"
    , Option "y" ["y-chromosome"]  (NoArg  set_ychrom       ) "Analyze only Y chromsome"
    , Option "R" ["regions"]       (ReqArg set_rgns   "FILE") "Restrict to regions in bed-file FILE" ]
  where
    set_output "-" c =                    return $ c { conf_output     =                ($ stdout) }
    set_output  a  c =                    return $ c { conf_output     =      withFile a WriteMode }
    set_ref     a  c =                    return $ c { conf_reference  =                    Just a }
    set_indiv   a  c =                    return $ c { conf_indivs     =                    Just a }
    set_nout    a  c = readIO a >>= \n -> return $ c { conf_noutgroups =                         n }
    set_transv     c =                    return $ c { conf_transv     =                      True }
    set_chrs    a  c =                    return $ c { conf_chroms     =                      re a }
    set_autosomes  c =                    return $ c { conf_chroms     = re "^(chr)?[0-9]+[a-z]?$" }
    set_xchrom     c =                    return $ c { conf_chroms     =            re "^(chr)?X$" }
    set_ychrom     c =                    return $ c { conf_chroms     =            re "^(chr)?Y$" }
    set_rgns    a  c =                    return $ c { conf_regions    =                    Just a }

    re :: String -> Maybe Regex
    re = Just . makeRegex

main_treemix :: [String] -> IO ()
main_treemix args = do
    ( hefs, ConfTmx{..} ) <- parseFile1Opts defaultConfTmx "treemix" [] opts_treemix args

    -- We read and merge all the HEF files (shell trickery is suggested
    -- to assemble the horrible command line).  We use the optional IND
    -- file to map them to populations.

    (pops, npops, popixs) <- case conf_indivs of
        Just fp -> do popmap <- readIndFile <$> B.readFile fp
                      return . toSymtab $ map (lookupHef popmap) hefs
        Nothing -> return (map (B.pack . takeBaseName) hefs, length hefs, U.enumFromN 0 (length hefs))

    decodeMany conf_reference hefs $ \ref inps -> do
      region_filter <- mkBedFilter conf_regions (either error rss_chroms ref)
      conf_output $ \hdl ->
        toHandle hdl $ gzip $ toStreamingByteString $

        (<>) (foldr (\a k -> byteString a <> char7 ' ' <> k) (char7 '\n') pops) $
        concatBuilders $ Q.map
            (\Variant{..} -> let AC nref nalt = U.foldl' (<>) mempty $ U.drop conf_noutgroups v_calls
                                 is_ti = not conf_transv || isTransversion v_alt

                                 refcounts = U.accumulate_ (+) (U.replicate npops 0) popixs $
                                             U.map (fromIntegral . ac_num_ref) v_calls
                                 altcounts = U.accumulate_ (+) (U.replicate npops 0) popixs $
                                             U.map (fromIntegral . ac_num_alt) v_calls

                                 show1 :: Int -> Int -> Writer Builder ()
                                 show1 a b = tell $ intDec a <> char7 ',' <> intDec b <> char7 ' '

                             -- samples (not outgroups) must show ref and alt allele at least once
                             in if nref /= 0 && nalt /= 0 && is_ti
                               then execWriter (U.zipWithM_ show1 refcounts altcounts) <> char7 '\n'
                               else mempty) $
        region_filter $ chrom_filter (either error rss_chroms ref) conf_chroms $
        Q.concat $ mergeLumps conf_noutgroups inps

chrom_filter :: Monad m => [B.ByteString] -> Maybe Regex -> Stream (Of Variant) m r -> Stream (Of Variant) m r
chrom_filter _chroms  Nothing  = id
chrom_filter  chroms (Just re) = go [ i | (i,chrom) <- zip [0..] chroms, matchTest re chrom ]
  where
    go [    ] = lift . Q.effects
    go (c:cs) = go1
      where
        go1 = lift . Q.next >=> \case
            Left r -> pure r
            Right (v,vs)
                | c < v_chr v -> go cs (Q.cons v vs)
                | c > v_chr v -> go1 vs
                | otherwise   -> v `Q.cons` go1 vs

-- | Reads an individual file.  Returns a map from individual to pop
-- population number.
readIndFile :: B.ByteString -> [(B.ByteString, B.ByteString)]
readIndFile = mapMaybe (get1 . B.words) . filter (not . B.isPrefixOf "#") . B.lines
  where
    get1 (x:_:y:_) = Just (x,y)
    get1     _     = Nothing

-- | Finds the population for a file.  If there is exactly one
-- individual whose name is a prefix of the basename of the filepath,
-- the result is its associated population.
lookupHef :: [(B.ByteString, B.ByteString)] -> FilePath -> B.ByteString
lookupHef assocs fp = case matches of
    (_ ,y) : [    ]               -> y
    (x1,y) : (x2,_) : _ | x1 > x2 -> y
    _:_:_ -> error $ "Multiple populations match " ++ takeBaseName fp
    [   ] -> error $ "No population matches " ++ takeBaseName fp
  where
    matches = sortBy (flip compare) . map (first B.length) .
              filter (\(x,_) -> x `B.isPrefixOf` B.pack (takeBaseName fp)) $ assocs

-- | Takes a list of stuff, returns the list without dups, its length,
-- and a list of symbols.
toSymtab :: [B.ByteString] -> ([B.ByteString], Int, U.Vector Int)
toSymtab = go [] [] H.empty 0
  where
    go ss is _ n [    ] = (reverse ss, n, U.reverse $ U.fromList is)
    go ss is h n (x:xs) = case H.lookup x h of
            Just  i -> go    ss  (i:is)               h        n  xs
            Nothing -> go (x:ss) (n:is) (H.insert x n h) (succ n) xs

