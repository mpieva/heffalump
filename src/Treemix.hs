module Treemix( main_treemix ) where

import Bio.Prelude
import Codec.Compression.GZip           ( compress )
import Data.ByteString.Builder          ( toLazyByteString, intDec, char7, byteString )
import System.Console.GetOpt
import System.FilePath                  ( takeBaseName )
import System.IO                        ( stdout )

import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.HashMap.Strict             as H
import qualified Data.Vector.Unboxed             as U

import Bed
import Lump
import NewRef
import Util

data ConfTmx = ConfTmx {
    conf_noutgroups :: Int,
    conf_indivs     :: Maybe FilePath,
    conf_reference  :: Maybe FilePath,
    conf_regions    :: Maybe FilePath,
    conf_transv     :: Bool,
    conf_chroms     :: (Int,Int),
    conf_output     :: L.ByteString -> IO () }

defaultConfTmx :: ConfTmx
defaultConfTmx = ConfTmx 1 Nothing Nothing Nothing False (0,maxBound) (L.hPut stdout)

opts_treemix :: [ OptDescr (ConfTmx -> IO ConfTmx) ]
opts_treemix =
    [ Option "o" ["output"]        (ReqArg set_output "FILE") "Write output to FILE (.tmx.gz)"
    , Option "r" ["reference"]     (ReqArg set_ref    "FILE") "Read reference from FILE (.2bit)"
    , Option "i" ["individuals"]   (ReqArg set_indiv  "FILE") "Read individuals from FILE (.ind)"
    , Option "n" ["numoutgroups"]  (ReqArg set_nout    "NUM") "The first NUM individuals are outgroups (1)"
    , Option "t" ["transversions"] (NoArg  set_transv       ) "Restrict to transversion variants"
    , Option "a" ["autosomes"]     (NoArg  set_autosomes    ) "Analyze only (human) autosomes"
    , Option "x" ["x-chromosome"]  (NoArg  set_xchrom       ) "Analyze only (human) X chromsome"
    , Option "y" ["y-chromosome"]  (NoArg  set_ychrom       ) "Analyze only (human) Y chromsome"
    , Option "R" ["regions"]       (ReqArg set_rgns   "FILE") "Restrict to regions in bed-file FILE" ]
  where
    set_output "-" c =                    return $ c { conf_output     = L.hPut stdout }
    set_output  a  c =                    return $ c { conf_output     = L.writeFile a }
    set_ref     a  c =                    return $ c { conf_reference  =        Just a }
    set_indiv   a  c =                    return $ c { conf_indivs     =        Just a }
    set_nout    a  c = readIO a >>= \n -> return $ c { conf_noutgroups =             n }
    set_transv     c =                    return $ c { conf_transv     =          True }
    set_autosomes  c =                    return $ c { conf_chroms     =        (0,21) }
    set_xchrom     c =                    return $ c { conf_chroms     =       (22,22) }
    set_ychrom     c =                    return $ c { conf_chroms     =       (23,23) }
    set_rgns    a  c =                    return $ c { conf_regions    =        Just a }

main_treemix :: [String] -> IO ()
main_treemix args = do
    ( hefs, ConfTmx{..} ) <- parseOpts True defaultConfTmx (mk_opts "treemix" [] opts_treemix) args

    -- We read and merge all the HEF files (shell trickery is suggested
    -- to assemble the horrible command line).  We use the IND file to
    -- map them to populations.

    (pops, npops, popixs) <- case conf_indivs of
        Just fp -> do popmap <- readIndFile <$> B.readFile fp
                      return . toSymtab $ map (lookupHef popmap) hefs
        Nothing -> return (map (B.pack . takeBaseName) hefs, length hefs, U.enumFromN 0 (length hefs))

    (ref,inps) <- decodeMany conf_reference hefs
    region_filter <- mkBedFilter conf_regions (either error nrss_chroms ref)

    conf_output $ compress $ toLazyByteString $
        foldr (\a k -> byteString a <> char7 ' ' <> k) (char7 '\n') pops <>
        foldMap
            (\Variant{..} -> let ve = U.foldl' (.|.) 0 $ U.drop conf_noutgroups v_calls
                                 is_ti = not conf_transv || isTransversion v_alt

                                 refcounts = U.accumulate (+) (U.replicate npops 0) $
                                             U.zip popixs $ U.map (fromIntegral . (3 .&.)) v_calls
                                 altcounts = U.accumulate (+) (U.replicate npops 0) $
                                             U.zip popixs $ U.map (fromIntegral . (`shiftR` 2)) v_calls

                                 show1 (a,b) k = intDec a <> char7 ',' <> intDec b <> char7 ' ' <> k

                             -- samples (not outgroups) must show ref and alt allele at least once
                             in if inRange conf_chroms v_chr && ve .&. 3 /= 0 && ve .&. 12 /= 0 && is_ti
                               then U.foldr show1 (char7 '\n') $ U.zip refcounts altcounts
                               else mempty)
            (region_filter $ concat $ mergeLumps conf_noutgroups inps)


-- | Reads an individual file.  Returns a map from individual to pop
-- population number.
readIndFile :: B.ByteString -> [(B.ByteString, B.ByteString)]
readIndFile = mapMaybe get1 . map B.words . filter (not . B.isPrefixOf "#") . B.lines
  where
    get1 (x:_:y:_) = Just (x,y)
    get1     _     = Nothing

-- | Finds the population for a file.  If there is exactly one
-- individual whose name is a prefix of the basename of the filepath,
-- the result is its associated population.
lookupHef :: [(B.ByteString, B.ByteString)] -> FilePath -> B.ByteString
lookupHef assocs fp = case matches of
    (_ ,y) : []                   -> y
    (x1,y) : (x2,_) : _ | x1 > x2 -> y
    _:_:_ -> error $ "Multiple populations match " ++ takeBaseName fp
    [   ] -> error $ "No population matches " ++ takeBaseName fp
  where
    matches = reverse. sort . map (\(x,y) -> (B.length x, y)) .
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
