module Treemix( main_treemix ) where

import BasePrelude
import Data.ByteString.Builder          ( hPutBuilder, intDec, char7 )
import System.Console.GetOpt
import System.FilePath                  ( takeBaseName )
import System.IO

import qualified Data.ByteString.Char8           as B
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.HashMap.Strict             as H
import qualified Data.Vector.Unboxed             as U

import Merge
import Stretch
import Util

data ConfTmx = ConfTmx {
    conf_noutgroups :: Int,
    conf_indivs     :: Maybe FilePath,
    conf_output     :: FilePath,
    conf_reference  :: FilePath }
  deriving Show

defaultConfTmx :: ConfTmx
defaultConfTmx = ConfTmx 1 Nothing
                         (error "no output file specified")
                         (error "no reference file specified")

opts_treemix :: [ OptDescr (ConfTmx -> IO ConfTmx) ]
opts_treemix =
    [ Option "o" ["output"] (ReqArg set_output "FILE") "Write output to FILE (.tmx)"
    , Option "r" ["reference"] (ReqArg set_ref "FILE") "Read reference from FILE (.fa)"
    , Option "i" ["individuals"] (ReqArg set_indiv "FILE") "Read individuals from FILE (.ind)"
    , Option "n" ["numoutgroups"] (ReqArg set_nout "NUM") "The first NUM individuals are outgroups (1)" ]
  where
    set_output a c =                    return $ c { conf_output     =      a }
    set_ref    a c =                    return $ c { conf_reference  =      a }
    set_indiv  a c =                    return $ c { conf_indivs     = Just a }
    set_nout   a c = readIO a >>= \n -> return $ c { conf_noutgroups =      n }

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

    ref <- readReference conf_reference
    inps <- mapM (fmap (decode . decomp) . L.readFile) hefs

    withFile conf_output WriteMode $ \hout -> do
        B.hPut hout $ B.unwords pops `B.snoc` '\n'
        forM_ (merge_hefs False conf_noutgroups ref inps) $ \Variant{..} -> do
            -- XXX  Some variants are probably useless.  Should we
            -- attempt to remove private variants?  Variants where
            -- everyone is different from the reference?

            let refcounts = U.accumulate (+) (U.replicate npops 0) $
                            U.zip popixs $ U.map (fromIntegral . (3 .&.)) v_calls
                altcounts = U.accumulate (+) (U.replicate npops 0) $
                            U.zip popixs $ U.map (fromIntegral . (`shiftR` 2)) v_calls

            let show1 (a,b) k = intDec a <> char7 ',' <> intDec b <> char7 ' ' <> k
            hPutBuilder hout $ U.foldr show1 (char7 '\n') $ U.zip refcounts altcounts


-- | Reads an individual file.  Returns a map from individual to pop
-- population number.
-- into the population list, and the population list.
readIndFile :: B.ByteString -> [(B.ByteString, B.ByteString)]
readIndFile = mapMaybe get1 . map B.words . filter (not . B.isPrefixOf "#") . B.lines
  where
    get1 (x:_:y:_) = Just (x,y)
    get1     _     = Nothing

-- | Finds the population for a file.  If there is exactly one
-- individual whose name is a prefix of the basename of the filepath,
-- the result is its associated population.
lookupHef :: [(B.ByteString, B.ByteString)] -> FilePath -> B.ByteString
lookupHef assocs fp = case filter (\(x,_) -> x `B.isPrefixOf` B.pack (takeBaseName fp)) assocs of
    [ y ] -> snd y
    _:_:_ -> error $ "Multiple populations match " ++ takeBaseName fp
    [   ] -> error $ "No population matches " ++ takeBaseName fp

-- | Takes a list of stuff, returns the list without dups, its length,
-- and a list of symbols.
toSymtab :: [B.ByteString] -> ([B.ByteString], Int, U.Vector Int)
toSymtab = go [] [] H.empty 0
  where
    go ss is _ n [    ] = (reverse ss, n, U.reverse $ U.fromList is)
    go ss is h n (x:xs) = case H.lookup x h of
            Just  i -> go    ss  (i:is)               h        n  xs
            Nothing -> go (x:ss) (n:is) (H.insert x n h) (succ n) xs

