module Eigenstrat where

import BasePrelude

import qualified Data.Foldable                  as F
import qualified Data.ByteString.Char8          as B
import qualified Data.ByteString.Lazy           as LB
import qualified Data.ByteString.Lazy.Char8     as L

-- ^ The Eigenstrat and Ancestrymap exporters.  Includes silly
-- generation of names and hash functions.



-- We wrote Eigenstrat in the first iteration.  Eigenstrat is one
-- character per genotype.  Alternatively, we could write packed
-- Ancestrymap, then it's one char per four genotypes.  It starts with a
-- header record, which contains text and is then padded to the record
-- size as defined below:
-- printf("GENO %d %d %x %x", n_individuals, n_snps, individual_hash, snp_hash);
--
-- That's of course retarded, we can't know n_snps or snp_hash before we
-- generate the snp file.  So all we can do is reserve space and print
-- the header later.  m(
--
-- To generate individual_hash, we apply nick_hasharr to the vector of
-- individual names; to generate snp_hash we apply nick_hasharr to the
-- vector of SNP identifiers.
--
-- The remainder is binary:  one record per snp, each record uses 2 bits
-- per individual, padded to a full byte, padded to at least 48 bytes.


-- Nickhash for strings.  XXX It's not clear if we got the signs right.
nick_hashit :: L.ByteString -> Int32
nick_hashit = LB.foldl (\h c -> 23 * h + fromIntegral c) 0

-- Nickhash for vectors of strings.  (This needs incremental updates,
-- but that's easy.
nick_hasharr :: F.Foldable v => v L.ByteString -> Int32
nick_hasharr = F.foldl (\h s -> 17 * h `xor` nick_hashit s) 0

mkname :: Int -> Int -> Char -> String
mkname x y z = enc (4*y + numOf z) [ B.index chars x ]
  where
    numOf 'A' = 0 ; numOf 'C' = 1 ; numOf 'G' = 2 ; numOf 'T' = 3
    numOf 'a' = 0 ; numOf 'c' = 1 ; numOf 'g' = 2 ; numOf  _  = 3

    chars = "0123456789ABCDEFGHKLMNPQRSTUWXYZ"

    enc 0 = id
    enc u = (:) (B.index chars (u .&. 31)) . enc (u `shiftR` 5)

