# Heffalump #

## What is Heffalump? ##

It's both a file format that stores genetic variation and the tool to work with it.  Compared to VCF, it is much more compact and faster to work with.  It also remembers which of the sites that didn't show a variant were called at all, allowing for a clean workflow.

Heffalump imports common file formats (VCF, BCF, Bam by picking from an arbitrary read at each site, HETFA (the files used by SGDP)), storing one file per sample.  Conversely, it merges multiple heffalump files and exports them into common formats (VCF, Eigenstrat) or computes simple statistics on it.

## How do I get set up? ##

* install GHC (see http://haskell.org/ghc) or have your admin install it,
* run `cabal update` (takes a while to download the current package list),
* run `cabal install
  https://bitbucket.org/ustenzel/heffalump/get/0.1.tar.gz`

## Homepage ##

https://bioinf.eva.mpg.de/heffalump/