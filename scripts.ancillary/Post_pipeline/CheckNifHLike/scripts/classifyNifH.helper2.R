#!/usr/bin/env Rscript

## Copyright (C) 2024 Jonathan D. Magasin

## Usage:
##   classifyNifH.helper2.R  queryFasta  positivesFasta
##
## Identify sequences in queryFasta that exactly match sequences in positivesFasta.
##

args <- commandArgs(trailingOnly=T)
queryFasta <- args[1]
positivesFasta <- args[2]
cat("Will use positive nifH sequences in", positivesFasta, "to quickly identify known nifH",
    "in", queryFasta, ".\n")
stopifnot(file.exists(queryFasta))
stopifnot(file.exists(positivesFasta))

suppressMessages(library(ShortRead))
qSeqs <- readFasta(queryFasta)
stopifnot(length(qSeqs) > 0)
pSeqs <- readFasta(positivesFasta)
stopifnot(length(pSeqs) > 0)

## Determine which queries have exact matches.  Careful since match() only finds
## 1st match and queries could be non-unique (in theory -- won't happen in the
## workflow).  Code below works if every (or no) qSeqs match pSeqs.
idx <- match(sread(qSeqs), sread(pSeqs))
idx.match   <- which(!is.na(idx))
idx.nomatch <- which(is.na(idx))
cat(length(idx.match), "ASVs have exact matches within", length(pSeqs), "known nifH sequences.\n")
cat(length(idx.nomatch), "ASVs do not and will have to be checked against",
    "results from ARBitrator.\n")
ids <- sub(' .*$', '', as.character(id(qSeqs)[idx.match]))
writeLines(ids, "positives_dontBlastx.ids")
ids <- sub(' .*$', '', as.character(id(qSeqs)[idx.nomatch]))
writeLines(ids, "unsure_needBlastx.ids")
