#!/usr/local/bin/Rscript

suppressMessages(library(ShortRead))
args <- commandArgs(trailingOnly=T)
idsFile   <- args[1]
fastqFile <- args[2]

fastq <- readFastq(fastqFile)
ids <- read.table(idsFile)[,1]

## Snip fastq read id's at whitespace. Garr.  This script is ~general, but also
## must snip so I can use the FragGeneScan results. (Unless I grep...)
## Note use of as.character() to extract the ID's from the BStringIO object.
ids.snip <- gsub(' +.*$', '', as.character(id(fastq)))

## Export a fastq of just the target reads.
idx <- which(ids.snip %in% ids)
stopifnot(length(idx) == length(ids))

writeFastq(fastq[idx], 'readsExtracted.fastq.gz')
