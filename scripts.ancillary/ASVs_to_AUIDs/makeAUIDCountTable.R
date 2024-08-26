#!/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin

##
## Create AUID count tables by returning to the ASV count tables and
## substituting in the AUID for the ASV.  The ASV-to-AUID mapping is within the
## table assignAUIDs2ASVs.R. The table also includes a Source column for the ASV
## FASTA, which should be alongside the ASV count table.  So we should have
## enough to make the relabled tables.
##
## By hand verified for 3-4 ASVs, from different DADA2 runs, that the counts in
## the original ASV abundance table match those for the AUID in the relevant
## count.[tag].tsv.
##

usageStr <- "
Usage:
    makeAUIDcountTable.R  auidInfo.tsv [options]

Convert a table of ASV counts into AUID counts, i.e. relabel ASVs by their AUIDs
as indicated in the auidInfo.tsv table (created by assignAUIDs2ASVs.R). The ASV
count tables are assumed to be in the same directory as the ASV FASTA files,
which are in the Source column of auidInfo.  Note that only the ASVs that can be
mapped (i.e. that appear in auidInfo) will be in the output tables. Output
tables are created for each distinct ASV count table, but they are given names
prefixed with 'count' followed by the Tag column in auidInfo, since each Tags
should uniquely describe the DADA2 run that produced the ASV table,
e.g. 'Arctic2016:DNA:Filt0.2'

Options:
    -md5  Relabel ASVs with their md5 checksums, not their AUIDs.

"

args <- commandArgs(trailingOnly=T)
auidInfoTsv <- args[1]
if (length(auidInfoTsv) == 0 || !file.exists(auidInfoTsv)) {
    cat(usageStr)
    quit(save='no')
}
useMd5 <- (length(args) >1 && args[2]=='-md5')


auidTab <- read.table(auidInfoTsv, header=T, sep='\t')[,c('Source','ASVid','AUID','md5','Tag','Length')]

## Create a list of ASV count tables.  Assume that the ASV count tables are
## named identically to the FASTA except that they have suffix .tsv.
x <- unique(auidTab[c('Source','Tag')])
x$Source <- sub('\\.(fasta|fna)$','.tsv',x$Source)
asvTsvs <- x$Source
names(asvTsvs) <- x$Tag
x <- which(!file.exists(asvTsvs))
if (length(x) > 0) {
    stop("Could not find some of the ASV tables:\n", paste(asvTsvs[x], collapse='\n'))
}
if (any(table(names(asvTsvs)) > 1)) {
    stop("Some of the Tags are associated with more than one ASV table.")
}


## Helper to remap one ASV table. Uses global auidTab
MapOneAsvTable <- function(tsv, tag, useMd5=F, verbose=T)
{
    stopifnot(tag %in% auidTab$Tag)
    auidSS <- subset(auidTab, Tag==tag)
    asvTab <- read.table(tsv, sep='\t', header=T, row.names=1, check.names=F)
    idx <- match(rownames(asvTab), auidSS$ASVid)            # Which ASVs have AUIDs
    if (verbose) {
        cat("The ASV table", tsv, "has tag", tag, "and contains", nrow(asvTab), "ASVs. ",
	    sum(!is.na(idx)), "of the ASVs were previously mapped by assignAUIDs2ASVs.R to AUIDs",
	        "and will be included in the new AUID count table.\n")
    }
    asvTab <- asvTab[!is.na(idx),,drop=F]                  # Throw out unmapped ASVs
    newLab <- c('AUID','md5')[useMd5+1]
    rownames(asvTab) <- auidSS[ idx[!is.na(idx)], newLab]  # Relabel using the AUIDs
    asvTab
}

for (tag in names(asvTsvs)) {
    ## Make a valid R name, which should also be a valid file name. At least
    ## slashes and other wierd characters will be turned to periods.
    outfile <- paste0('counts.',make.names(tag),'.tsv')
    cat("### Working on",outfile,"###\n")
    if (file.exists(outfile)) { stop(outfile," already exits!") }
    newTab <- MapOneAsvTable(asvTsvs[tag], tag, useMd5)
    write.table(newTab, file=outfile, sep='\t', quote=F)
    cat("\n")
}
