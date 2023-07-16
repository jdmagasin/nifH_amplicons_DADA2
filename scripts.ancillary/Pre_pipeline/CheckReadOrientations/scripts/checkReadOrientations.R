#!/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin

usage <- '
Usage:

  reorientReadsDir.R  fwdPrimer  revPrimer  fastqDir

Older Illumina paired read libraries were sometimes prepared such that amplicons
appeared in both orientations between the adapters. Consequently the forward
("R1") FASTQs had a mix of the target amplicon and its reverse complement, and
similar for the reverse ("R2") FASTQs. This is a problem for tools that expect
just one orientation in a FASTQ (e.g. cutadapt and DADA2).  This script
identifies reads that are in the unexpected orientation because they lack the
expected primer at the 5\' end.

The script searches recursively within fastqDir for FASTQs. Files with
"_R1.fastq" or "_1.fastq" in their name (optionally with a final ".gz") are
searched for reads that lack the fwdPrimer but have the revPrimer (in the 5\' to
3\' orientation for both primers). "_R2" files are searched for reads that lack
the revPrimer but have the fwdPrimer (again, 5\' to 3\').

Output:
  * Files that list the IDs of reads in the unexpected orientation, with suffix
    ".misoriented.ids"

'

args <- commandArgs(T)
if (length(args)==0) { cat(usage); stop("Need parameters.") }

primers <- c(fwd=args[1], rev=args[2])
stopifnot(is.character(primers))

fastqDir <- args[3]
stopifnot(dir.exists(fastqDir))

cat("Loading libraries...\n")
suppressMessages(library(ShortRead))


CheckReadOrientations <- function(fastqFile, fwdPrimer, revPrimer, minCheckLen = 50)
{
    cat("Reading",fastqFile,"...")
    fastq <- readFastq(fastqFile)
    cat("Got",length(fastq),"reads.")

    ## Truncate reads b/c we only want to look for primers towards the 5'. (Also
    ## faster.)
    w <- width(fastq)
    w[w > minCheckLen] <- minCheckLen
    firstNt <- narrow(sread(fastq), start=1, end=w)

    ## Search the first minCheckLen bp for the forward primer. Allow IUPAC
    ## wildcards to match appropriate bases. Otherwise use defaults: no
    ## mismatches, no indels.
    haveFwdPrimer <- which(vcountPattern(fwdPrimer, firstNt, fixed=F) > 0)
    if (length(haveFwdPrimer) == 0) {
        stop("Could not find the forward primer in the first ",minCheckLen,"nt of any read. Aborting.")
    }

    ## Check for the reverse primer (in 5' to 3' orientation) in the reads that
    ## lack a forward primer. Still searching just the first nt.
    haveRevPrimer <- which(vcountPattern(revPrimer, firstNt[-haveFwdPrimer], fixed=F) > 0)
    ids <- ""
    if (length(haveRevPrimer)>0) {
        ## Convert to index into fastq
        haveRevPrimer <- ((1:length(firstNt))[-haveFwdPrimer])[haveRevPrimer]
        ## First version of this script reverse-complemented the backwards 
        ## reads and wrote out a new FASTQ. That's not helpful b/c it puts low
        ## quality base calls at the 5' which makes trimming (e.g. in DADA2)
        ## unworkable.
        ##fastq[haveRevPrimer] <- reverseComplement(fastq[haveRevPrimer])
        ids <- as.character(id(fastq)[haveRevPrimer])
    }
    cat("\n")
    ##return( list(fastq=fastq, ids=ids) )
    return( ids )
}


fastqList <- list.files(fastqDir, pattern="*_R{0,1}[1,2]\\.fastq(.gz|)*", recursive=T, full.names=T)
for (fq in fastqList) {
    fwdrev <- primers
    if (grepl('_R{0,1}2\\.fastq',fq)) { fwdrev <- rev(primers) }
    ids <- CheckReadOrientations(fq, fwdrev[1], fwdrev[2])
    outfile <- file.path(dirname(fq), paste0(sub('\\.fastq.*','',basename(fq)), '.misoriented.ids'))
    writeLines(ids, outfile)
}
cat("Done!\n")
