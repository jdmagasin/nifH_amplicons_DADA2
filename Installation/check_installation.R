#!/usr/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin
##
## Check for required R packages and external tools used by the pipeline.
##
## Usage:
##   To check for just the pipeline use:
##      check_installation.R
##   To check for the pipeline as well as the ancillary scripts use:
##      check_installation.R  ancillary
##

args <- commandArgs(T)
ancillary <- (length(args) >= 1) && (args[1] == 'ancillary')

cat("\nChecking for R packages required by the DADA2 nifH pipeline:  ")
## List below based on grep'ing for libary() in all R scripts.  Some packages
## are required by dada2 (e.g. ShortRead) and not listed in
## environment_DADA2_nifH.yml
needPacks <- c('dada2','ShortRead','MASS','ggplot2','reshape2','ggrepel','vegan','digest')
installedPacks <- data.frame(installed.packages())
missing <- which(! needPacks %in% installedPacks$Package)
if (length(missing) == 0) {
    cat("Good news -- no missing R packages.\n")
} else {
    cat("\n   The following R packages are missing: ", needPacks[missing], "\n")
}


cat("\nChecking for external tools required by the pipeline:  ")
needTools <- c('cutadapt','FragGeneScan','hmmalign')
missing <- sapply(paste('which',needTools), system, ignore.stdout=T, ignore.stderr=T)
missing <- which(missing != 0)
if (length(missing) == 0) {
    cat("Good news -- no missing external tools.\n")
} else {
    cat("\n   The following tools are missing: ", needTools[missing], "\n")
}


##------------------------------------------------------------------------------
## Remainder of script checks for required R packages and external tools used by
## the ancillary scripts.
##

if (ancillary) {
    cat("\nChecking for external tools and packages required by the pipeline ancillary scripts:  ")
    needTools <- c('fastqc','vsearch','blastn','blastx','python')
    missing <- sapply(paste('which',needTools), system, ignore.stdout=T, ignore.stderr=T)
    missing <- which(missing != 0)
    if (length(missing) == 0) {
        cat("Good news -- no missing external tools.\n")
    } else {
        cat("\n   The following tools are missing: ", needTools[missing], "\n")
    }

    cat("\nChecking for Biopython:  ")
    if (system('pip list | grep biopython', ignore.stdout=T, ignore.stderr=T) == 0) {
        cat("Good, you have it.\n")
    } else {
        cat("\n   Biopython is not among the results of 'pip list'\n")
    }
}

quit(save='no')
