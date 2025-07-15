#!//usr/local/bin/rscript

## Copyright (C) 2023 Jonathan D. Magasin

## Recurisvely search a directory for FASTQ files. Use all of them
## to create an error model using DADA2's learnErrors().  The intent
## is to use this script on FASTQ's that have been filtered to have
## just one kind of biological sequence.  For example,
##    -- Contaminants removed, somehow.  (Next step probably will
##       address anyway.)
##
##    -- Only sequences that are sufficiently like the target gene, identified
##       using blast or HMMER3 with the appropriate HMM.
##

## Load utils. run_DADA2_pipeline.sh is either the one in "bin" found by the user's
## PATH, or it is in the current directory. Either way, "scripts/" is close by.
x <- suppressWarnings(system2('which', 'run_DADA2_pipeline.sh', stdout=T))
if (length(x) == 0) { stop("Failed to find scripts/Utils/Rutils.R.") }
x <- file.path(sub('bin$','',dirname(x)), "scripts","Utils","Rutils.R")
if (!file.exists(x)) { stop("Failed to find scripts/Utils/Rutils.R.") }
source(x); rm(x)

args <- commandArgs(trailingOnly=T)
fqDir <- args[1]
if (!dir.exists(fqDir)) {
    stop("Need the directory to search for fastq.gz files.")
}
fastqList <- list.files(fqDir, 'fastq\\.gz', recursive=T, full.names=T)
##fastqList <- fastqList[1:3]  # Development: test with just 3.
if (length(fastqList) == 0) {
    stop("Found 0 files below directory '",fqDir,"' that have extension fastq.gz")
}
cat("Found",length(fastqList),"FASTQ files under",fqDir,"that will be used for learning",
    "base calling error rates.\n")

## Bail if any of the main output files exist. (Below check intermediate files.)
if (any(sapply(c('errorModel.rds','errorModel.pdf'), file.exists))) {
   stop("Looks like you already have results from this script in this directory. ",
         "Aborting!")
}

cat("Loading DADA2 and ggplot2...")
suppressMessages(library(dada2))
suppressMessages(library(ggplot2))
cat("done.\n\n")

cat("First filter and trim the reads using reasonable parameters. The goal is to\n",
    "find parameters that produce a good enough error model.  No optimization.\n")
fastqList.filtered <- gsub('\\.fastq.gz$', '.trimForErrorModel.fastq.gz', fastqList)
idx <- which(sapply(fastqList.filtered, file.exists))
if (length(idx) > 0) {
    stop("Aborting because these files already exist:\n",
         paste(fastqList.filtered[idx], collapse=', '), "\n")
}
track <- Graceful_filterAndTrim(fwd=fastqList, filt=fastqList.filtered,
                                truncQ=10, maxEE=2, maxN=0, minLen=80)
if (is.null(track)) {
    stop("Aborting because filterAndTrim() failed.")
}
df <- data.frame(FASTQ        = sapply(fastqList, dirname),
                 reads.in     = track[,'reads.in'], 
                 reads.out    = track[,'reads.out'],
                 PctReadsKept = round(100*track[,'reads.out'] / track[,'reads.in'], 1))
rownames(df) <- NULL
print(df)
rm(df)
cat("\n\n")


cat("Now the main event, learning the error rates that best predict the data.",
    "If there are lots of data sets, then you have time get some coffee and maybe",
    "even to roast the beans.\n")

## Mainly default params (except randomize and multithread, and use 50% more for
## nbases) but let's be clear.  Note that learnErrors() will stop reading in
## data sets once it has seen nbases.  My pipeline has a nifty way to force
## dada2 to use all the data.  But probably I won't need this.
errors <- learnErrors(fastqList.filtered, nbases = 1.5e+08,
                      errorEstimationFunction = loessErrfun,
                      randomize=TRUE, multithread=TRUE, verbose=TRUE)
cat("Done learning errors!\n\n")

cat("Saving error model to errorModel.rds\n")
saveRDS(errors, "errorModel.rds")

cat("Plotting error model and saving to errorModel.pdf.\n")
plotErrors(errors, nominalQ=T)
ggsave('errorModel.pdf')

cat("All done!\n")
quit(save='no')
