#!/usr/bin/env Rscript


################################################################################
##
## Adaptation of the DADA2 script I wrote for Ana's data (PacBio reads for Ana's
## 18S ITS regions) to nifH amplicons.  The original script (and this one) was
## adapted from the markdown from Callahan et al. 2018 (a **PacBio** example):
##   https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html
##
## Some of the big data tips for speeding up processing have been adopted from:
##    https://benjjneb.github.io/dada2/bigdata.html
##
## Usage
##    dada2_nifH_amplicons.R [dataDir] [outputDir] [params] [precalcErrorModelsDir]
##
## where dataDir defaults to Data and outputDir defaults to Dada2Out.  Note that
## non-biological sequences should already have been removed: barcodes,
## adapters, and primers. (The primer is synthetic and does not capture
## biological variability of nifH.)  See the script runCutadapt.sh which trims
## nifH primers (and thus adapters and barcodes) from paired-end Illumina reads.
##
## This script works with DADA2 v1.16.0 and R v4.0.3 (at least). I recommend
## that you install both in a conda environment as described in the INSTALL.txt
## document distributed with the DADA2 nifH pipeline.
##
## History
## -------
##  2021 Mar 12      Initial version.
##  2021 Apr  5      Optionally use a pre-specified error model.
##  2021 Jun         Take a parameters file
##  2021 Sep  2      Allow params file, error model, or both.
##  2021 Oct  7      NMDS now with metaMDS and only color if <=15 samples.
##                   More flexible fastq name expectations (Fastq2Samp).
##  2021 Oct 26      Support ASVs that are based only on the R1 reads.  Use
##                   if the R2s are poor quality (not uncommon) & won't merge.
##  2021 Oct 30      Support truncLen parameter for filterAndTrim()
##
################################################################################

################################################################################
##
## Some of the parameters that we vary
##
NBASES_FOR_LEARNING_ERRORS = 1e+08  # This is the default number of bases that
                                    # learnErrors() wants to examine for making
                                    # error models.  This script distributes
                                    # these bases among randomly picked reads
                                    # from each sample.

## Parameters for filtering and trimming. (Many we leave as default.)  One value
## means use that value for forward and reverse strands; two values are for
## forward and reverse (resp).
## truncQ   = Truncate read at 1st position with Q <= truncQ
##
## truncLen = Truncate reads at a specified bp which you should select based on
##            quality profiles such as from FastQC. 0 by default means do not
##            truncate. This parameter treats all reads the same (ignoring
##            quality for the specific read).  It should make reads tend toward
##            the same length (if truncQ is low)) and thus could be set so that
##            dada2 sees more examples of basecalls at later positions in the
##            reads which might improve ASV inference (vs. more variable read
##            lengths if one chops each based on truncQ).  I suggest running the
##            pipeline using different approaches and checking which ASVs
##            persist.
## maxEE    = Max allowed (expected) errors in final read
##
## minLen   = Final read must be at least this long. Empirically 162 produces
##            ASVs >= 300nt
##
filterAndTrimParams <- list(truncQ = 2, maxEE.fwd = 3, maxEE.rev = 5, minLen = 20,
                            truncLen.fwd = 0, truncLen.rev = 0)

## These are the default values.  Allow overriding in the params file.
mergePairsParams <- list(minOverlap=12, maxMismatch=0, justConcatenate=F)

## Special parameters:
specialParams <- list(useOnlyR1Reads = FALSE)   # Base ASVs only on the R1 reads.
                                                # Ignore R2's at filterAndTrim and
                                                # do not mergePairs.

################################################################################
##
## Set up
##
cat("Starting on", date(),"\n")

## Assume have this many columns for tables printed to stdout
options(width = 110)

## Params:
##  - dataDir will be searched recursively for fastq's
##  - dataDir and data2Outdir must be specified if other params are.
##  - params or error models or both can be present.
##
args <- commandArgs(trailingOnly=T)
stopifnot(length(args) <= 4)
dataDir     = ifelse(!is.na(args[1]), args[1], 'Data')
dada2OutDir = ifelse(!is.na(args[2]), args[2], 'Dada2OutDir')
paramsFile  = args[3] # NA if unspecified
preCalculatedErrorModel = args[4]
if (is.na(preCalculatedErrorModel) && grepl('\\.rds$',paramsFile,ignore.case=T)) {
    ## Fixup: There was no params file, just an error models file.
    preCalculatedErrorModel <- args[3]
    paramsFile <- NA
}

if (!is.na(paramsFile)) {
    ## Parameters file. Tabular. Very simple for now.
    stopifnot(file.exists(paramsFile))
    ptab <-read.csv(paramsFile, header=F, row.names=1, comment.char = "#")
    p <- union(names(filterAndTrimParams), names(mergePairsParams))
    p <- union(p, names(specialParams))
    plist <- intersect(p, rownames(ptab))
    cat("Will use parameters file",paramsFile,"for",paste(plist,collapse=','),"\n")
    for (p in plist) {
        if (p %in% names(filterAndTrimParams)) {
            filterAndTrimParams[[p]] <- as.integer(ptab[p,1])
        }
        if (p == 'justConcatenate') {
            mergePairsParams[[p]] <- ifelse(ptab[p,1]=='TRUE',TRUE,FALSE)
        } else if (p %in% names(mergePairsParams)) {
            mergePairsParams[[p]] <- as.integer(ptab[p,1])
        }
        if (p == 'useOnlyR1Reads') {
            ## Impacts filterAndTrim and skips mergePairs.
            cat("NOTE!  You have asked for ASVs to be based just on the R1 reads.\n")
            specialParams$useOnlyR1Reads <- ifelse(ptab[p,1]=='TRUE',TRUE,FALSE)
        }
    }
    rm(ptab,plist,p)
}
stopifnot(is.na(preCalculatedErrorModel) | file.exists(preCalculatedErrorModel))

################################################################################

cat("Fastq files expected in ", dataDir, "\n")
cat("Output directory will be ", dada2OutDir, "\n")


## Extract the sample identity that is encoded at the start of a fastq file name.
## Presumes names are structured like this: {Samp}{_stuff1_}R{1,2}{stuff2}.fastq.gz
## where:
##   - Samp can have any character other than "_".
##   - stuff1, if present, can have any character but must be flanked by "_"
##     Usually stuff1 will be/include the sequencing lane (L001).  stuff1 is
##     not interesting to the study and is dropped.
##   - stuff2 can be anything, or absent
##   - The .fastq.gz could be absent, but that would be bad style.
Fastq2Samp <- function(fq) {
    fq <- basename(fq)
    x <- gsub('\\.fastq\\.gz',    '',    fq)        # drop fastq.gz if any
    gsub('^([^_]+).*_R[0-9]+.*$', '\\1', x)         # keep just {Samp}
}
Fastq2Index <- function(fq) {
    fq <- basename(fq)
    x <- gsub('\\.fastq\\.gz',          '',    fq)  # drop fastq.gz if any
    as.numeric(gsub('^.+_R([0-9]+).*$', '\\1', x))  # keep just the direction index
}


## Make sure we can find dada2, and others.  In a cluster environment (like hummingbird)
## the script might need some help finding locally installed packages.
needPkgs <- c('dada2','ShortRead',
              'MASS',  # for NMDS
              'ggplot2','reshape2','ggrepel',  # plots
              'vegan')
havePkgs <- list.files(.libPaths())
havePkgs <- sapply(needPkgs, function (p) any(grepl(p,havePkgs)))
if (!all(havePkgs)) {
    ## Retry with local path that works for me **on hummingbird.**
    .libPaths("~/R/x86_64-pc-linux-gnu-library/3.6")
    cat("Added ~/R/x86_64-pc-linux-gnu-library/3.6 to .libPaths() to find missing packages.\n")
    havePkgs <- list.files(.libPaths())
    havePkgs <- sapply(needPkgs, function (p) any(grepl(p,havePkgs)))
    if (!all(havePkgs)) {
        stop("Aborting because missing packages:  ",
             paste(names(which(havePkgs==F)), collapse=' '),
             "\nPerhaps they are not installed. They do not appear in ",
             ".libPaths()")
    }
}
x <- sapply(needPkgs, function (p) { library(p, character.only=T) })
rm(needPkgs, havePkgs, x)

cat("Will use the following dada2 version:\n")
packageVersion("dada2")


## Directories and files.  As in DADA2 examples and from the 2018 Callahan
## paper, store successively processed fastq files in deeper directories (but
## retain the original fastq names).
if (!dir.exists(dada2OutDir)) {
     dir.create(dada2OutDir)
} else {
    cat("The output directory", dada2OutDir, "already exists.  Previous results",
        "will be used when possible.\n")
}
plotsDir = file.path(dada2OutDir,'Plots');    if (!dir.exists(plotsDir)) { dir.create(plotsDir) }
rObjsDir = file.path(dada2OutDir,'RObjects'); if (!dir.exists(rObjsDir)) { dir.create(rObjsDir) }
txtsDir  = file.path(dada2OutDir,'TextData'); if (!dir.exists(txtsDir))  { dir.create(txtsDir) }


## The trimmed raw reads.  Assumptions:
##   1. Paired fastq have reads in corresponding order.  NOT CHECKED BY THIS CODE.
##   2. Naming convention of dirs and samples is such that sorting the fastq file
##      paths will make 'fastqs' entries i and i+1 be the pair for one sample.
fastqs <- sort(list.files(dataDir, pattern="fastq.gz", full.names=TRUE, recursive=TRUE))
CheckFastq <- function(fastqs)
{
    ## Non-portable check that fastq's i and i+1 are for the same sample
    cat("Checking that have paired FASTQ for each sample.\n")
    stopifnot(length(fastqs) >= 1)
    for (i in seq(1,(length(fastqs)/2),2)) {
        ## Same sample?
        stopifnot( Fastq2Samp(fastqs[i]) == Fastq2Samp(fastqs[i+1]) )
        stopifnot( Fastq2Index(fastqs[i]) == 1  &&  Fastq2Index(fastqs[i+1]) == 2)
    }
    cat("*Assuming* that reads in paired FASTQ are in the same order.\n")
}
CheckFastq(fastqs)

## Helpers to get the forward (R1) and the reverse (R2) fastqs.
FwdFastqIdx <- function(fastqs) { seq(1,length(fastqs),2) }
RevFastqIdx <- function(fastqs) { seq(1,length(fastqs),2)+1 }

if (specialParams$useOnlyR1Reads) {
    ## If using only the R1s, remove the R2s.  Will propagate to filteredFastqs.
    fastqs <- sort(list.files(dataDir, pattern=".*_R1.*fastq.gz", full.names=TRUE, recursive=TRUE))
    ## Also redefine Fwd/RevFastqIdx.
    FwdFastqIdx <- function(fastqs) { 1:length(fastqs) }
    RevFastqIdx <- function(fastqs) { NULL }
}

filteredFastqs <- file.path(dada2OutDir, "filtered", basename(fastqs))
stopifnot(length(filteredFastqs) >= 1)

dereplicatedRdsFile <- file.path(rObjsDir, "dereplicated.rds")
errorsRdsFile       <- file.path(rObjsDir, "errors.rds")
ddsRdsFile          <- file.path(rObjsDir, "dds.rds")
mergedRdsFile       <- file.path(rObjsDir, "merged.rds")
trackRdsFile        <- file.path(rObjsDir, "track.rds")
seqTabRdsFile       <- file.path(rObjsDir, "sequenceTab.rds")
seqTabNoChimRdsFile <- file.path(rObjsDir, "sequenceTab.noChimera.rds")
bimsRds             <- file.path(rObjsDir, "chimera.rds")
asvsAbundTxt        <- file.path(txtsDir,  "asvs.tsv")
asvsNoChimAbundTxt  <- file.path(txtsDir,  "asvs.noChimera.tsv")
asvsFastaTxt        <- file.path(txtsDir,  "asvs.fasta")
asvsNoChimFastaTxt  <- file.path(txtsDir,  "asvs.noChimera.fasta")
asvsNoChimAbundLenFreqTxt <- file.path(txtsDir, "asvs.noChimera.len_freq.tsv")

## Maps fwd or rev fastq file to sample name.
sampleNames <- sapply(fastqs, Fastq2Samp)

cat("Will analyze",length(sampleNames),"paired-end fastqs from these samples:\n")
unique(sampleNames)


## Makes the logs look nicer.
BigStep = function(desc)
 {
    cat("\n")
    cat("################################################################################\n")
    cat("##\n")
    cat("## ",desc,"\n")
    cat("##\n")
}


################################################################################
BigStep("Quality plots and trim/filter the reads...")

## Note that removePrimers() was for CCS data and the docs suggest using
## external tools for primer trimming such as cutadapt or trimmomatic.  Also,
## the DADA2 tutorial presumes that primers have already been removed. So I
## created runCutadapt.sh.

## Quality score vs. position for the reads. We could use this info for trimming
## (filterAndTrim(...trunLen=...)) reads when the quality crashes,
qpdir <- file.path(plotsDir,'ReadQuality')
if (!dir.exists(qpdir)) { dir.create(qpdir) }
for (i in FwdFastqIdx(fastqs)) {
    samp <- Fastq2Samp(fastqs[i])
    plotFile = file.path(qpdir, paste0(samp,'_quality.pdf'))
    if (!file.exists(plotFile)) {
        cat(" -- Plotting quality for",samp,"fastq files.\n")
        g <- plotQualityProfile(fastqs[i:(i+1)])
        cap <- paste0('At each read position: Gray heat map = freq of QS; ',
                      'Green line = mean QS; Orange line (dashed) = median QS (quartiles)')
        g <- g + labs(caption = cap)
        ggsave(plotFile, units='in')
        rm(g,cap)
    }
}
rm(qpdir,i,samp,plotFile)


## In this post Brian Callahan said our filtering parameters [for Ana's CCS
## data] looked reasonable and that the main goal is to not pick *bad*
## parameters (rather than to optimize):
##   https://github.com/benjjneb/dada2/issues/897
##
## FIXME: Kendra uses PEAR as with parameter "-q 19" so that two consecutive
## bases with QS < 19 cause truncation of "the rest of the read."  She also
## requires *assembled* reads to be between 300 and 400 nt for the Guam data
## (assembled, so don't do that here).
##
## Okay, let's try filtering here with the goal of later producing reads between
## around 300-400nt.  Reason is that I think the error models might be bad right
## now, and I wonder if the reason is that we are getting lots of non-nifH,
## based on the many ASVs (merged reads) that are much shorter than 300nt.  Just
## as for Ana's data I tried to make learnErrors() see only the 18S and not the
## ITS, perhaps the short sequences are non-nifH and are thus introducing true
## variability that learnErrors() could mistake as sequencing error.
## So... mergePairs() tells me the number of overlapping bases (nmatch, since
## nmismatch and nindel seem to always be 0), and the ASV which lets me estimate
## the length of the paired reads that produced the ASV.  Just assume that the
## merge happened as follows with R1/2 only of equal length:
##        R1 only            overlap              R2 only
##        ----------------[-------------]----------------
## Looking at the first four samples, if I want ASV with lengths 300-400nt, then
## I should pick reads that have 162 <= len <= 284.
##     readLen  =  nmatch + (asvLen - nmatch)/2
##
if (!all(file.exists(filteredFastqs))) {
    ## Missing some filtered fastqs (whether R1 or R2 -- see above).
    cat("Quality filtering the reads. Will truncate reads at the first position",
        "with a quality score Q <=",filterAndTrimParams$truncQ,".\n")
    cat("Will drop reads that have any uncalled bases or that have fewer than",
        filterAndTrimParams$minLen,"nt.\n")
    cat("Will drop reads with more than", filterAndTrimParams$maxEE.fwd,
        "expected errors in the forward read.\n")
    if (filterAndTrimParams$truncLen.fwd > 0) {
        cat("Will chop forward reads at position",
            filterAndTrimParams$truncLen.fwd, "\n")
    }
    if (!specialParams$useOnlyR1Reads) {
        cat("Will drop reads with more than", filterAndTrimParams$maxEE.rev,
            "expected errors in the reverse read.\n",
            "Note that both reads are rejected if either the forward or",
            "reverse is rejected.\n")
        if (filterAndTrimParams$truncLen.rev > 0) {
            cat("Will chop reverse reads at position",
                filterAndTrimParams$truncLen.rev, "\n")
        }
    }
    fwdIdx <- FwdFastqIdx(filteredFastqs)
    filtFs <- filteredFastqs[fwdIdx]
    names(filtFs) <- Fastq2Samp(filtFs)
    if (!specialParams$useOnlyR1Reads) {
        ## Use the forward and reverse reads.
        revIdx <- RevFastqIdx(filteredFastqs)
        filtRs <- filteredFastqs[revIdx]
        names(filtRs) <- Fastq2Samp(filtRs)
        track <- filterAndTrim(fwd = fastqs[fwdIdx], filt     = filtFs,
                               rev = fastqs[revIdx], filt.rev = filtRs,
                               truncQ = filterAndTrimParams$truncQ,
                               truncLen = c(filterAndTrimParams$truncLen.fwd,
                                          filterAndTrimParams$truncLen.rev),
                               maxN   = 0,        # Default. Tolerate no uncalled positions.
                               maxEE  = c(filterAndTrimParams$maxEE.fwd,
                                          filterAndTrimParams$maxEE.rev),
                               minLen = filterAndTrimParams$minLen,
                               matchIDs = TRUE,   # only output reads that are paired.  FIXME FIXME FIXME: Ids do not match for Danish straits data.  Need a flexible solution...
                               multithread=TRUE)  # Set F if errors occur (see Details)
    } else {
        ## Only filter based on the forward reads.
        track <- filterAndTrim(fwd      = fastqs[fwdIdx],
                               filt     = filtFs,
                               truncQ   = filterAndTrimParams$truncQ,
                               truncLen = filterAndTrimParams$truncLen.fwd,
                               maxN     = 0,
                               maxEE    = filterAndTrimParams$maxEE.fwd,
                               minLen   = filterAndTrimParams$minLen,
                               multithread=TRUE)
    }

    cat("Here are the numbers of reads input and retained:\n")
    print( data.frame(track,
                      PctRetained = round(100*track[,'reads.out']/track[,'reads.in'],1)) )
    saveRDS(track, file=trackRdsFile)
} else {
    cat("Already filtered.\n")
}


plotFile = file.path(plotsDir, 'hist.readLengthHistFiltered.pdf')
if (!file.exists(plotFile)) {
    ## These are the lengths of reads after filtering.
    lens.fn <- lapply(filteredFastqs, function(fn) nchar(getSequences(fn)))
    lens <- do.call(c, lens.fn)
    qplot(lens, geom="histogram", binwidth = round(sd(lens)),
          main='Quality-trimmed reads', xlab = 'read length (nt)') +
        theme_bw()
    ggsave(plotFile, width=5, height=3, units='in')
    cat("Stats for read lengths after filtering:\n")
    print(summary(lens))
    rm(lens.fn,lens,plotFile)
}


################################################################################
BigStep("Dereplicating")

## Note that this should not impact error learning because learnErrors()
## dereplicates each sample (peek at the code).  This ~answers my question on
## why the fecal markdown script dereplicates and then learns errors but the big
## data example learns errors without 1st dereplicating.
##
if (!file.exists(dereplicatedRdsFile)) {
    cat("Dereplicating...\n")
    system.time({  dereplicated <- derepFastq(filteredFastqs, verbose=TRUE)  })
    saveRDS(dereplicated, file=dereplicatedRdsFile)
} else {
    cat("Already dereplicated, loading...\n")
    dereplicated <- readRDS(dereplicatedRdsFile)
}

## 'dereplicated' is only used by LearnErrorsFromSubsampledReads().  This script
## changed to do derepFastq() in the loop that does sample inference, so that we
## can free up the huge amount of memory needed for 'dereplicated', at a very
## small performance hit for derepliating each sample again (very fast).
##
cat("Note that dereplicated reads at this point helps us make the error models",
    "in the next step.  Later, we dereplicate again during sample inference.\n")


################################################################################
BigStep("Making error models")

## This function was created to speed up error modeling by limiting the total bp
## fed to learnErrors().  Basically the idea worked -- we got similiar models
## for Ana's PacBio data with only small improvements (visually in the error
## model panel plot) if we learned from 1K, 2K, 3K, 6K, or 10K PacBio reads.
## I changed the function to take nBasesPerSample rather than seqsPerSample
##
## The DADA2 big data example says 100M total base pairs is more than adequate
## for learning.  And looking at the learnErrors() code, it breaks out of the
## samples loop once it has seen 'nbases' which defaults to 100M.  So
## LearnErrorsFromSubsampledReads() has the added benefit of distributing the
## learning over more the samples (until the total nt wanted for learning are
## seen).
##
## Note that the actual call to learnErrors() was changed for MiSeq. Various
## PacBio specific params were dropped.
##
## Caution: LearnErrorsFromSubsampledReads() relies on a few globals (including
## 'dereplicated').
##
LearnErrorsFromSubsampledReads = function(nBasesForLearning)
{
    ## For each sample, pick a random subset of the dereplicated reads.
    set.seed(123)  # For repeatability
    forLearningErrorsFastqs <- file.path(dada2OutDir, "forlearningerrors",
                                         basename(filteredFastqs))
    stopifnot(length(forLearningErrorsFastqs) >= 1)
    names(forLearningErrorsFastqs) = basename(filteredFastqs)
    if (!dir.exists(dirname(forLearningErrorsFastqs[1]))) {
        dir.create(dirname(forLearningErrorsFastqs[1]))
    }

    ## Use mean read length in the first fastq to estimate how many reads we
    ## need from each sample to attain nBasesForLearning.
    seqsPerSample <- (nBasesForLearning / mean(width(readFastq(filteredFastqs[1])))) / length(filteredFastqs)
    seqsPerSample <- round(1.1 * seqsPerSample)
    cat('Will randomly select',seqsPerSample,'reads from each sample so that in total we',
        'feed provide about',nBasesForLearning,'nt to learnErrors().\n')
    for (fq in filteredFastqs) {
        bNam <- basename(fq)
        ## Pick from unique sequences.  Uniques vector has names that are sequences and
        ## values that are abundances (docs online).  If too few reads, sample() errs.
        uniqs <- unique(dereplicated[[basename(fq)]][['map']])
        if (length(uniqs) >= seqsPerSample) {
            cat("Sampling",seqsPerSample,"random sequences from the",length(uniqs),"uniques in",fq,"\n")
            idx <- sample(uniqs, seqsPerSample, replace=F)
        } else {
            cat("Using all",length(uniqs),"sequences from",fq,"\n")
            idx <- uniqs
        }
        srq = readFastq(dirname(fq), pattern=basename(fq))[idx]  # Class ShortReadQ
        ofile = forLearningErrorsFastqs[bNam]
        if (file.exists(ofile)) { unlink(ofile) }
        writeFastq(srq, file=forLearningErrorsFastqs[bNam], compress=T)
    }

    ## For MiSeq.  (Adapted from Ana's script which had special handling for
    ## PacBio including errorEstimationFunction=PacBioErrFun and BAND_SIZE=32
    ## [which has no documentation?], as in Callahan 2018).
    errors <- learnErrors(forLearningErrorsFastqs, nbases=nBasesForLearning,
                          multithread=TRUE, verbose=TRUE)
    return(errors)
}


## Removed several functions for learning error rates from 18S or 5.8S parts of
## Ana's PacBio reads -- so that highly variable ITS regions would not confuse
## DADA2 about errors vs. true variation.  Right now we don't need this approach
## for nifH, but let's keep it in mind... (since as of Mar 24, 2021 I don't like
## the error modes, they rates never get much below 1%.).



## Learn errors from a random sample of the reads.
if (!is.na(preCalculatedErrorModel)) {
    cat("Will use pre-calculated error model in",preCalculatedErrorModel,"\n")
    stopifnot(file.exists(preCalculatedErrorModel))
    errors <- readRDS(preCalculatedErrorModel)
    stopifnot(is.list(errors) && is.matrix(errors$trans))
    saveRDS(errors, errorsRdsFile)
    cat("Saved the pre-calculated model to",errorsRdsFile,".\n")
} else if (!file.exists(errorsRdsFile)) {
    cat("Start learning errors at",date(),"\n")
    ## Apr 2020:  See older script for code that used to learn error models based on
    ## sampling 100, 500, ... MAX_READS_TO_SAMPLE reads, when we were skeptical about
    ## whether our error models could be trusted [for Ana's PacBio data].
    system.time({  errors <- LearnErrorsFromSubsampledReads(NBASES_FOR_LEARNING_ERRORS)  })
    msg <- paste('Error models based on',NBASES_FOR_LEARNING_ERRORS,
                 'nt from each sample')
    plotErrors(errors, nominalQ=T) + ggtitle(msg)
    saveRDS(errors, errorsRdsFile)
    ggsave(file.path(plotsDir,'errorModels.pdf'))
    cat("Finshed learning errors at",date(),"\n")
} else {
    errors <- readRDS(errorsRdsFile)
    cat("Already had error model.\n")
}


################################################################################
BigStep("Denoising -- identifying ASV's")

cat('Deleting object dereplicated to free some memory.\n')
cat('Will recalculate the dereplicated reads for each sample as needed.\n')
rm(dereplicated)

if (!file.exists(ddsRdsFile)) {
    ## Per big data suggestions page/example at DADA2 site, do by sample.  (Note
    ## that they suggest by sample despite 'pool=FALSE' by default.)
    dds <- vector("list", length(filteredFastqs))
    names(dds) <- filteredFastqs
    tmpDdsFile <- file.path(rObjsDir, 'ddsInProgress.rds')  # In case crash, have something
    for(sam in filteredFastqs) {
        cat("Processing", sam, "at",date(),"\n")
        cat("First, dereplicating (again)...\n")
        derep <- derepFastq(sam, verbose=T)
        ## derep knows the abundances of the unique sequences, so don't worry
        ## that dada() won't know relative abundances.
        cat("Now, denoising...\n")
        ## FIXME: (1) Dropped BAND_SIZE; (2) Should I 'pool'?
        cpuTime = system.time({ dds[[sam]] <- dada(derep, err=errors, multithread=TRUE) })
        cat("Time spent in dada():\n")
        print(cpuTime)
        rm(derep) # so not in memory during next derepFastq()
        saveRDS(dds, file=tmpDdsFile)
        cat("Finished ", sam, "at",date(),"\n")
    }
    saveRDS(dds, ddsRdsFile)
    unlink(tmpDdsFile)
} else {
    dds <- readRDS(ddsRdsFile)
    cat("Already had DADA results.\n")
}



track <- readRDS(trackRdsFile)
cat("Here's how the number of reads have changed after each of the previous steps.\n")
## FIXME: Would be nice to show the initial read count.  Unlike for Ana I
## don't have data file for primer removal though.
##df = cbind(initial=???, filtered=track[,2],
df = cbind(filtered=track[,2], denoised=sapply(dds, function(x) sum(x$denoised)))
rownames(df) <- sapply(rownames(df), function(fq) {
    paste0(Fastq2Samp(fq), '_R', Fastq2Index(fq))
})

processedTable <- df  # Need this for final table.
print(processedTable)
df = melt(df)
colnames(df) = c('Sample','step','reads')
ggplot(df, aes(x=step, y=reads, color=Sample)) + geom_line(aes(group=Sample)) + geom_point() +
    ggtitle('Reads processed by dada2 by sample') +
    xlab('Processing step') + ylab("Number of reads") +
    geom_text(aes(label=Sample), data=subset(df, step %in% c('initial.ccs')))
#    theme(legend.position="none")
ggsave(file.path(plotsDir,'readsProcessed.pdf'), width=7.5, height=5, units='in')


################################################################################
BigStep("Merging denoised paired reads")

if (specialParams$useOnlyR1Reads) {
    cat("Skipping read merging because useOnlyR1Reads is TRUE.\n")
} else {
    ## The names of dds are the file faths to the filtered reads.
    fidx <- FwdFastqIdx(names(dds))
    ridx <- RevFastqIdx(names(dds))
    ## Verify correspondence.
    stopifnot(Fastq2Samp(names(dds)[fidx]) == Fastq2Samp(names(dds)[ridx]))
    stopifnot(Fastq2Index(names(dds)[fidx])+1 == Fastq2Index(names(dds)[ridx]))
    if (!file.exists(mergedRdsFile)) {
        mergers <- mergePairs(dds[fidx], names(dds)[fidx],
                              dds[ridx], names(dds)[ridx],
                              minOverlap      = mergePairsParams$minOverlap,
                              maxMismatch     = mergePairsParams$maxMismatch,
                              justConcatenate = mergePairsParams$justConcatenate,
                              verbose=T)
        ## Convert mergers to a list of data.frames if necessary (b/c just 1 seq run).
        if (is.data.frame(mergers)) {
            stopifnot(length(fidx)==1)
            mergers <- list(mergers)
            names(mergers) <- names(dds)[1]
        }
        ## Make the names not mention R1.  (Currently they are the R1 fastq paths.)
        names(mergers) <- Fastq2Samp(names(mergers))
        saveRDS(mergers, mergedRdsFile)

        ## Make a list of unique ASV sequences by enumerating the df's in the list.
        ## (Could skip -- see makeSequenceTable() below.)
        numAsvs <- length(unique(unlist(lapply(mergers, function(v) v$sequence))))
        cat("mergePairs() created",numAsvs,"ASVs in",length(mergers),"samples.\n")
        if (mergePairsParams$justConcatenate) {
            cat("mergePairs() was asked to concatenate, so each ASV will have a 5' part\n",
                "and a 3' part that are joined by 10 N's.  The 3' part will have been reverse\n",
                "complemented and it might or might not share sequence with the 5' part.\n\n")
        } else {
            ## Provide some stats on the merge.
            pre  <- sum(sapply(mergers, function(v) sum(v[,'abundance'])))
            post <- sum(sapply(1:length(dds), function(i) sum(dds[[i]][['denoised']])))
            cat("The merge tolerated at most",mergePairsParams$maxMismatch,"mismatches between\n",
                "the forward and reverse ASVs.  In total", post, "read pairs correspond to\n",
                "post-merge ASVs compared to", post, "reads (fwd and rev separate) in pre-merge\n",
                "ASVs.\n")
            cat("Below for each sample are stats for the top 10 post-merge ASVs:\n",
                "  Num reads (abundance), and num matches, mismatches, and indels in\n",
                "  the overlapping regions.\n")
            for (i in 1:length(mergers)) {
                cat("\nSample",names(mergers)[i],"\n")
                print(head(mergers[[i]][,-1], n=10))
                cat("\n")
            }
            cat("\n\n")
        }
    } else {
        mergers <- readRDS(mergedRdsFile)
        cat("Already had merged reads.\n")    
    }
    rm(fidx,ridx)
}


################################################################################
BigStep("Making sequence table") # used to be part of denoising

cat("Making sequence table.  Here are the dimensions (samples X ASVs)\n")
if (!specialParams$useOnlyR1Reads) {
    sequenceTab <- makeSequenceTable(mergers)
} else {
    ## There was no merge so make the sequence table from the dds.
    ## As in the merge step, simplify the sample names (but afterwards).
    sequenceTab <- makeSequenceTable(dds)
    rownames(sequenceTab) <- Fastq2Samp(rownames(sequenceTab))
}
dim(sequenceTab)
saveRDS(sequenceTab, file=seqTabRdsFile)

## Caution. I want to assign ID to the ASVs (rather than use their DNA sequence
## as DADA2 does). But I must be careful for the ID's to persist, specifically
## across the removal of chimera below.  So here make a persistent map of DNA
## sequences to ASV IDs, and the reverse map.
asvSeq2Id <- paste0('ASV.', 1:ncol(sequenceTab))
names(asvSeq2Id) <- colnames(sequenceTab)
asvId2Seq <- names(asvSeq2Id)
names(asvId2Seq) <- as.character(asvSeq2Id)

cat("Before removing chimeras, the sequence table has", ncol(sequenceTab),
    "ASV's and", nrow(sequenceTab), "samples.\n")

cat("Saving ASV abundances (", asvsAbundTxt,") ",
    "and fasta file (", asvsFastaTxt, ")\n")
df = data.frame(t(sequenceTab))                       # Make rows the ASV ID's
rownames(df) = as.character(asvSeq2Id[rownames(df)])
colnames(df) <- rownames(sequenceTab)                 # Use original (R-unfriendly) sample names
write.table(df, file=asvsAbundTxt, sep="\t", quote=F)
write(paste0('>',as.character(asvSeq2Id),"\n",names(asvSeq2Id)), file=asvsFastaTxt)

if (nrow(df) < 5 || ncol(df) <= 1) {
    cat("Not doing NMDS. Too few ASVs and/or samples.\n")
} else {
    ## NMDS of samples by their ASV profiles.
    df <- df[,which(apply(df, 2, function(v) !all(v==0)))]   # remove empty samples
    df <- df[which(apply(df,  1, function(v) !all(v==0))),]  # remove empty ASVs
    ## Normalize sequencing depths. At least need to if we use Bray-Curtis dissimilarities
    ## because B-C is affected by sampling sizes.
    df <- decostand(t(df), method='total')
    dmeth <- 'bray'  # euclidean, canberra, jaccard, gower, ...
    nmds <- metaMDS(df, dmeth, k=2, autotransform=F)         # autotransform=F b/c used decostand()
    stress = round(nmds$stress,3)
    cat("NMDS of",dmeth,"distances between sample ASV abundances had stress", stress, "\n")
    df = data.frame(nmds$points, rownames(nmds$points))
    colnames(df) = c('x','y','Sample')
    if (length(unique(df$Sample)) <= 15) {
        ## Not so many samples that the legend will explode.
        g <- ggplot(df, aes(x,y,color=Sample))
    } else {
        g <- ggplot(df, aes(x,y))
    }
    g <- g + geom_point(size=3) + coord_fixed(ratio=1) +
        labs(title = "Samples represented by ASV abundances",
             caption = paste("NMDS on",dmeth,"distances; stress =", stress),
             x=NULL, y=NULL) +
        ##geom_text_repel(aes(label=Sample), size=3) +  # FIXME: Can crash here with viewport of 0 dimensions.
        theme(legend.position="none") + theme_bw()
    ggsave(file.path(plotsDir,'asvsNMDS.pdf'), width=4.5, height=4.5, units='in')
    rm(nmds, dmeth)
}
rm(df)

## Info for the most abundant ASV's in each sample
asvTab = t(sequenceTab)
stopifnot(rownames(asvTab) %in% names(asvSeq2Id)) # every sequence has an ID
rownames(asvTab) <- as.character(asvSeq2Id[rownames(asvTab)])   # ASV.n rather than ASV sequence
colnames(asvTab) = sapply(strsplit(basename(colnames(asvTab)), "_"), `[`, 1)
maxAsv = apply(asvTab, 2, max)
maxAsvId <- sapply(names(maxAsv),
                   function (samp) {
                       idx <- match(maxAsv[samp], asvTab[,samp])[1] # unlikely ties
                       rownames(asvTab)[idx]
                   })
maxAsvPct <- round(100 * maxAsv / colSums(asvTab), 1)
cat("Here are the most abundant ASV's in each sample (ignoring ties in a sample)\n")
data.frame(ASV.id=maxAsvId, Abund=maxAsv, Pct=maxAsvPct)

## Pretty density plot of ASV log abundances
df = read.table(asvsAbundTxt, sep="\t")
if (nrow(df) < 5) {
    cat("Skipping ASV density plot. Too few ASVs.\n")
} else {
    df = melt(df, id=NULL)
    colnames(df) = c('Sample','Abundance')
    df = df[df$Abundance > 100,]
    ggplot(df, aes(x=Abundance, color=Sample)) + geom_density() + # aes(alpha=0.15)
      scale_x_continuous(trans='log10') + ggtitle('ASV abundances') +
      xlab('Abundance for ASVs with >100 amplicons')
    ggsave(file.path(plotsDir,'asvsDensity.pdf'), width=7.5, height=5, units='in')
}
rm(df)


################################################################################
BigStep("Assigning taxonomy -- DROPPED for nifH")

cat("Not supported.  We could use DADA2's assignTaxonomy() which uses a naive Bayes classifier.\n",
    "But for nifH we have better ways.\n")


################################################################################
BigStep("Identifying chimera")

if (!file.exists(bimsRds)) {
    cat("Looking for chimeric ASV's, as in Callahan et al. 2018...")
    system.time({  bims <- isBimeraDenovo(sequenceTab,
                                          minFoldParentOverAbundance=3.5,
                                          multithread=TRUE)  })
    saveRDS(bims, file=bimsRds)
    cat("Done!\n")
} else {
    cat("Loading previously saved taxonomic assignments.")
    bims <- readRDS(bimsRds)
}

cat("How many of the ASVs are probably chimeric (TRUE) vs. not (FALSE):")
print(table(bims))

cat("What % of the total abundance is from chimera? ",
    100*sum(sequenceTab[,bims])/sum(sequenceTab), "%\n")

if (exists('taxa') && !is.null(taxa)) {
    cat("This table summarizes the genera associated with chimeric ASVs")
    x = sort(table(taxa[bims,6]), decreasing=T);           x.u = sum(is.na(taxa[bims,6]))
    y = sort(table(taxa[    ,6]), decreasing=T)[names(x)]; y.u = sum(is.na(taxa[    ,6]))
    df = data.frame(as.vector(x), as.vector(y))
    df = rbind(df, c(x.u, y.u))
    rownames(df) = c(names(x), 'unknown')
    colnames(df) = c('chimeric.ASVs','total.ASVs')
    print(df)
    rm(x,y,df)
}

cat("Throwing out chimeric ASVs\n")  # Not sure why Callahan did not.
sequenceTab <- sequenceTab[,!bims,drop=F]
cat("After removing chimeras, the sequence table has", ncol(sequenceTab),
    "ASV's and", nrow(sequenceTab), "samples. Saving to",seqTabNoChimRdsFile,"\n")
saveRDS(sequenceTab, file=seqTabNoChimRdsFile)
cat("Saving FASTA of non-chimeric ASVs as",asvsNoChimFastaTxt,"\n")
## Make sure we retain the original ASV ID's, even if there will be holes for
## the chimeric ASV's we throw out.
stopifnot(colnames(sequenceTab) %in% names(asvSeq2Id))
write(paste0('>',as.character(asvSeq2Id[colnames(sequenceTab)]),"\n",colnames(sequenceTab)),
      file=asvsNoChimFastaTxt)
cat("Saving abundance table of non-chimeric ASVs as",asvsNoChimAbundTxt,"\n")
df = data.frame(t(sequenceTab))  # rows are now ASV sequences
rownames(df) <- as.character(asvSeq2Id[rownames(df)])
colnames(df) <- rownames(sequenceTab) # insist on sample names even if not R-friendly
write.table(df, file=asvsNoChimAbundTxt, sep="\t", quote=F)

## 2021 Mar: At this point we have high quality ASVs (that had both primers,
## passed the quality filtering, etc.), but there is a wide range of lengths.
## What are the mysterious ASVs that start/end with nifH primers but have
## non-typical nifH length?  At least make a table of these (and all) ASVs so
## that we can follow up on them.
## Pass in the ASV abundance tsv and its corresponding the ASV fasta.
MakeTableOfASV_Len_Freq <- function(asvAbund, asvFasta)
{
    tab <- read.table(asvAbund, sep='\t')
    asvIds <- rownames(tab)
    fas <- readFasta(asvFasta) # class ShortRead
    if (!all(id(fas) == asvIds)) {
        stop('The asvAbund table and asvFasta do not have the matched ASV identifiers.',
             'ASVs might be missing or they are not in the same order.')
    }
    df <- data.frame(asvLen=width(fas), totReads=rowSums(tab),
                     numSamps=apply(tab,1,function(v) sum(v>0)))
    rownames(df) <- asvIds
    df
}

cat("Writing table of ASVs that shows their lengths, frequencies, and num samples in",
    "which they appear: ", asvsNoChimAbundLenFreqTxt,"\n",
    "This will help you follow up on abundant ASVs that have non-nifH-like lengths.\n")
df <- MakeTableOfASV_Len_Freq(asvsNoChimAbundTxt, asvsNoChimFastaTxt)
write.table(df, file=asvsNoChimAbundLenFreqTxt, sep="\t", quote=F)


## Show earlier table but now with ASV abundances.
## FIXME: processedTable has separate R1 and R2 lines but 'df'
## has the merged counts, so this step fails.

cat("Here again are the numbers of reads after each step, no including after chimera removal.\n")
df <- data.frame(processedTable, deChimeraed = colSums(df)[rownames(processedTable)])
print(df)
## Overwrite the previous plot of reads processed.
df$Sample = rownames(df)
df <- melt(df, id='Sample')
colnames(df) = c('Sample','step','reads')
ggplot(df, aes(x=step, y=reads, color=Sample)) + geom_line(aes(group=Sample)) + geom_point() +
    ggtitle('Reads processed by dada2 by sample') +
    xlab('Processing step') + ylab("Number of reads") +
    geom_text(aes(label=Sample), data=subset(df, step %in% c('initial.ccs')))
#    theme(legend.position="none")
ggsave(file.path(plotsDir,'readsProcessed.pdf'), width=7.5, height=5, units='in')

if (exists("taxa") && !is.null(taxa)) {
    taxa <- taxa[!bims,]
    cat("Saving taxa table after chimera removal to",taxaNoChimeraRds,"\n")
    saveRDS(taxa, file=taxaNoChimeraRds)
}

cat("Finished on", date(),"\n")
quit(save='no')
