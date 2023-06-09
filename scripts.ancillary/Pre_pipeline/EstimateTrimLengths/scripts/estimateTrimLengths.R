#!/bin/env Rscript

usage <- '
Usage:

  estimateTrimLengths.R  fwdPrimer  revPrimer  fastqDir  [sampSize]

For a data set with possibly many paired FASTQ files, show the expected
proportions of paired Illumina reads that will be utilized by DADA2 in the ASVs
using different combinations of forward and reverse read trim lengths. This helps
one select trim lengths to use with DADA2\'s filterAndTrim() in order to maximize
the number of reads used to create ASVs.

In typical usage of DADA2, one selects trim lengths based on inspection of FastQC
reports for forward and reverse reads. This is laborious and unsystematic if
there are dozens to hundreds of FASTQ files. The trim lengths are passed to
DADA2\'s filterAndTrim(). Trimming impacts the later merge step because the tails
of the forward and reverse ASVs must overlap when DADA2 merges them to form the
full length ASVs.  Usually one requires the overlap to be nearly perfect (via
parameters to mergePairs()). Under-trimming leads to mismatches in the overlap
while overly agressive trimming prevents overlap entirely.

Notes:
 1. This tool requires a minimal overlap of 20nt and tolerates at most 2
    mismatches in the overlap.  Small samples (<1K reads) are ignored.

 2. DADA2 runs after trimming of primers and adapters by e.g. cutadapt, but the
    trim lengths one would manually choose by looking at FastQC reports include
    the primers (and possibly adapters).  This tool evaluates trim lengths for
    the primer-less reads. I.e. it helps one pick parameters for filterAndTrim().

Inputs:

  -- Forward and reverse primers in 5\' to 3\' orientation.  IUPAC characters.

  -- Directory to search recursively for paired FASTQs.  Paired FASTQs must have
     identical names except for their suffixes which must end in "_R1.fastq" and
     "_R2.fastq", or "_R1_stuff.fastq" and "_R2_stuff.fastq".  "stuff" can be
     anything.  The "R" before the "1" or "2" can be absent. Optionally the
     FASTQs can be compressed and the names can end with ".fastq.gz".

  -- Number of paired reads to randomly sample when estimating overlaps.
     Default 1000.

Outputs: 
  -- Messages to stdout that include the top 10 combinations of R1 and R2 trim
     lengths.  The expected proportions of reads that will overlap is reported.
     The proportions are weighted by the sample sizes.

  -- trimLengthsHeatMap.png:  Heat map where each cell indicates the estimated
     percentage of reads that will overlap if the trim lengths for R1 (row)
     and R2 (column) are used.
'

args <- commandArgs(T)
if (length(args)==0) { cat(usage); stop("Need parameters.") }

primers <- c(fwd=args[1], rev=args[2])
stopifnot(is.character(primers))

fastqDir <- args[3]
stopifnot(dir.exists(fastqDir))

sampSizeParam <- args[4]
if (is.na(sampSizeParam)) { sampSizeParam <- 1000 }
cat("Estimates will be based on",sampSizeParam,"reads randomly picked from each sample.\n")

suffixPat <- '_R{0,1}[1,2](_.*|)\\.fastq(\\.gz|)$'
fastqList <- list.files(fastqDir, pattern=paste0("*",suffixPat), recursive=T, full.names=T)
## fastqList should have R1's in the odd slots and R2's in the evens.
n <- length(fastqList)
cat("Will look at",n/2,"paired FASTQs.\n")
stopifnot(n %%2 == 0) # Unpaired!
n <- n/2
stopifnot(n > 0)

cat("Loading libraries...")
suppressMessages(library(ShortRead))
suppressMessages(library(dada2))     # Just for nwalign
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))


## Estimate the average number of miscalled bases in a ShortRead, looking
## at just the first n bases, n taken from lengthVec.
AverageErrorsAcrossAllReads <- function(sr, lengthVec)
{
    ## Get qualities into a matrix of error probabilities for each base (using
    ## Phred defn). For each trim length, the expected num errors for a read is
    ## the sum of the probabilities 1:num, and then average over the reads.
    qm <- as(quality(sr), "matrix")
    qm <- 10^(-qm/10)
    x <- as.data.frame(sapply(lengthVec, function(num) sum(qm[,1:num], na.rm=T)) / nrow(qm))
    colnames(x) <- 'ExpectedErrors'
    x$TrimLen <- lengthVec
    x
}


## Count the number of fwd and rev reads that will overlap under various
## combinations of trimming lengths for the forward and reverse FASTQs.
## ** The trimming lengths are what should be passed to DADA2::filterAndTrim(),
##    so they are for reads without primers (and adapters).
##   fq1 is the forward (R1) FASTQ, fq2 the reverse.
##   fq1 should have primerFwd in 5' to 3' orientation near its start (adapters okay).
##   fq2 should have primerRev in 5' to 3' orientation near its start (adapters okay).
##   Fwd and rev reads must overlap by >= minOverlap nt and with <= maxMismatch'es
##   To save time sample sampSize paired reads.
##   
EvalOverlaps <- function(fq1, fq2,
                         trimLensFwd = seq(150,250,25), trimLensRev = seq(150,250,25),
                         minOverlap=20, maxMismatch=2, sampSize=sampSizeParam)
{
    stopifnot(minOverlap <= min(c(trimLensFwd,trimLensRev)))

    ## Get the proportion of reads that exceed each trim length so that we can
    ## scale the returned counts.  (There was a bug before where % of
    ## overlapping reads were w.r.t. sufficiently long reads, rather than
    ## w.r.t. all reads.)  Must sapply in same order as 'for' loop below.
    longR1scal <- sapply(trimLensFwd, function(tl) sum(width(fq1) > tl, na.rm=T))/length(fq1)
    longR2scal <- sapply(trimLensRev, function(tl) sum(width(fq2) > tl, na.rm=T))/length(fq2)

    ## Drop short reads. Then gather the paired reads that remain.
    fq1 <- fq1[width(fq1) > max(trimLensFwd)]
    fq2 <- fq2[width(fq2) > max(trimLensRev)]
    x <- PairUpTheReads(fq1, fq2)
    fq1 <- x$fq1;  fq2 <- x$fq2;  rm(x)

    ## Save time by working with a sample of the reads. Do not first identify
    ## unique fwd/rev combinations (like DADA2's mergePairs source), because (1)
    ## complicated; (2) quick check suggests limited data reduction -> speed up.
    widx <- sample.int(n=length(fq1), size = min(sampSize,length(fq1)), replace=F)
    fq1 <- fq1[widx];  fq2 <- fq2[widx]

    ## Revcomp the R2s. Required by CountAcceptableOverlaps().
    fq2 <- reverseComplement(fq2)

    ## Determine overlaps for every combination of fwd/rev trim lengths.
    mat <- matrix(0, nrow = length(trimLensFwd), ncol = length(trimLensRev),
                  dimnames = list(paste0('r1.',trimLensFwd),paste0('r2.',trimLensRev)))
    for (i in 1:length(trimLensFwd)) {
        len1 <- trimLensFwd[i]
        fq1.trim <- narrow(sread(fq1), start=1, end=len1)  # Will fail if read is < len1
        for (j in 1:length(trimLensRev)) {
            len2 <- trimLensRev[j]
            cat("\rCounting overlaps if trim R1 to",len1,"and R2 to",len2)
            fq2.trim <- narrow(sread(fq2), start = width(fq2) - len2 +1, end = width(fq2))
            cnt <- CountAcceptableOverlaps(fq1.trim, fq2.trim, minOverlap, maxMismatch)
            ## Scale to num reads in sample.
            mat[i,j] <- cnt * min(longR1scal[i],longR2scal[j])
        }
    }
    ## Return the counts matrix and the number of paired reads evaluated,
    ## which could be less than sampSize.
    return( list(counts=as.matrix(mat), numPairedReads=as.numeric(length(fq1))) )
}


## Helper used only by EvalOverlaps().
PairUpTheReads <- function(fq1,fq2)
{
    ## Get reads that are paired and order them identically.
    ids1 <- sub(' +.*$','',as.character(id(fq1)))  # ID precedes any whitespace
    ids2 <- sub(' +.*$','',as.character(id(fq2)))
    pairedIds <- intersect(ids1,ids2)
    widx <- which(ids1 %in% pairedIds)
    x <- 100*(length(widx)/length(ids1))
    if (x < 50) { warning("Only ",round(x,1),"% of the R1 reads are paired.") }
    fq1 <- fq1[widx]
    widx <- which(ids2 %in% pairedIds)
    x <- 100*(length(widx)/length(ids2))
    if (x < 50) { warning("Only ",round(x,1),"% of the R2 reads are paired.") }
    fq2 <- fq2[widx]
    list(fq1=fq1, fq2=fq2)
}


## Pairwise align the R1 (q) and R2 (s) reads and count the alignments
## that are long enough and with few mismatches.
CountAcceptableOverlaps <- function(q,s,minOverlap,maxMismatch)
{
    stopifnot(minOverlap > 0 && maxMismatch >= 0)

    ## Use DADA2 to pairwise align the forward (q) reverse (s, revcomp'd) reads.
    ## mapply/nwalign similar to source for DADA2::mergePairs(), but here we
    ## align raw reads (not unique fwd/rev ASVs).  Use default match, mismatch,
    ## and gap penalties, and no banding. Do not penalize for gaps at the end of
    ## the alignment!
    alns <- mapply(nwalign, s1=as.character(q), s2=as.character(s), endsfree=T, vec=T,
                   SIMPLIFY = FALSE)

    ## Calculate overlap lengths and mismatches. The overlap ends when R1's
    ## alignment sequence has only "-", and it starts where the R2 alignment's
    ## "-"'s end.
    res <- as.data.frame(t(sapply(alns, function(aln) {
    start <- attr(regexpr('^-+', aln[2]), 'match.length')
        end   <- as.numeric(regexpr('\\-+$',   aln[1]))
        len <- max(end - start, 0)
        fseq <- substr(aln[1],start+1,end-1)
        rseq <- substr(aln[2],start+1,end-1)
        ## Counts mismatches and indels.
        mm <- sum(strsplit(fseq,'')[[1]] != strsplit(rseq,'')[[1]])
        c(len=len, mm=mm)
    })))
    rownames(res) <- NULL
    ## Count overlaps that were long enough and with few mismatches.
    res <- subset(res, len >= minOverlap & mm <= maxMismatch)
    nrow(res)
}

## Might recycle this code to report the length(s) at which median quality
## of a window drops to Q.
    ## FIXME: Think this doesn't help (or make sense).
    ## We only want to consider trim lengths that exclude really low quality
    ## bases. Note that this can drop some reads entired.
    ###pc <- as(PhredQuality(trimErrorRate),"CharacterList")[[1]]  # char assoc. with error rate
    ###fq1 <- trimTailw(fq1, k=3, a=pc, halfwidth=10)
    ###fq2 <- trimTailw(fq2, k=3, a=pc, halfwidth=10)
    

## Assume same number of cycles for all FASTQs (and R1 and R2).

## FIXME: Guess a useful set of trim lengths to check. Max trim is such that
## only 55% of each read is retained, and min keeps 85% of each read. These are
## guesses that should ensure R1 and R2 reads ~always overlap.  The 85% should
## keep enough of the reads so that quality can deteriorate, because we want the
## heat map to show some bad trim length combinations in addition to good ones.
## Also assuming %'s of overlapping reads changes smoothly if we step by 15nt.
## (Too small steps and the script will take too long.)
## Can improve setting of trimLengths based on trimTailw()!!

fq1 <- readFastq(fastqList[1])
x <- round(c(0.55,0.85) * median(width(fq1)))
trimLens <- round(seq(x[1],x[2], by=15))
rm(fq1,x)
cat("\nWill try the following trim lengths: ", paste(trimLens),"\n")
cat("These lengths span 55% to 85% of the median read length for",fastqList[1],"\n")

resList <- list()  # results for each of the paired fastqs
numReadsList <- list()
eeR1 <- list()   # expected errors.
eeR2 <- list()   # expected errors.
for (i in seq(1,length(fastqList),2)) {

    nam1 <- fastqList[i]
    nam2 <- fastqList[i+1]
    stopifnot( sub(suffixPat,'',nam1) == sub(suffixPat,'',nam2) )  # List order misassumption

    cat("\rLoading paired reads for", sub(suffixPat,'',nam1))
    fq1 <- readFastq(nam1)
    fq2 <- readFastq(nam2)
    cat("\n")

    ## Drop samples with <1K reads. We don't want them influencing meanFracs.
    if (length(fq1) < 1000 || length(fq2) < 1000) {
        cat("Ignoring",nam1,"and",nam2,"because < 1K reads.\n")
        next
    }

    ## Remove primers and don't worry if not found. Apparently this also removes
    ## adapters. Docs unclear.
    fq1 <- trimLRPatterns(Lpattern=primers['fwd'], subject=fq1, Lfixed=F)
    fq2 <- trimLRPatterns(Lpattern=primers['rev'], subject=fq2, Lfixed=F)
    numReads <- mean(length(fq1),length(fq2))

    res <- EvalOverlaps(fq1, fq2, trimLensFwd=trimLens, trimLensRev=trimLens)
    ## Save matrix of the expected proportion of paired reads that will overlap.
    resList[[length(resList)+1]] <- res$counts/res$numPairedReads
    numReadsList[[length(numReadsList)+1]] <- numReads

    ## Estimate the expected number of errors in the trimmed R1 and R2 reads.
    eeR1[[length(eeR1)+1]] <- AverageErrorsAcrossAllReads(fq1, trimLens)
    eeR2[[length(eeR2)+1]] <- AverageErrorsAcrossAllReads(fq2, trimLens)
}
rm(i,nam1,nam2,fq1,fq2,numReads,res)

## Weight the expected proportion of overlapping reads for each sample by the
## sample's size (out of total reads).
totReads <- sum(unlist(numReadsList))
wts <- unlist(numReadsList)/totReads
stopifnot(sum(wts)==1)
wresList <- lapply(1:length(resList), function(i) resList[[i]]*wts[i])

## Make a data frame with each row corresponding to an R1/R2 fastqs pair (a
## 'sample') trimmed as indicated by 'TrimR1/2'. Values are the weighted
## proportion of reads from the sample that are expected to overlap.
dat <- do.call(rbind, lapply(wresList,melt))
colnames(dat) <- c('TrimR1','TrimR2','Wfrac')

## Get the expected proportion of reads that will overlap for each R1-R2 trim
## length combination. This is just an expected value: 'sum' over the weighted
## values, where the values are proporitions (from resList) and the weights are
## proportions of total reads.
expectedFracs <- aggregate(Wfrac ~ TrimR1 + TrimR2, dat, sum)

expectedFracs$Percentage <- as.character(round(100*expectedFracs$Wfrac,1))
oidx <- order(expectedFracs$Wfrac, decreasing=T)
expectedFracs$Percentage[oidx[1]] <- paste(expectedFracs$Percentage[oidx[1]],'*')

cat("\n\nThis table shows the top 10 combinations of trim lengths for maximizing the % of\n",
    "reads that will overlap.  Also shown are the expected numbers sequencing errors when\n",
    "R1 and R2 reads are trimmed as indicated. These stats will guide your selection of\n",
    "DADA2::filterAndTrim() parameters truncLen and maxEE.\n")
x <- expectedFracs[oidx[1:10],c('TrimR1','TrimR2','Percentage')]
rownames(x) <- NULL
x$TrimR1 <- sub('^r.\\.','',x$TrimR1)
x$TrimR2 <- sub('^r.\\.','',x$TrimR2)
## Tack on the expected errors for R1 and R2 reads if trimmed at each indicated combination.
## Matrixify the espected errors for R1 and R2 reads.
eeR1 <- cbind(eeR1[[1]][,'TrimLen'], do.call(cbind,lapply(eeR1, function(v) v[,1])))
eeR2 <- cbind(eeR2[[1]][,'TrimLen'], do.call(cbind,lapply(eeR2, function(v) v[,1])))
stopifnot(dim(eeR1) == dim(eeR2))
x$expErrR1.mean <- NA
x$expErrR1.sd   <- NA
x$expErrR2.mean <- NA
x$expErrR2.sd   <- NA
for (i in 1:nrow(x)) {
    evec <- eeR1[ eeR1[,1] == x[i,'TrimR1'], -1]
    x[i,'expErrR1.mean'] <- round(mean(evec),2)   # Not weighted by sample size...
    x[i,'expErrR1.sd']   <- round(sd(evec),2)
    evec <- eeR2[ eeR2[,1] == x[i,'TrimR2'], -1]
    x[i,'expErrR2.mean'] <- round(mean(evec),2)
    x[i,'expErrR2.sd']   <- round(sd(evec),2)
}
print(x)
rm(x)

##save.image('debug.RData')

g <- ggplot(expectedFracs, aes(TrimR1,TrimR2)) +
     geom_tile(aes(fill=Wfrac)) + scale_fill_gradient(low='gray', high='red') +
     geom_text(aes(label=Percentage)) +
     theme_bw() +
     labs(title='Estimated % of paired reads that will overlap',
          x = 'R1 trim length', y = 'R2 trim length')
ggsave('trimLengthsHeatMap.png',g)

cat("Done!\n")
