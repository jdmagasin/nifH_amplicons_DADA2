#!/usr/bin/env Rscript

## Copyright (C) 2025 Jonathan D. Magasin

## Compare outputs from different runs of the DADA2 pipeline.

usageStr <- "
Make a Venn diagram showing the ASVs shared by multiple runs (at most 7) of the DADA2 pipeline.

Usage:
    comparePipelineOutputs.R  lab1=asvFasta1 lab2=asvFasta2 ... {filterType}

Paths to any ASV FASTA files can be used and the labels provided will be used in logged messages and
the Venn diagram.  (Use short labels with no white space.)  Usually that FASTAs will be
asvs.noChimera.fasta, the final FASTA outputs from the pipeline.  The script will create a Venn
diagram that shows ASVs shared among the pipeline runs.  The diagram will count ASVs but will ignore
their abundance.  However, if .tsv files are found (e.g., asvs.noChimera.tsv), then the script will
compare abundances from ASVs that are shared among the FASTAs.

Optionally the last argument can specify how to filter ASVs.  Options are:

  A:  Filter as in nifH-ASV-workflow.  Only count ASVs that are detected (>= 1 read) in at least two
      samples or that have >= 1000 reads in at least one sample.  Discard ASVs if not 281-359 nt.
  B:  Only count ASVs with >= 10 reads in at least one sample.  Discard ASVs if not 281-359 nt.
  C:  Only count ASVs with >= 100 reads in at least one sample.  Discard ASVs if not 281-359 nt.
  
"

args <- commandArgs(T)
if (length(args) == 0 || args[1] %in% c("-h","--help")) {
    cat(usageStr)
    stop("Please provide at least one directory.")
}

library(venn)

FiltLikeWorkflow <- function(abundVec) { (sum(abundVec >= 1) >= 2) || any(abundVec >= 1000) }
Filt10           <- function(abundVec) { any(abundVec >= 10) }
Filt100          <- function(abundVec) { any(abundVec >= 100) }
LenFilt <- function(asvSeqs) {
    len <- nchar(asvSeqs)
    which((280 < len) & (len < 360))
}

ftype2func <- list(A = FiltLikeWorkflow, B = Filt10, C = Filt100)
filterFunc <- NULL
if (args[length(args)] %in% names(ftype2func)) {
    cat("Filter type", args[length(args)], "selected.\n")
    filterFunc <- ftype2func[[args[length(args)]]]
    args <- args[-length(args)]
}
nicknames <- sapply(strsplit(args, '='), '[', 1)
if (length(nicknames) <= 1) {
    cat(usageStr)
    stop("Please provide file paths to at least two FASTAs.")
}
fastas <- sapply(strsplit(args, '='), '[', 2)
if (!all(file.exists(fastas))) {
    if (any(file.exists(nicknames))) {
        stop("You must provide a label for every FASTA, e.g., ",
             "\"RunA=path/to/run/A.fasta\".")
    }
    stop("The following FASTAs do not exist:\n\t",
         paste(fastas[!file.exists(fastas)], collapse = "\n\t"))
}
names(fastas) <- nicknames


## Load indicated FASTA and abundance tsv. Filter both using the passed function.
## Return a list with the sequences and abundance table, the former limited to sequences
## in the latter.
LoadAsvsAndAbunds <- function(fas, tsv, filtFunc = NULL)
{
    abundTab <- NULL
    ## Okay if caller passes a tsv that does not exist.
    if (file.exists(tsv)) {
        cat("Found abundances (.tsv) for", fas, "\n")
        abundTab <- read.table(tsv)
    }
    if (!is.null(abundTab) && !is.null(filtFunc)) {
        idx <- names(which(apply(abundTab, 1, filtFunc)))
        abundTab <- abundTab[idx,]
    }

    asvs <- readLines(fas)   # single-line FASTAs only!, no 
    idx <- grep("^[ \t]*$", asvs)
    if (length(idx) > 0) { asvs <- asvs[idx] }
    idx <- grep('^>', asvs)
    asvs <- data.frame(id = sub('^>', '', asvs[idx]), seq = asvs[-idx])
    if (!is.null(abundTab)) {
        asvs <- asvs[which(asvs$id %in% rownames(abundTab)),]
        idx <- match(asvs$id, rownames(abundTab))
        stopifnot(!is.na(idx))
        asvs <- asvs[idx,]
    }
    if (!is.null(filtFunc)) {
        ## If there is a filter function, then we should also filter ASVs based
        ## on length.  See usageStr.
        asvs <- asvs[LenFilt(asvs$seq),]
        if (!is.null(abundTab)) {
            abundTab <- abundTab[match(asvs$id, rownames(abundTab)),]
        }
    }
    ## Return asv and abundTab with same ASVs in same order
    list(asvs = asvs, abunds = abundTab)
}


## For very intersection in vennObj, get the % of total abundance comprised by the ASVs in
## the intersection relative to each sample in the intersection.
GetPctAbundInEachIntersection <- function(vennObj, datList)
{
    isects <- attr(vennObj, 'intersection')
    sapply (names(isects), function (inam) {
        asvSeqs <- isects[[inam]]
        ## For each sample snam, get % abundances for the ASVs in this intersection.
        sapply(strsplit(inam, ":")[[1]], function (snam) {
            d <- datList[[snam]]
            idx <- match(asvSeqs, d$asvs[,"seq"])
            stopifnot(!is.na(idx))
            asvIds <- d$asvs[idx,"id"]
            ifelse(length(asvIds) == 0, 0, 100*sum(d$abunds[asvIds,])/sum(d$abunds))
        })
    }, simplify = F, USE.NAMES = T)
}


datList <- lapply(fastas, function (fas) {
    tsv <- sub('\\.fasta$', '.tsv', fas)  # okay if does not exist
    LoadAsvsAndAbunds(fas, tsv, filterFunc)
})
names(datList) <- names(fastas)

cat("\nSummary of pipeline runs:\n")
x <- lapply(datList, function (d) {
    tot <- sum(d$abunds)  # NA if NULL
    list(ASVs = nrow(d$asvs), tot_abund = tot)
})
print(do.call(rbind, x))
cat("\n\n")

x <- lapply(datList, function (x) x$asv$seq)
png("ASV_Venn.png", width = 7, height = 7, units = "in", res = 144)
v <- venn(x, ilabels = "counts", ilcs = 1, sncs = 1, zcolor = "style", box = F)
sink("/dev/null"); dev.off(); sink()
cat("Saved ASV_Venn.png which shows how many ASVs are shared by each output from the DADA2 pipeline.\n")
cat("Abundances of the ASVs are *not* shown in the Venn diagram.  Abundances (% of total reads in each\n",
    "pipeline output) are summarized as follows:\n")

x <- GetPctAbundInEachIntersection(v, datList)
names(x) <- paste("% total abundance from ASVs that are in", sub(':', ' and ', names(x)))
print(x, digits = 2)

quit(save="no")
