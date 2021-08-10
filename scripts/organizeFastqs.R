#!/bin/env Rscript

##
## Make a directory structure that organizes FASTQs in a useful way, specified
## in the input .tsv file.
##
## Additional info (mainly covered in usageStr):
##
## It is useful to organize sequencing runs input to the DADA2 pipeline (for
## example) by sample type (DNA vs. RNA, size fraction, etc):
##
##   1. Makes it easier to browse the data from the command line, including
##      short shell scripts that recover data by type. 
##
##   2. Allows the pipeline to partition the FASTQs by sample type and do
##      completely independent runs.  Then you can see whether identical ASVs
##      are found (e.g. in DNA vs. RNA pipeline runs) which is very helpful when
##      trying to determine what parameters work well with the data at hand.
##      More importantly, DNA and RNA samples have very different relative
##      abundances (e.g. rare members can be transcriptionally very active) and
##      perhaps different sequencing error characteristics, so it seems safer
##      for DADA2 to process them separately.
##
## - Each line is a path to a FASTQ file.  The paths can be full or relative, as
##   long as they are reachable from the working directory in which you run the
##   pipeline.  If "ls some/path/to/my/R1.FASTQ" works, you're fine.
##
## - Columns >=2 describe the folder hierarchy that the pipeline will create
##   to hold each FASTQ.  For example column 2 could hold "DNA" or "RNA" as
##   appropriate for each FASTQ, column 3 could hold the size fraction, etc.
##   For paired sequences the final column should be a descriptive name for
##   a directory to hold the R1 and R2 FASTQ files.  This can be whatever
##   you like, but probably you will reuse the folder name assigned by the
##   sequencing center, optionally simplified to exclude checksums or other
##   numbers that aren't biologically meaningful.
##

usageStr <- "

organizeFastqs.R

Make a directory structure that meaningfully organizes FASTQs as specfied in the
input table (.tsv).  Usually sequencing centers provide a flat directory
hierarchy e.g. with one directory that contains subdirectories for each pair of
FASTQs (R1 and R2).  Use this script to organize the FASTQs meaningfully. For
example, you can partition them into groups with specific sequence type (DNA
vs. mRNA) and/or size fraction (example below).  This is a good choice when
using this script before running the DADA2 pipeline.

Usage:
    organizeFastqs.R   fastqMap.tsv

Each line of the fastqMap.tsv maps one FASTQ file (in column 1) to a location in
a directory structure (hierarchy in columns 2 to N).  This script creates the
directory structure, and the leaf folders contain symbolic links to the original
FASTQs.

Example: Suppose the sequencing center provided our paired FASTQs in
sudirectories for each MiSeq run such as MyMiSeqRun100_65647-43208182,
MyMiSeqRun100_65647-43208888, etc. Suppose we have DNA and RNA samples
and multiple size fractions and want to reorganize first by sample type
(DNA, mRNA) and then by size fraction.  The following fastqMap.tsv (tabs
separate columns) will do that.

     some/path/to/MyMiSeqRun100_65647-43208182/MyMiSeqRun100-65647_S100_L001_R1_001.fastq.gz  DNA  Frac02  MyMiSeqRun100-65647
     some/path/to/MyMiSeqRun100_65647-43208182/MyMiSeqRun100-65647_S100_L001_R2_001.fastq.gz  DNA  Frac02  MyMiSeqRun100-65647
     some/path/to/MyMiSeqRun101_65647-43208888/MyMiSeqRun101-65647_S101_L001_R1_001.fastq.gz  RNA  Frac02  MyMiSeqRun101-65647
     some/path/to/MyMiSeqRun101_65647-43208888/MyMiSeqRun101-65647_S101_L001_R2_001.fastq.gz  RNA  Frac02  MyMiSeqRun101-65647
     ...
 
  In the current directory this script will create top directory LinksToFASTQs
  with the following structure:
     LinksToFASTQs
      +--DNA
          +--Frac02
              +--MyMiSeqRun100-65647
                  +--[link to MyMiSeqRun100-65647_S100_L001_R1_001.fastq.gz]
                  +--[link to MyMiSeqRun100-65647_S100_L001_R2_001.fastq.gz]
      +--RNA
          +--Frac02
              +--MyMiSeqRun101-65647
                  +--[link to MyMiSeqRun100-65647_S101_L001_R1_001.fastq.gz]
                  +--[link to MyMiSeqRun100-65647_S101_L001_R2_001.fastq.gz]

  This example also simplifies the names of the folders that contain the paired
  FASTQs since we do not care about the \"43208182\" attached by the sequencing
  center.

"

args <- commandArgs(trailingOnly=T)
if (length(args)==0 || grepl("-h",args[1])) {
    cat(usageStr)
    quit(save='no')
}
fastqMapTsv <- args[1]
if (!file.exists(fastqMapTsv)) {
    stop("ERROR: Missing the fastqMap.tsv.")
}

topDir <- 'LinksToFastqs'
if (dir.exists(topDir)) {
    stop("ERROR: ",topDir," exists. Rename or remove it.\n")
}

###fastqMapTsv <- 'fastqMap.tsv' # developing
fastqMap <- read.table(fastqMapTsv, header=F, comment.char='#')
stopifnot(ncol(fastqMap) >=2)  # col 1 has source and the remainder have directory parts


## Basic checks of the tsv, for paired FASTQs only.
if (grepl('_R[1,2]_',fastqMap[1,1])) {
    ## Data seems to have paired FASTQs so make sure that they all come in pairs.
    x <- table(sub('_R[1,2]_','',fastqMap[,1]))
    if (any(x!=2)) {
        cat("The following data sets seem to not have exactly one R1 and R2 file.\n")
	print( names(x)[x!=2] )
    }

    ## Make sure that for each pair of FASTQs the holding directory (final
    ## column) is the same.
    x <- paste(sub('_R[1,2]_','',fastqMap[,1]), ':', fastqMap[,ncol(fastqMap)])
    x <- table(x)
    if (any(x!=2)) {
        cat("The following FASTQs seem not to be with their paired FASTQ in ",
            "the final destination directory:\n")
        print(x[x!=2])
    }
}


cat("Making directory structure rooted at",topDir,"with symbolic links to the FASTQs.\n")

paths <- c()
for (i in 1:nrow(fastqMap)) {
    p <- fastqMap[i,-1,drop=F] # All but the source path
    p <- do.call(file.path, p)
    paths <- c(paths,p)
}
paths <- file.path(topDir, paths)
for (p in unique(paths)) {
    if (!dir.exists(p)) { dir.create(p, recursive=T) }
}

## paths still is 1:1 with the fastqs which is what we need for setting up the
## symlinks.
fpaths <- file.path(paths, basename(fastqMap[,1]))  # full symlink paths (dest)
idx <- which(!file.exists(fpaths))
cat("Will make",length(idx),"symbolic links to FASTQs.\n")
if (length(idx) > 0) {
    ## This converts relative paths to absolute or aborts. Safer so that the
    ## DADA2 pipeline can switch the current working directory and still find
    ## the FASTQs.
    srcFqPaths <- normalizePath(fastqMap[idx,1], mustWork=TRUE)
    ok <- sum(file.symlink(srcFqPaths, fpaths[idx]))
    cat("Successfully created",ok,"symlinks.\n")
}
