#!/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin

##
## Overview.  See the command line docs (or usageStr).
##
## Make a directory structure (in LinksToFASTQs) that organizes FASTQs into
## groups for processing by the pipeline.  The structure is specified in a .tsv
## file.  For example, one might want to organize sequencing runs by whether
## they are DNA vs. RNA, collection station, size fraction, etc.  Reasons for
## organizing this way include:
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
##      for DADA2 to process them separately. (Granted, the pipeline runs dada()
##      with pool=F so ASV inference happens separately for each sample. However,
##      DADA2 does not create an error model for each sample and isBimeraDenovo()
##      takes all samples at once.)
##
## Overview of input tsv file:
## - Column 1 of each line has a path to a FASTQ file.  The paths can be full or
##   relative, as long as they are reachable from the working directory in which
##   you run the pipeline.  If "ls some/path/to/my/R1.FASTQ" works, you're fine.
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

Often the samples from a study represent different locations, seasons, depths,
size fractions, and/or nucleic acid type (DNA or RNA), and consequently
different microbial communities and sequencing error profiles.  We suggest
partitioning samples along these differences and processing them separately by
DADA2.

Our DADA2 nifH pipeline supports the separate handling of different sample types
through \"processing groups\", each with its own error model and ASV inference by
dada().  Before running the pipeline, one runs organizeFastqs.R to create a
directory structure (LinksToFastqs) that describes the processing groups. Each
processing group corresponds to a path through the LinksToFastqs hierarchy. For
example, one path might be LinksToFastqs/DNA/Station23/SizeFract0.2.  (See
example below.)

Usage:
    organizeFastqs.R   fastqMap.tsv

Each line of fastqMap.tsv describes one FASTQ file.  Column 1 has the path (full
or relative) to the source FASTQ file, with separate lines for the forward and
the reverse FASTQs.  Columns 2 to N describe intended path to the FASTQ within
LinksToFastqs, with column N the name of a directory that holds symbolic links
to the source FASTQs to save disk space.

The input tsv file can be named whatever you like. Within the tsv you can
comment out (with '#' at the line start) or remove samples, e.g. bad sequencing
runs, then rerun the script.

Example: Suppose the sequencing center provided our paired FASTQs in
sudirectories for each MiSeq run such as MyMiSeqRun100_65647-43208182,
MyMiSeqRun100_65647-43208888, etc. Suppose we have DNA and RNA samples and
multiple size fractions and want to reorganize first by sample type (DNA, mRNA)
and then by size fraction.  The following input file will do that. (Note that
*tabs* separate the columns.)

     some/path/to/MyMiSeqRun100_65647-43208182/MyMiSeqRun100-65647-S100_L001_R1_001.fastq.gz  DNA  Frac02  MyMiSeqRun100-65647
     some/path/to/MyMiSeqRun100_65647-43208182/MyMiSeqRun100-65647-S100_L001_R2_001.fastq.gz  DNA  Frac02  MyMiSeqRun100-65647
     some/path/to/MyMiSeqRun101_65647-43208888/MyMiSeqRun101-65647-S101_L001_R1_001.fastq.gz  RNA  Frac02  MyMiSeqRun101-65647
     some/path/to/MyMiSeqRun101_65647-43208888/MyMiSeqRun101-65647-S101_L001_R2_001.fastq.gz  RNA  Frac02  MyMiSeqRun101-65647
     ...
 
  In the current directory this script will create top directory LinksToFASTQs
  with the following structure:
     LinksToFASTQs
      +--DNA
          +--Frac02
              +--MyMiSeqRun100-65647
                  +--[link to MyMiSeqRun100-65647-S100_L001_R1_001.fastq.gz]
                  +--[link to MyMiSeqRun100-65647-S100_L001_R2_001.fastq.gz]
      +--RNA
          +--Frac02
              +--MyMiSeqRun101-65647
                  +--[link to MyMiSeqRun100-65647-S101_L001_R1_001.fastq.gz]
                  +--[link to MyMiSeqRun100-65647-S101_L001_R2_001.fastq.gz]

  This example also simplifies the names of the folders that contain the paired
  FASTQs (in the last column) since we do not care about the \"43208182\"
  attached by the sequencing center.


** FASTQ naming conventions **

The DADA2 nifH pipeline must be able to extract sample names and sequencing read
direction ('R1' or 'R2') from the FASTQ file names.  FASTQ names should follow
this format:

    {Samp}{_stuff1}_R{1,2}{stuff2}{.fastq.gz}

where:
   - Samp is required.  It can have any character other than '_'.  Do not use
     '_' to delimit parts of your sample names because those parts will not
     appear e.g. in the sample names in the ASV abundance table.

   - stuff1, if present, can have any character (including '_') but must be
     flanked by '_'. Usually stuff1 will be the sequencing lane e.g. L001 in the
     examples above.

   - stuff2 can be anything, or absent. Most likely it will begin with '_' as in
     the example.

   - The '.fastq.gz' can be absent, but that is bad style.

Sometimes you must write a small script to change your FASTQ names to follow the
format, usually to convert '_' that appear in the sample name to some other
character.  FASTQs from the Sequencing Read Archive use the format <run
name>_1.fastq.gz.  For convenience organizeFastqs.R automatically changes _1 and
_2 to _R1 and _R2 for the symbolic links.

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
    ## Offer to delete topDir so we don't have to suggest that they 'rm -rf' it.
    cat(topDir, "exists.  You must rename or delete it to proceed.\n")
    cat("Delete",topDir,"('y' or 'n')?: ")
    if ('y' == readLines('stdin', 1)) {
        if (0 == unlink(topDir, recursive=TRUE, force=FALSE)) {
            cat(paste0("Deleted.  Will create a new ", topDir,".\n"))
        } else {
            stop("ERROR: Failed to delete", topDir,"\n")
        }
    } else {
        stop("Aborting.  You must rename or delete", topDir, "\n")
    }
}

fastqMap <- read.table(fastqMapTsv, header=F, comment.char='#')
stopifnot(ncol(fastqMap) >=2)  # col 1 has source and the remainder have directory parts

## Check for pipeline-friendly FASTQ names.
ValidFastqNames <- function(namVec)
{
    badNams <- which(!sapply(basename(namVec), grepl,
                             pattern='^[^_]+(|.*[^_])_R[12].*$'))
    if (length(badNams) > 0) {
        cat("The following FASTQ names do not follow the format required by the",
            "DADA2 nifH pipeline.  See organizeFastqs.R docs for description.\n",
            namVec[badNams], "\n")
        stop("Abborting due to bad FASTQ names.\n")
    }
}

## SRA names are not pipeline-friendly but can be easily fixed
Check_SRA_fastq_names <- function(namVec) { any(grepl("_[12]\\.fastq\\.gz$", namVec)) }
Fixup_SRA_fastq_names <- function(namVec) { sub('_([12])\\.fastq\\.gz$', '_R\\1.fastq.gz', namVec) }

## Abort if the FASTQ names are invalid -- but tolerate SRA names
## since we will convert them to valid when we make the symlinks.
ValidFastqNames(Fixup_SRA_fastq_names(fastqMap[,1]))


## Basic checks of the tsv, for paired FASTQs only.
if (grepl('_R[1,2]',fastqMap[1,1])) {
    ## Data seems to have paired FASTQs so make sure that they all come in pairs.
    x <- table(sub('_R[1,2]','',fastqMap[,1]))
    if (any(x!=2)) {
        cat("The following data sets seem to not have exactly one R1 and R2 file.\n")
        print( names(x)[x!=2] )
    }

    ## Make sure that for each pair of FASTQs the holding directory (final
    ## column) is the same.
    x <- paste(sub('_R[1,2]','',fastqMap[,1]), ':', fastqMap[,ncol(fastqMap)])
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
cat("Will make",length(fpaths),"symbolic links to FASTQs.\n")

## This converts relative paths to absolute or aborts. Safer so that the
## DADA2 pipeline can switch the current working directory and still find
## the FASTQs.
srcFqPaths <- normalizePath(fastqMap[,1], mustWork=TRUE)
pipeFriendlyFastqNames <- fpaths
if (Check_SRA_fastq_names(pipeFriendlyFastqNames)) {
    cat("Note:  Symbolic links will use \"_R1.fastq.gz\" rather than \"_1.fastq.gz\"\n",
        "      which appears in the source FASTQs, and similiar for the reverse reads.\n",
        "      Perhaps you are using data from the Sequencing Read Archive.  If so\n",
        "      then should verify that your read IDs used Illumina (CASAVA) format.  If\n",
        "      not, then you should set \"id.field\" in your pipeline parameters file.\n")
    pipeFriendlyFastqNames <- Fixup_SRA_fastq_names(pipeFriendlyFastqNames)
}
ok <- sum(file.symlink(srcFqPaths, pipeFriendlyFastqNames))
cat("Successfully created",ok,"symlinks.\n")
