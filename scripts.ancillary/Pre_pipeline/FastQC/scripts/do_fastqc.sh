#!/bin/bash

##
## Helper to run FastQC.

## Currently uses the default fastqc parameters in particluar for
## --contaminants, --adapters, and --limits.  Possibly we could ammend
## the default lists for each of these to include lab/nifH-specific info.
##   However, my check_nifH_contaminants.sh will check ASVs after running the
##   pipeline.  This has a possible small advantage of letting nifH-like
##   contaminants impact DADA2 error models.  For the other
##
## Gotcha: If you pass an empty fastqs.tmp to fastqc, then fastqc will try to
## launch in GUI mode, which will fail if you did not log in with X11 forwarding
## enabled.  This script insists on non-GUI-mode:  do not launch fastqc if
## fastqs.tmp is empty.
##

OUTDIR="FastQC.out"

usageStr="

Usage:
    do_fastqc.sh topDir

Run FastQC on all FASTQ files (with extension .fastq.gz) that can be found
recursively within topDir.  The output directory is ${OUTDIR}.  Only the HTML
quality reports are retained, not the zip files.

"


topDir=$1
if [ ! -d "$topDir" ] ; then
    echo
    echo "Need a directory within which to find fastq files."
    echo "$usageStr"
    exit -1
fi
if [ -d "$OUTDIR" ] ; then
    echo "$OUTDIR already exists. Please rename or remove it."
    exit -1
fi


find "$topDir" -name "*.fastq.gz" > fastqs.tmp
if [ ! -s fastqs.tmp ] ; then
    echo "Found no fastq.gz files."
    rm fastqs.tmp
    exit -1
fi

mkdir "$OUTDIR"
fastqc -o "$OUTDIR" --noextract -t 5 -f fastq `cat fastqs.tmp`
rm fastqs.tmp ${OUTDIR}/*.zip
