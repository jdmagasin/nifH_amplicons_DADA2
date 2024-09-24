#!/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

##
## Runs findReadsWithNifHDomain.sh on the cutadapt-trimmed fastqs.  Just run on
## the R1 reads.  That should be sufficient, and findReadsWithNifHDomain.sh (and
## the scripts it calls) used fixed output file names, so we cannot have R1 and
## R2 results in the same directory.
##

DATATRIMMEDDIR=$1
if [ ! -d "$DATATRIMMEDDIR" ] ; then
    echo "Pass the path to the root directory that holds the trimmed fastq files for R1."
    echo "Optionally also pass in the min number of residues and then the min bit score"
    echo "at which predicted ORFs must align to NifH."
    exit -1
fi

MINLEN=$2
MINBITS=$3

SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
if [ ! -d "$SDIR" ] ; then echo "Cannot determine directory of $0"; exit -1 ; fi

for FQ in `find -L "$DATATRIMMEDDIR" -name '*_R1*.fastq.gz'`; do
    OUTDIR=`echo $(dirname $FQ) | sed -e "s:$DATATRIMMEDDIR::" -e 's:\/$::'`
    OUTDIR="Data.NifH_prefilter/$OUTDIR"
    if [ ! -d "$OUTDIR" ] ; then
        echo "Searching for NifH in reads in $FQ..."
          ## ~Ugly but the script outputs are to the current directory
          ## and moved at completion. Save the log to cwd too (since the
          ## script expects the outdir not to exist).
        $SDIR/findReadsWithNifHDomain.sh  "$FQ"  "$OUTDIR" "$MINLEN" "$MINBITS" \
            > log.nifScan.txt 2>&1
        if [ "$?" -ne 0 ] ; then
            echo "findReadsWithNifHDomain.sh exited abnormally."
            exit -1
        fi
        mv log.nifScan.txt "$OUTDIR"
    else
        echo "Already have $OUTDIR"
    fi
done
exit 0
