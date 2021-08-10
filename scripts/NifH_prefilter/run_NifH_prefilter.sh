#!/bin/bash

##
## Runs findReadsWithNifHDomain.sh on the cutadapt-trimmed fastqs.  Just run on
## the R1 reads.  That should be sufficient, and findReadsWithNifHDomain.sh (and
## the scripts it calls) used fixed output file names, so we cannot have R1 and
## R2 results in the same directory.
##

DATATRIMMEDDIR=$1
if [ ! -d "$DATATRIMMEDDIR" ] ; then
    echo "Pass the path to the root directory that holds the trimmed fastq files for R1."
    exit -1
fi

## Expect findReadsWithNifHDomain.sh to be in the same directory as this script.
## From https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel
SDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

for FQ in `find -L "$DATATRIMMEDDIR" -name '*_R1_*.fastq.gz'`; do
    OUTDIR=`echo $(dirname $FQ) | sed -e "s:$DATATRIMMEDDIR::" -e 's:\/$::'`
    OUTDIR="Data.NifH_prefilter/$OUTDIR"
    if [ ! -d "$OUTDIR" ] ; then
        echo "Searching for NifH in reads in $FQ..."
          ## ~Ugly but the script outputs are to the current directory
          ## and moved at completion. Save the log to cwd too (since the
          ## script expects the outdir not to exist).
        $SDIR/findReadsWithNifHDomain.sh  "$FQ"  "$OUTDIR" \
            > log.nifScan.txt 2>&1
        mv log.nifScan.txt "$OUTDIR"
    else
        echo "Already have $OUTDIR"
    fi
done
