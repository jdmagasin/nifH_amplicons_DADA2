#!/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

## Search a FASTQ for reads that contain NifH domains and extract them to
## readsExtracted.fastq.gz

usage="Usage:
    findReadsWithNifHDomain.sh  fastq  outdir  [minLen] [minBitScore]
The fastq may be gzip'd.  The outdir must not exist.  The script runs in the
current directory but moves everything to outdir when complete.  Reads that
align to the NifH PFAM model PF00142 at trusted cutoffs (encoded in PF00142) are
retained by HMMER and post-filtered for a minLen residues (default >33) and
minBitScore (default >150).
"

SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
if [ ! -d "$SDIR" ] ; then echo "Cannot determine directory of $0"; exit -1 ; fi

FASTQGZ=$1
OUTDIR=$2
MINLEN=$3
MINBITS=$4
if [ ! -f "$FASTQGZ" ] ; then
    echo "Fastq does not exist."
    echo "$usage"
    exit -1
fi
if [ -d "$OUTDIR" ] ; then
    echo "Outdir exists already."
    echo "$usage"
    exit -1
fi
if [ -z "$MINLEN" ] ; then
    MINLEN=33
fi
if [ -z "$MINBITS" ] ; then
    MINBITS=150
fi
echo "Running on $FASTQGZ.  Results will be in $OUTDIR when script completes."
mkdir -p $OUTDIR

## Easier (but uglier) to run my scripts from the current directory and then to
## move results to OUTDIR.  Look for a few (not all) output files in the current
## directory (incomplete eariler result or parallel jobs, both problematic).
for oldfile in "orfs.faa.gz" "readsWithNifH.ids.gz" "readsExtracted.fastq.gz" "hmmsearch.Fer_NifH.domtab.gz" ; do
    if [ -f "$oldfile" ] ; then
        echo "Deleting old $oldfile"
        rm "$oldfile"
    fi
done


echo "Predicting ORFs using FragGeneScan. Get some coffee..."
gunzip -t $FASTQGZ 2> /dev/null
if [ "$?" -eq 0 ] ; then
    cat $FASTQGZ | gunzip | $SDIR/fastq2orf.sh
else
    cat $FASTQGZ | $SDIR/fastq2orf.sh
fi
gzip orfs.faa
## Save some space since I only want the .faa.
rm orfs.{gff,ffn,out}
echo "Done predicting ORFs!"
echo

## Now HMMER3. Then browse using the fields in the cut below and see the crappy
## stuff at the end (by bit score, and seem short).
echo -n "Searching ORFs for Pfam domain PF00142 (Fer4_NifH)..."
gunzip -c orfs.faa.gz | $SDIR/searchOrfsForPF00142.sh
gzip hmmsearch.Fer_NifH.domtab
echo "done searching for PF00142."
echo


## This magic retains the reads with bit score >MINBITS and hmm coords that
## span >MINLEN residues.
echo -n "Identifying reads that have Fer4_NifH with bit score > ${MINBITS} and > ${MINLEN} residues aligned..."
cat hmmsearch.Fer_NifH.domtab.gz | gunzip \
  | grep -v '^#' | tr -s " " "\t" | cut -f1,8,16-19 \
  | awk -F"\t" -v minBits="$MINBITS" -v minLen="$MINLEN" \
        '{if ($2 > minBits  &&  $4-$3+1 > minLen) print $1}' \
  | sed 's/_[^_]*_[^_]*_[+-]$//' \
  | sort | uniq \
  | gzip \
  > readsWithNifH.ids.gz
numNifHReads=`cat readsWithNifH.ids.gz | gunzip | wc -l`
if [ "$numNifHReads" -eq 0 ] ; then
    echo
    echo "PROBLEM!  No nifH-like reads were found by the NifH prefilter for:"
    echo "    $FASTQGZ"
    echo "This will cause a later failure when the pipeline has no reads for building error models."
    echo "The usual reason for not finding any nifH-like reads is that your NifH_minBits is too high"
    echo "for the quality of your sequencing data.  Yours is ${MINBITS}.  Try lowering NifH_minBits."
    echo "Set it to 0 if you want to use the trusted bit cut off defined in PF00142."
    exit -1
fi
echo "done."
echo


## Then probably easiest to use R ShortRead package to extract the good reads.
## Then DADA2 to check out the error models.
echo -n "Running R script to extract the NifH-containing reads into a fastq..."
Rscript $SDIR/extractReadsFromFastq.R  readsWithNifH.ids.gz  $FASTQGZ
echo "done."
echo

echo "Moving everything to the $OUTDIR"
mv orfs* hmmsearch.Fer_NifH.* readsWithNifH.* readsExtracted.fastq.gz $OUTDIR
echo "Predicted" `gunzip -c $OUTDIR/orfs.faa.gz | grep -c '^>'` "orfs."
echo "Identified" `gunzip -c $OUTDIR/readsWithNifH.ids.gz* | wc -l` "reads with NifH domains."

exit 0
