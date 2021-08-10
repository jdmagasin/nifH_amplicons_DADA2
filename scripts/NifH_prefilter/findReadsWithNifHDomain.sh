#!/bin/bash

## Search a FASTQ for reads that contain NifH domains and extract them to
## readsExtracted.fastq.gz

usage="Usage:
    findReadsWithNifHDomain.sh  fastq  outdir
The fastq may be gzip'd.  The outdir must not exist.
The script runs in the current directory but moves everything to outdir
when complete.
"

SDIR="$(dirname $0)"

FASTQGZ=$1
OUTDIR=$2
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
echo "Runnning on $FASTQGZ.  Results will be in $OUTDIR when script completes."
mkdir -p $OUTDIR

## Easier (but uglier) to run my scripts from the current directory and then to
## move resultws to OUTDIR.  Look for a few (not all) output files in the current
## directory (incomplete eariler result or parallel jobs, both problematic).
if [ -f "orfs.faa.gz" ] || [ -f readsWithNifH.ids.gz ] || [ -f readsExtracted.fastq.gz ] ; then
    echo "orfs.faa.gz,  readsWithNifH.ids.gz, or readsExtracted.fastq.gz already exists"
    echo "in the current directory. Aborting."
    exit -1
fi


echo -n "Predicting ORFs using FragGeneScan. Get some coffee..."
gunzip -t $FASTQGZ 2> /dev/null
if [ "$?" -eq 0 ] ; then
    cat $FASTQGZ | gunzip | $SDIR/fastq2orf.sh
else
    cat $FASTQGZ | $SDIR/fastq2orf.sh
fi
gzip orfs.faa
## Save some space since I only want the .faa.
rm orfs.{gff,ffn,out}
echo "done predicting ORFs!"
echo

## Now HMMER3. Then browse using the fields in the cut below and see the crappy
## stuff at the end (by bit score, and seem short).
echo -n "Searching ORFs for Pfam domain PF00142 (Fer4_NifH)..."
gunzip -c orfs.faa.gz | $SDIR/searchOrfsForPF00142.sh
gzip hmmsearch.Fer_NifH.domtab
echo "done searching for PF00142."
echo


## This magic retains the reads with bit score >150 and hmm coords taht span > 33 residues.
echo -n "Identifying reads that have Fer4_NifH with bit score > 150 and > 33 residues aligned..."
cat hmmsearch.Fer_NifH.domtab.gz | gunzip \
  | grep -v '^#' | tr -s " " "\t" | cut -f1,8,16-19 \
  | awk -F"\t" '{if ($2>150 && $4-$3+1 > 33) print $1}' | sed 's/_[^_]*_[^_]*_[+-]$//' \
  | gzip \
  > readsWithNifH.ids.gz
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
