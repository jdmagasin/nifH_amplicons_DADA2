#!/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

## Requires HMMER3 so you probably have already activated your DADA2_nifH
## environment.  Search specified ASV sequences (or any FASTA) for PF00142, the
## NifH domain.  This script is a simplified version of the scripts in
## NifH_prefilter.

ASVFASTA=$1
OUTDIR=$2

SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
SEARCHPF00142=`realpath "$SDIR/NifH_prefilter/searchOrfsForPF00142.sh"`
if [ ! -x "$SEARCHPF00142" ] ; then
    echo "Hmmm, I cannot find searchOrfsForPF00142.sh"
    exit -1
fi

if [ ! -f "$ASVFASTA" ] ; then
    echo "I need a FASTA file (uncompressed)."
    exit -1
fi
if [ -z "$OUTDIR" ] || [ -d "$OUTDIR" ] ; then
    echo "I need a nonexistent output directory name."
    exit -1
fi

mkdir $OUTDIR
if [ "$?" -ne 0 ] ; then
    echo "Failed to create $OUTDIR. It must not be hierarchical."
    exit -1
fi
CWD=`pwd`
cp "$ASVFASTA" "$OUTDIR/asvs.fasta"
cd "$OUTDIR"
ASVFASTA="asvs.fasta"

## To predict ORFs for the ASV's.  (Cannot use my fastq2orf.sh because I have a fasta.)
run_FragGeneScan.pl -genome="$ASVFASTA" -out=orfs -complete=0 -train=illumina_10 -thread=6 &> log.fgs.txt


## Now search for PF00142 using trusted cut offs (and not my >150 bits and >33
## residues post-HMMER filter).  338 ORFs pass.
cat orfs.faa | "$SEARCHPF00142"


## Get the ASVs for which we did, and did not, find NifH (whether the ASV had an
## ORF or not).  First get ids of ASVs with NifH.
grep -v '^#' hmmsearch.Fer_NifH.domtab | cut -d ' ' -f1 \
  | sed 's/_[0-9][0-9]*_[0-9][0-9]*_[-+]$//' \
  | uniq > ids.NifH.txt
extractFasta.pl ids.NifH.txt "$ASVFASTA" > asvsWithNifH.fasta
## Then infer those without.
cat asvs.fasta | grep '^>' | tr -d '>' \
  | cat - ids.NifH.txt | sort | uniq -c | grep '^  *1 ' \
  | sed 's/^  *1 //' > tmp
mv tmp ids.noNifH.txt
extractFasta.pl ids.noNifH.txt "$ASVFASTA" > asvsWithoutNifH.fasta


## Here are the most abundant ASVs that seem to have, or not have, NifH.
echo
echo "Here are the top 10 ASV's that have NifH domains:"
grep '^>' asvsWithNifH.fasta | sed 's/>ASV\.//' | sort -n | head -n 10

echo
echo "Here are the top 10 ASV's that do not have NifH domains:"
grep '^>' asvsWithoutNifH.fasta | sed 's/>ASV\.//' | sort -n | head -n 10

cd "$CWD"

echo "Done!"
exit 0
