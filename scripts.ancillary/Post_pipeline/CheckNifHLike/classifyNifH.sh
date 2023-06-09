#!/bin/bash

## Usage:
##    classifyNifH.sh fasta
## where fasta is a nucleotide fasta file of sequences you want to
## classify as nifH (e.g. amplicon ASVs).
## You must have NCBI blastx in your path!

NUCLFASTA=$1
if [ ! -f "$NUCLFASTA" ] ; then
    echo "Where is my nucleotide fasta of nifH-like sequences to classify?"
    exit -1
fi

## Make sure the positive and negative NifH databases exist alongside this script.
SDIR="$(dirname $0)"
DBDIR=`readlink -f $SDIR/DBs`
if [ ! -f "$DBDIR/NifH.pos.pdb" ] || [ ! -f "$DBDIR/NifH.neg.pdb" ] ; then
    echo "The positive and/or negative NifH example databases do not exist"
    echo "in directory DBs.  These are supposed to be included with this tool.\n"
    exit -1
fi
RHELP="$SDIR/scripts/classifyNifH.helper.R"
if [ ! -f "$RHELP" ] ; then
    echo "Missing $RHELP which should have been distributed with this tool."
    exit -1
fi

## This is how the positive and negative DBs are created:
## zcat ../Private/Data/positives.fasta.gz | makeblastdb -dbtype prot  -title NifH.pos -out NifH.pos
## zcat ../Private/Data/negatives.fasta | makeblastdb -dbtype prot  -title NifH.neg -out NifH.neg


## Make the databases be the same effective size so that e-values can be
## compared.  1E6 increases the pos DB slightly and doubles the neg DB, and the
## changes to evalues are, respectively, slight increases and doubling.
echo "Searching sequences in '$NUCLFASTA' against the positive NifH examples."
if [ ! -f posHits.tab ] ; then
    cat $NUCLFASTA \
      | blastx -db "$DBDIR/NifH.pos" -dbsize 1000000 -max_target_seqs=10 -out posHits.tab -outfmt 6 -num_threads 6
else
    echo "Already have posHits.tab. Reusing."
fi

echo "Searching sequences in '$NUCLFASTA' against the negative NifH examples."
if [ ! -f negHits.tab ] ; then
    cat $NUCLFASTA \
      | blastx -db "$DBDIR/NifH.neg" -dbsize 1000000 -max_target_seqs=10 -out negHits.tab -outfmt 6 -num_threads 6
else
    echo "Already have negHits.tab. Reusing."
fi

echo "Now classifying the sequences based on the blast results using the helper\n"
echo "script classifyNifH.helper.R\n"
$RHELP posHits.tab negHits.tab
echo "classifyNifH.sh has finished."
