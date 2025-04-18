#!/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

## Usage:
##    classifyNifH.sh fasta [positives fasta]
## where fasta is a nucleotide fasta file of sequences you want to
## classify as nifH (e.g. amplicon ASVs).
## You must have NCBI blastx in your path!
##
## Optional performance boost:  Pass 'positives fasta' that has
##   sequences (ASVs) that are known to be nifH.  Sequences in
##   'fasta' that exactly (full length) match the positives will
##   not be blastx'd.
##

NUCLFASTA=$1
if [ ! -f "$NUCLFASTA" ] ; then
    echo "Where is my nucleotide fasta of nifH-like sequences to classify?"
    exit -1
fi

POSFASTA=$2
if [ -n "$POSFASTA" ] && [ ! -f "$POSFASTA" ] ; then
    echo "Positives fasta file '$POSFASTA' does not exist."
    exit -1
fi
## Use link so that can swap in subset below. Caution!  Don't modify the
## original fasta via this link!
ln -f -s $NUCLFASTA ./blastxQueries.fasta

## Make sure the positive and negative NifH databases exist alongside this script.
SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
DBDIR=`readlink -f $SDIR/DBs`
if [ ! -f "$DBDIR/NifH.pos.pdb" ] || [ ! -f "$DBDIR/NifH.neg.pdb" ] ; then
    echo "The positive and/or negative NifH example databases do not exist"
    echo "in directory DBs.  These are supposed to be included with this tool.\n"
    exit -1
fi
RHELP="$SDIR/scripts/classifyNifH.helper.R"
RHELP2="$SDIR/scripts/classifyNifH.helper2.R"
if [ ! -f "$RHELP" ] || [ ! -f "$RHELP2" ] ; then
    echo "Missing $RHELP or $RHELP2 which should have been distributed with this tool."
    exit -1
fi

## If positives were provided, use them now.  Files created are positives_dontBlastx.ids
## and unsure_needBlastx.ids.  Either could be empty.
if [ -n "$POSFASTA" ] ; then
    rm -f positives_dontBlastx.ids unsure_needBlastx.ids
    $RHELP2 blastxQueries.fasta $POSFASTA
    if [ `cat positives_dontBlastx.ids | wc -l` -gt 0 ] ; then
        ## Some positives were found. Remove symlink and make a subset for blastx.
        rm blastxQueries.fasta
        extractFasta.pl unsure_needBlastx.ids $NUCLFASTA > blastxQueries.fasta
    fi
fi


## This is how the positive and negative DBs are created:
## zcat ../Private/Data/positives.fasta.gz | makeblastdb -dbtype prot  -title NifH.pos -out NifH.pos
## zcat ../Private/Data/negatives.fasta | makeblastdb -dbtype prot  -title NifH.neg -out NifH.neg


## Make the databases be the same effective size so that E-values can be
## compared.  1E6 increases the pos DB slightly and doubles the neg DB, and the
## changes to evalues are, respectively, slight increases and doubling.
echo "Searching sequences in '$NUCLFASTA' against the positive NifH examples."
if [ ! -f posHits.tab ] ; then
    cat blastxQueries.fasta \
      | blastx -db "$DBDIR/NifH.pos" -dbsize 1000000 -max_target_seqs=10 -out posHits.tab -outfmt 6 -num_threads 6
else
    echo "Already have posHits.tab. Reusing."
fi

echo "Searching sequences in '$NUCLFASTA' against the negative NifH examples."
if [ ! -f negHits.tab ] ; then
    cat blastxQueries.fasta \
      | blastx -db "$DBDIR/NifH.neg" -dbsize 1000000 -max_target_seqs=10 -out negHits.tab -outfmt 6 -num_threads 6
else
    echo "Already have negHits.tab. Reusing."
fi
## No longer need this (whether a symlink or file)
rm blastxQueries.fasta

echo "Now classifying the sequences based on the blast results using the helper\n"
echo "script classifyNifH.helper.R\n"
$RHELP posHits.tab negHits.tab
## Blend in ASVs determined to be nifH using the postives FASTA if provided.
if [ -n "$POSFASTA" ] ; then
    ## First .ids is from RHELP2 (and could be empty), and second .ids is from RHELP.
    cat positives_dontBlastx.ids positives.ids > posIds.tmp
    mv posIds.tmp positives.ids
    rm -f positives_dontBlastx.ids unsure_needBlastx.ids
fi
cat $NUCLFASTA | grep '^>' | sed 's/^>//' | cut -d' ' -f1 | sort | uniq > allseqs.tmp
cat allseqs.tmp positives.ids negatives.ids unsure.ids | sort | uniq -c \
  | grep '^ *1 ' | sed 's/^ *1 *//' > nohits.ids
rm allseqs.tmp
NH=`cat nohits.ids | wc -l`
echo "Identified $NH sequences with no hits to either database (in nohits.ids)."
if [ "$NH" -gt 0 ] ; then
    echo "Probably these sequences are not nifH."
fi
echo "classifyNifH.sh has finished."
