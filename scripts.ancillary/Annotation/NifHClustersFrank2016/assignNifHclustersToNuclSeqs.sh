#!/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

##
## Usage:
##     assignNifHclustersToNuclSeqs.sh  fasta  [minNtLen]  [maxNtLen]
## 
## Assign NifH clusters to the nucleotide sequences in fasta, e.g. to nifH ASV's
## output by the DADA2 pipeline.
## - Use the CART decision tree method by Frank et al. 2016 (script downloaded
##   from lab web page and then modified for automated runs by jmagasin).
## - Since CART works on protein sequences, first we predict ORFs. Keep just the
##   NifH-like ORFs b/c we need to align them.  Also, optionally restrict to
##   to ORFs with nt length minNtLen to maxNtLen.
##
##
## Maintainer notes:
## 1. This script was originally developed for ASV's from DADA2, reflected in
##    some of the comments and variable names.  However, the script will work on
##    any nifH nucelotide sequences that capture (after ORF calling) residues
##    used by CART.  That is why the script name includes "NuclSeqs" rather than
##    "ASVs."
## 2. Reuse most of the ORF calling scripts that the pipeline uses to identify
##    NifH-like reads to be used in error models.  That's why there is a symlink
##    to NifH_prefilter.
##

asvFasta=$1
if [ ! -f "$asvFasta" ] ; then
    echo "Usage:"
    echo "    assignNifHclustersToNuclSeqs.sh  fasta  [minNtLen]  [maxNtLen]"
    echo
    echo "Pass a FASTA file with nifH nucleotide sequences.  These could be ASV's"
    echo "from DADA2, or similar.  Optionally pass integers for the minimum and"
    echo "number of nucleotides for open reading frames (ORFs).  Only ASV's with"
    echo "ORFs in the specified range will be used.  By default there is no length"
    echo "restriction."
    exit -1
fi

orfMinLen="$2"
orfMaxLen="$3"
if [ -z "$orfMinLen" ] ; then orfMinLen=1; fi
if [ -z "$orfMaxLen" ] ; then orfMaxLen=4294967295; fi
[ -z "${orfMinLen//[0-9]}" ] && [ -n "$orfMinLen" ] || exit -1
[ -z "${orfMaxLen//[0-9]}" ] && [ -n "$orfMaxLen" ] || exit -1

## This script calls several others that should be in the same directory.
SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"

## Define output files. As noted above, script works with any nifH nt seqs.
## Left the var names as "ASV" but updated the actual file names so they do not
## mention ASVs.
outbase=`echo $(basename ${asvFasta%.*})`
echo "Output files will be prefixed with $outbase"
echo
ORFS="${outbase}.orfs"
ORFSFAA="${ORFS}.faa"
ORFSGFF="${ORFS}.gff"
ORFSSIZEFILTFAA="${ORFS}_size_filtered.faa"
ASVSWITHNIFHFAA="${outbase}.withNifH.faa"
ASVSWITHNIFHTAB="${outbase}.hmmsearch.Fer_NifH.domtab"
ASVSWITHNIFHA2M="${outbase}.withNifH.a2m"
ASVSWITHNIFHALN="${outbase}.withNifH.aln"
ASVNIFHCLUSTERSFASTA="${outbase}.withNifH_Clusters.fasta"
ASV2CLUSTER="${outbase}.clusters.map"

if [ ! -f "$ORFSFAA" ] ; then
    echo "### Step 1.  Orf prediction. ###"
    ## Hard-coded parameters for FragGeneScan. Note the 1% error model selected
    ## for our DADA2 ASVs.
    run_FragGeneScan.pl -genome=$asvFasta -out="$ORFS" -complete=0 -train=illumina_10 -thread=6
fi

if [ ! -f "$ORFSSIZEFILTFAA" ] ; then
    echo "### Step 2. Filtering for orfs between $orfMinLen and $orfMaxLen nucleotides. ###"
    cat "$ORFSGFF" \
      | awk -v omin="$orfMinLen" -v omax="$orfMaxLen" -F"\t" '{ if ($5>=omin && $5<=omax) print $9 }' \
      | sed -e 's/;.*$//' -e 's/^ID=//' \
      > ids.tmp
    extractFasta.pl ids.tmp "$ORFSFAA" > "$ORFSSIZEFILTFAA"
    rm ids.tmp
fi

if [ ! -f "$ASVSWITHNIFHTAB" ] || [ ! -f "$ASVSWITHNIFHFAA" ] ; then
    echo "### Step 3.  Identifying orfs with NifH (PF00142) using trusted cut offs. ###"
    ## Identify ASV's that look like NifH, since those are the only ones we can align
    ## and pass to CART.  Note that we do not call findReadsWithNifHDomain.sh, which
    ## additionally uses bit score and length requirements.  Want more flexibility here.
    cat "$ORFSSIZEFILTFAA" | "$SDIR/NifH_prefilter/searchOrfsForPF00142.sh"
    mv hmmsearch.Fer_NifH.domtab "$ASVSWITHNIFHTAB"
    cat "$ASVSWITHNIFHTAB" | grep -v '^#' | tr -s ' ' '\t' | cut -f 1 \
      > ids.tmp
    extractFasta.pl ids.tmp "$ORFSFAA" > "$ASVSWITHNIFHFAA"
    rm ids.tmp
fi

if [ ! -f "$ASVSWITHNIFHA2M" ] ; then
    echo "### Step 4.  Multiple alignment. Adding Azotobacter vinelandii NifH WP_012698955.1 first. ###"
    echo -e ">A.vinelandii WP_012698955.1\n"\
            "STRLILHSKAQGTVMEMAASAGSVEDLELEDVLQIGFGGVKCVESGGPEPGVGCAGRGVITAINFLEEEGAYSDDLDFV"\
            "FYDVLGDVVCGGFAMPIRENKAQEIYIVCS\n"\
      | cat - "$ASVSWITHNIFHFAA" \
      > temp.faa
    ## We could use hmmalign, but in initial testing I seem to get the wrong
    ## cluster and/or subcluster assignments for reference UCYN-A nifH (A1,A2)
    ##     cat temp.faa | hmmalign --outformat A2M Fer4_NifH.hmm - > asvsWithNifH.a2m
    ## So let's try MAFFT...
    ## By default MAFFT output order is same as input, which we need. And it's fasta (A2M).
    ## There are lots of options to explore (mafft --man) but many seem to be for refining
    ## if there are at most ~200 sequences.
    mafft --thread 6 --clustalout temp.faa > "$ASVSWITHNIFHALN"
    mafft --thread 6 temp.faa > "$ASVSWITHNIFHA2M"
    rm temp.faa
fi

if [ ! -f "$ASVNIFHCLUSTERSFASTA" ] ; then
    echo "### Step 5.  NifH cluster and subcluster classification.  ###"
    ## Have the script figure out where the A.vinelandii residues begin.
    python "${SDIR}/NifH_Clusters.py" "$ASVSWITHNIFHA2M"
fi

if [ ! -f "$ASV2CLUSTER" ] ; then
    echo "### Step 6.  Making a map file of sequences to clusters, subclusters. ###"
    ## Example:  >A.vinelandii WP_012698955.1 main cluster = 1 subcluster = 1G
    grep '^>' $ASVNIFHCLUSTERSFASTA \
      | sed -e 's/^>\(..*\) main cluster = \(..*\) subcluster = \(..*\)$/\1,\2,\3/' \
            -e 's/_[0-9]*_[0-9]*_[-+]//' \
      | grep -v "WP_012698955" \
      > "$ASV2CLUSTER"
    echo "How many sequences are in each NifH cluster and subcluster:"
    cat "$ASV2CLUSTER" | cut -d, -f2,3 | sort | uniq -c | sort -k1gr
fi

echo "Done!"
exit 0
