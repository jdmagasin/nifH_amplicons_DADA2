#!/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

##
## Blastn sequences (e.g. ASVs from nifH amplicon sequence data) against 2017
## ARB nifH sequences.
##

USAGE="blastnARB2017.sh nuclFasta outPath [queryPctId] [queryPctCvg]"

fasta=$1
if [ ! -f "$fasta" ] ; then
    echo "You did not specify a nucleotide FASTA file (nuclFasta)."
    echo -e "Usage:\n\t$USAGE\n"
    exit -1
fi
ODIR=$2
if [ -z "$ODIR" ] || [ -d "$ODIR" ] ; then
    echo "Please specify an output directory (outPath) that does not already exist."
    echo -e "Usage:\n\t$USAGE\n"
    exit -1
fi
pid=98
qcov=98
if [ ! -z "$3" ] ; then pid="$3" ; fi
if [ ! -z "$4" ] ; then qcov="$4"; fi

echo "Will blastn each sequence in $fasta against ARB 2017 nifH genes"
echo "and retain hits that are >${pid}%id and >${qcov}% query coverage."
echo "Alignment output will be limited to max 100 one line descriptions"
echo "and 10 alignments.  Outputs will be in $ODIR"

SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
## Pre-created.  See ARB_nifH_2017/readme.
NIFHDB="$SDIR/ARB_nifH_2017/NifHDB"

mkdir -p $ODIR

## Save as ASN.1 and then generate tabular and alignment formats.
##   https://www.ncbi.nlm.nih.gov/books/NBK569843/
ASN="${ODIR}/blastnNifH.${pid}id.${qcov}cvg.asn"
cat "$fasta" \
  | blastn -db "$NIFHDB"  -perc_identity "$pid"  -qcov_hsp_perc "$qcov"  -outfmt 11 \
  > "$ASN"

## Limit the number of alignments per query to 10.
blast_formatter -archive "$ASN" -outfmt 0 -num_alignments 10 -num_descriptions 100 \
  | gzip > "${ODIR}/blastnNifH.${pid}id.${qcov}cvg.aln.gz"

blast_formatter -archive "$ASN" -outfmt 6 > "${ODIR}/blastnNifH.${pid}id.${qcov}cvg.tab"

## The ASN file requires the DB to exist.  It might no longer (or have chagned)
## so delete the ASN.
rm "$ASN"

exit 0
