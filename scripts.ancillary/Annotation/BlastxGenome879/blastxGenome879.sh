#!/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

##
## Blastx sequences (e.g. ASVs from nifH amplicon sequence data) against genome879
## taken from the Zehr Lab pages (March 2019 version).
##

## blastx has no -perc_identity hit filter so we do ourselves with awk.

USAGE="blastxGenome879.sh nuclFasta outPath [queryPctId] [queryPctCvg]"

SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
## Pre-created.  
NIFHDB="$SDIR/Genome879DB/genome879DB"

fasta=$1
if [ ! -f "$fasta" ] ; then
    echo "You did not specify a nucleotide FASTA file (nuclFasta)."
    echo -e "Usage:\n\t$USAGE\n"
    exit -1
fi
ODIR=$2
if [ -d "$ODIR" ] ; then
    echo "Please specify an output directory (outPath) that does not already exist."
    echo -e "Usage:\n\t$USAGE\n"
    exit -1
fi
pid=92
qcov=98
if [ ! -z "$3" ] ; then pid="$3" ; fi
if [ ! -z "$4" ] ; then qcov="$4"; fi

echo "Will blastx each sequence in $fasta against Genome879 NifH proteins"
echo "and retain hits that are >${pid}% amino identical >${qcov}% query coverage,"
echo "The % amino identical test only applies to the tabular output.  Alignment"
echo "output will be limited to max 100 one line descriptions and 10 alignments."
echo "Outputs will be in $ODIR"

mkdir -p $ODIR

## Save as ASN.1 and then generate tabular and alignment formats.
##   https://www.ncbi.nlm.nih.gov/books/NBK569843/
## The above page suggests using max_target_seqs "control the number of matches
## recorded in the alignment."  However, there are issues with max_target_seqs
## namely it (and the num_ options) can be used *during* the search and could
## miss top hits:
##   https://blastedbio.blogspot.com/2018/11/blast-max-alignment-limits-repartee-one.html
## max_target_seqs has a default of 500 so its use is unavoidable (I think).
## To keep the asn file from blowing up (since it was multi-GB under the default),
## set max_target_seqs to ... 200.
ASN="${ODIR}/blastxGenome879.${qcov}cvg.asn"
cat "$fasta" | blastx -db "$NIFHDB" -qcov_hsp_perc "$qcov" -outfmt 11 -max_target_seqs 200 > "$ASN"

## Limit the number of alignments per query to 10 (rather than default 250) and
## then number of descriptions to 100 (rather than default 500).  Note (above)
## that *search* results (and search itself) were restricted using
## max_target_seqs.  I think that here the num_ args can only be w.r.t.  what is
## shown. (Vs. impacting the actual search -- see blog post above.)
blast_formatter -archive "$ASN" -outfmt 0 -num_alignments 100 -num_descriptions 100 \
  | gzip > "${ODIR}/blastxGenome879.${qcov}cvg.aln.gz"

blast_formatter -archive "$ASN" -outfmt 6 \
  | awk -F"\t" -v pid="$pid" '$3 > pid' \
  > "${ODIR}/blastxGenome879.${pid}pid.${qcov}cvg.tab"               

## The ASN file requires the DB to exist.  It might no longer (or have chagned)
## so delete the ASN.
rm "$ASN"

exit 0
