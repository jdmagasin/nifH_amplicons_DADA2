#!/bin/bash

## Run the vsearch implementation of uchime3_denovo to check for chimera.

usageStr="

Usage:
    check_chimera_denovo.sh  asvAbundances.tsv  asv.fasta

Use vsearch implementation of uchime3_denovo to check for chimeric sequences.
Sequence ID's in asv.fasta must correspond to the row names in the abundance
table asvAbundances.tsv (which has columns = samples).

Output files:
    chimera.out      - Stats/results for all ASVs in uchimeout5 format.
    nonchimera.fasta - Non-chimeric ASVs.
    chimera.fasta    - Chimeric ASVs.

asvAbundances.tsv should be the abundance table from a single study because that
is what uchime3_denovo expects to see.  Although the script could execute on a
combined abundance table, it is not clear how count data from different studies
with different sequencing depths would impact uchime3_denovo.

Dependencies: vsearch must be installed and so must the R package ShortRead.
Conda packages exist for both, respecitvely 'vsearch' and
'bioconductor-shortread'.  You can install both into an environment for
processing ASVs after running the DADA2 nifH pipeline. E.g. my environment
'DADA2_nifH_postFilter' includes these packages as well as others needed for
check_nifH_contaminants.sh.

"

asvAbundances=$1
asvFasta=$2
if [ ! -f "$asvAbundances" ] || [ ! -f "$asvFasta" ] ; then
    echo
    echo "Need an abundance table such as from the DADA2 pipeline"
    echo "and a fasta file."
    echo "$usageStr"
    exit -1
fi

SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
if [ ! -d "$SDIR" ] ; then echo "Cannot get directory for $0"; exit -1; fi
"$SDIR/prepareFastaForUchime3_denonovo.R" "$asvAbundances" "$asvFasta" u3d.fasta

CMD="vsearch --uchime3_denovo u3d.fasta --nonchimeras nonchimera.fasta --chimera chimera.fasta --uchimeout chimera.out --uchimeout5"
echo "Will run vsearch as follows: "
echo "$CMD"
$($CMD)
howdItGo="$?"
rm u3d.fasta
if [ ! "$howdItGo" -eq 0 ] ; then
    echo "vsearch --uchime3_denovo returned error $howdItGo"
    exit -1
fi

## Strip off the size info that we had to add for uchime3_denovo.
for x in chimera.fasta nonchimera.fasta; do
    cat $x | sed 's/^>size=[0-9]*;/>/' > tmp.fasta
    mv tmp.fasta $x
done
echo "Done!"
exit 0
