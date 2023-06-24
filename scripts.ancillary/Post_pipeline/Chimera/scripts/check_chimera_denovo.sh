#!/bin/bash

## Run the vsearch implementation of uchime3_denovo to check for chimera.

WORKDIR=Work_on_chimera

usageStr="

Usage:
    check_chimera_denovo.sh  asvAbundances.tsv  asv.fasta

Use vsearch implementation of uchime3_denovo to check for chimeric sequences.
Sequence ID's in asv.fasta must correspond to the row names in the abundance
table asvAbundances.tsv (which has columns = samples).  uchime3_denovo is run
separately for each sample in the table.  Then each of the samples in which
an ASV was detected votes on whether the ASV is a chimera (majority vote wins).

Output:
    nonchimera.fasta - Non-chimeric ASVs.
    chimera.fasta    - Chimeric ASVs.
    ${WORKDIR}  - Directory with uchime3_denovo outputs for each sample.
                       Includes uchime3_denovo stats (uchimeout5 format) and
                       chimeraVotes.csv which shows the number of samples
                       that voted for an ASV to be chimeric or not.  Only
                       ASVs with >0 chimera votes are included.

Dependencies: vsearch must be installed and so must the R package ShortRead.
vsearch can be installed as described in INSTALL_ancillary.txt.  ShortRead is
required for dada2 and gets installed as described in INSTALL.txt.

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

if [ -d "$WORKDIR" ] ; then echo "$WORKDIR already exists."; exit -1; fi

## Ensure that we have required scripts.
SDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
if [ ! -d "$SDIR" ] ; then echo "Cannot get directory for $0"; exit -1; fi
PREPAREUCHIME3="$SDIR/prepareFastaForUchime3_denonovo.R"
DECIDECHIMERA="$SDIR/decideChimeraHelper.R"
if [ ! -x "$PREPAREUCHIME3" ] || [ ! -x "$DECIDECHIMERA" ] ; then
    echo "Missing required helper scripts for check_chimera_denovo.sh.  Aborting"
    exit -1
fi
if [ -z `which vsearch` ] ; then
    echo "check_chimera_denovo.sh cannot find vsearch.  Aborting"
    exit -1
fi

## Prepare one FASTA for each column of the abundance table.
"$PREPAREUCHIME3" "$asvAbundances" "$asvFasta" "$WORKDIR"
if [ ! "$?" -eq 0 ] ; then
    echo "prepareFastaForUchime3_denonovo.R failed. Aborting."
    exit -1
fi

ccdPWD=$(realpath `pwd`)
cd $WORKDIR

echo "Will run vsearch on each sample individually.  Outputs are in $WORKDIR"
for fas in `find . -name "*.fasta"`; do
    desc="${fas%.*}"
    ## The 'chimera' and 'nonchimera' fastas are used by DECIDECHIMERA.
    CMD="vsearch --uchime3_denovo ${fas} --nonchimeras ${desc}.nonchimera.fasta --chimera ${desc}.chimera.fasta --uchimeout ${desc}.chimera.out --uchimeout5"
    echo "### Running uchime3 denovo on $fas using the following command:"
    echo "    $CMD"
    $($CMD)
    howdItGo="$?"
    if [ ! "$howdItGo" -eq 0 ] ; then
        echo "vsearch --uchime3_denovo returned error $howdItGo"
        exit -1
    fi
done

## Strip off the ASV "size" (counts) info that was added for uchime3_denovo so
## that the ids files output by DECIDECHIMERA will have IDs that extractFasta.pl
## can find in asvFasta.
for x in *.fasta; do
    cat $x | sed 's/^>size=[0-9]*;/>/' > tmp.fasta
    mv tmp.fasta $x
done

## Vote on chimera
"$DECIDECHIMERA"
if [ ! "$?" -eq 0 ] ; then
    echo "decideChimeraHelper.R failed. Aborting."
    exit -1
fi

## Wrap-up:  Delete the fastas (the per sample and the vsearch outputs) and
## pop up to initial directory to create the main output files.
rm *.fasta
cd $ccdPWD

## Create main output files.
extractFasta.pl "${WORKDIR}/chimera.ids"    "$asvFasta" > chimera.fasta
extractFasta.pl "${WORKDIR}/nonchimera.ids" "$asvFasta" > nonchimera.fasta

echo "Done!"
exit 0
