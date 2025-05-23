#!/bin/bash

## Copyright (C) 2023 Jonathan D. Magasin

##
## Search nifH ASVs for possible contaminants that are:
##   - known from previous work (pulbications) to exist in PCR reagents
##     or in primers as received from the vendors.
##   - observed in the current studies negative controls.
##

usageStr="

Usage:
    check_nifH_contaminants  asvs.fasta  pctid  [myContaminants.fasta]

Check a nifH ASV nucleotide fasta file for known contaminants in nifH amplicon
studies as well as for user/lab-specific contaminants identified in negative
controls.

Here are the steps and the [key output files]:
 - Translate the ASVs using FragGeneScan. [asvs.faa]

 - Create a DB of known contaminants (below) and those provided in
   myContaminants.fasta (optionally).  [nifHContamDB]

 - Search the DB using BLASTp output for hits at > the specified amino identity cut
   off [blast.out and blast.pctId.out]

The known contaminants include 26 nifH-like sequences known to haunt reagents /
primers from the following studies:
   - Zehr et a. 2003          [3 sequences]
   - Goto et al. 2005         [6 sequences]
   - Farnelid et al. 2009     [11 sequences, fewer than reported in the pub]
   - Turk-Kubo et al. 2011    [6 confirmed contaminants in Fig. 1]

Dependencies: This script requires FragGeneScan and the BLAST suite.
FragGeneScan should already be available as part of installing the
DADA2 pipeline.  If your system does not already have NCBI's BLAST+
tools, you can install them as described in INSTALL_ancillary.txt.


"

asvFasta=$1
aminoId=$2
myContam=$3
if [ ! -f "$asvFasta" ]; then
    echo "Need ASV nt sequences in fasta format."
    echo "$usageStr"
    exit -1
fi
if [ -z "$aminoId" ] ; then
    echo "Need % amino identity cut off in rangge [0,100)."
    echo "$usageStr"
    exit -1
fi


##
## Create contaminants.fasta that contains the user's contaminants followed by
## known.  The latter are included in this script as a heredoc because it avoids
## depending on an external FASTA.
##

CONTAMFASTA=contaminants_tmp.fasta
if [ ! -z "$myContam" ] ; then
    if [ ! -f "$myContam" ] ; then
        echo "I cannot find your contaminants fasta file."
        echo "$usageStr"
        exit -1
    fi
    cp "$myContam" "$CONTAMFASTA"
else
    echo '' > "$CONTAMFASTA"
fi

cat << _known_contaminants_ >> "$CONTAMFASTA"
>zehr2003_AY225105.1 Uncultured bacterium 400A2 nitrogenase Fe protein (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCCAAGGCCCAGGACACCATCCTGTCGCTGGCGGCCGAAGCCGGCTCGG
TGGAGGACCTGGAGCTCGAGGACGTGATGAAGATCGGCTACGAGGACATCCGTTGCGTCGAATCCGGTGG
CCCGGAGCCCGGAATGGGCTGCGCCGGCCGCGGCGTGATCACCTCGATCAACTTCCTTGAAGAAAACGGC
GCCTACGACGGTGTTGACTACGTCTCTTACGACGTGCTGGGCGACGTGGTGTGCGGCGGCTTCGCCATGC
CCATCCGCGAAAACAAGGCGCAAGAGATCTACATCGTCATGTCC

>zehr2003_AY225106.1 Uncultured bacterium 436A5-2 nitrogenase Fe protein (nifH) gene, partial cds
TCGACCCGCCTGATCTTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGAGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>zehr2003_AY225107.1 Uncultured bacterium 2580A2 nitrogenase Fe protein (nifH) gene, partial cds
TCGACCCGTCTTATCCTTCACTCGAAGGCCCAGGACACCATCCTCAGCCTCGCCGCTGCCGCTGGTTCCG
TTGAGGACCTCGAAATCGAAGACGTCATGAAGGTCGGTTACCTCGACATCCGTTGCGTCGAGTCGGGTGG
TCCGGAGCCGGGCGTTGGCTGCGCGGGTCGTGGTGTTATCACCTCGATCAACTTCCTCGAGGAAAACGGC
GCTTACGAAGACGTTGACTACGTTTCCTACGACGTTCTCGGCGACGTGGTCTGCGGCGGTTTCGCCATGC
CGATCCGTGAGAACAAGGCTCAGGAAATTTACATCGTCATGTCC


>goto2005_AB198366.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conI1-1
CTCCACGCGCCTCATCCTGCACGCCAAGGCTCAGGACACCATCCTCAGCCTCGCCGCCGAGCAGGGCAGC
GTCGAGGACCTCGAACTCGAAGACGTAATGAAGATCGGCTACCAAAACATCCGTTGTGTGGAATCCGGCG
GTCCGGAGCCGGGCGTCGGCTGCGCTGGCCGCGGCGTCATCACCTCGATCAACTTCCTCGAGGAAAACGG
CGCCTACGAGGACATCGACTACGTCTCCTACGACGTGCTGGGC

>goto2005_AB198390.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conT3-1
CTCCACCCGCCTGATCCTGCACGCAAAGGCTCAGGACACCATCCTGTCGCTGGCCGCTGAAGCCGGTTCG
GTGGAAGACCTCGAGATCGATGATGTGATGAAGGTGGGCTATCGCGACATCCGTTGCGTGGAGTCCGGTG
GTCCTGAGCCCGGCGTGGGCTGTGCCGGCCGTGGCGTGATCACCTCGATCAACTTCCTGGAAGAAAACGG
TGCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGC

>goto2005_AB198383.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conT1-1
CTCCACCCGCCTGATCCTGCACGCTAAAGCACAGAACACCATTATGGAGATGGCCGCGGAAGTCGGCTCG
GTCGAGGATCTCGAGCTCGAAGACGTGCTGCAAATTGGCTACGGCGACGTGCGCTGCGCGGAATCCGGCG
GCCCGGAGCCAGGCGTCGGCTGCGCGGGGCGCGGCGTGATCACGGCGATCAACTTTCTTGAAGAAGAAGG
CGCCTACGAGGACGATCTCGATTTTGTCTTCTATGACGTGCTCGGC

>goto2005_AB198374.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conI3-1
CTCCACCCGTCTGATCCTTCACGCTAAAGCCCAGAACACCATCATGGAGATGGCGGCGGAAGTGGGCTCG
GTCGAGGATCTGGAGCTCGAAGACGTTCTGCAAATCGGCTATGGCGATGTCCGTTGCGCCGAATCCGGCG
GCCCGGAGCCAGGCGTCGGCTGCGCCGGACGCGGGGTGATCACCGCCATCAACTTCCTCGAGGAAGAAGG
CGCCTATGAAGAAGATTTGGATTTCGTCTTCTATGACGTCCTCGGC

>goto2005_AB198376.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conI4-1
TTCTACAAGGCTCCTGCTTGGCGGCCTTGCCCAAAAAACAGTATTAGACACACTTCGCGAAGAAGGAGAA
GATGTAGAATTAGACGATGTTGTAAAACTCGGCTTCATGGGTTCACGTTGTGTAGAATCGGGTGGGCCTG
AACCCGGTGTTGGTTGTGCAGGCCGTGGTATCATTACTTCCATTAACCTGCTTGAACAATTAGGTGCTTA
TAACGAAGAACAATTTGAAACAGATTTTGTGTTTTATGATGTTTTGGGT

>goto2005_AB198387.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conT2-1
TTCCACGAAGAACCTAATGGAGGGCAGGCAGATTCCCACGGTTTTAGACCAGATTTTGAACAAGGGAAAG
GACTTGCGTCTAGAAGATGTGGTCTTTCCCAGCAAGACCGGGGTTCTGTGTGTGGAGGCCGGAGGCCCTG
CCCCGGGTGTCGGTTGTGCAGGAAGAGAGATTATTTCCGCCTTTGAGAAGCTGGAAGAATTGAAAGCCTT
TGAGGCCTACAAGCCGGACATTGTCATCTATGACGTATTGGGG


>farnelid2009_EU916690.1 Uncultured bacterium clone balntcK49 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGAGGACCTCGAACTGGTAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>farnelid2009_EU916675.1 Uncultured bacterium clone balntcK3 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGGAGCTGGTTCGG
TGGAGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>farnelid2009_EU916666.1 Uncultured bacterium clone balntcK15 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGATGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGGGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>farnelid2009_EU916694.1 Uncultured bacterium clone balntcK54 nitorgenase reductase (nifH) gene, partial cds
GCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGGGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>farnelid2009_EU916686.1 Uncultured bacterium clone balntcK41 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGGGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>farnelid2009_EU916671.1 Uncultured bacterium clone balntcK21 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGAGGACCTCGAACTGGAAGACGTGATGAGGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGTAATCTACATCGTCATGTCC

>farnelid2009_EU916668.1 Uncultured bacterium clone balntcK18 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGAGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGTGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCGCTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>farnelid2009_EU916368.1 Uncultured bacterium clone bal19mayC48 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCCAAGGCACAGGACACCATCCTGTCGCTGGCCGCTGAAGCCGGTTCGG
TGGAAGACCTCGAGCTCGAGGACGTCATGAAGATCGGCTACAGGGACATCCGCTGCGTCGAATCCGGTGG
CCCGGAACCCGGAGTGGGTTGTGCCGGCCGTGGCGTCATCACCTCGATCAACTTCCTCGAGGAAAACGGT
GCCTATGACGGCGTCGACTACGTCAGCTACGACGTGCTGGGCGACGTGGTGTGCGGCGGCTTTGCCATGC
CGATCCGCGAAAACAAGGCGCAGGAAATCTACATCGTGATGTCG

>farnelid2009_EU916480.1 Uncultured bacterium clone bal30julG27 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGAGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>farnelid2009_EU916514.1 Uncultured bacterium clone bal30julG60 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGAGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGCGGCTTTGCCATGC
CGATCCGCGAAAACAAGGCGCAGGAAATCTACATCGTGATGTCG

>farnelid2009_EU916706.1 Uncultured bacterium clone balntcK9 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCGCTGAAGCTGGTTCGG
TGGAGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC


>turkkubo2011_AB198382.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conI6-2
CTCCACCCGTCTGATCCTTCACGCTAAAGCCCAGAGCACCATCATGGAGATGGCGGCGGAAGTGGGCTCG
GTCGAGGATCTGGAGCTCGAAGACGTTCTGCAAATCGGCTATGGCGATGTCCGTTGCGCCGAATCCGGCG
GCCCGGAGCCAGGCGTCGGCTGCGCCGGACGCGGGGTGATCACCGCCATCAACTTCCTCGAGGAAGAAGG
CGCCTATGAAGAAGATTTGGATTTCGTCTTCTATGACGTCCTCGGC

>turkkubo2011_AB198384.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conT1-2
CTCCACCCGCCTGATCCTGCACGCTAAAGCACAGAACACCATTATGGAGATGGCCGCGGAAGTCGGCTCG
GTCGAGGATCTCGAGCTCGAAGACGTGCTGCAAATTGGCTACGGCGACGTGCGCTGCGCGGAATCCGGCG
GCCCGGAGCCAGGCGTCGGCTGCGCGGGGCGCGGCGTGATCACGGCGATCAACTTTCTTGAAGAAGAAGG
CGCCTACGAGGACGATCTCGATTTTGTCTTCTATGACGTGCTCGGC

##>turkkubo2011_EU916368.1 Uncultured bacterium clone bal19mayC48 nitorgenase reductase (nifH) gene, partial cds
##TCGACCCGCCTGATCCTGCACGCCAAGGCACAGGACACCATCCTGTCGCTGGCCGCTGAAGCCGGTTCGG
##TGGAAGACCTCGAGCTCGAGGACGTCATGAAGATCGGCTACAGGGACATCCGCTGCGTCGAATCCGGTGG
##CCCGGAACCCGGAGTGGGTTGTGCCGGCCGTGGCGTCATCACCTCGATCAACTTCCTCGAGGAAAACGGT
##GCCTATGACGGCGTCGACTACGTCAGCTACGACGTGCTGGGCGACGTGGTGTGCGGCGGCTTTGCCATGC
##CGATCCGCGAAAACAAGGCGCAGGAAATCTACATCGTGATGTCG

>turkkubo2011_AB198391.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conT3-2
CTCCACCCGCCTGATCCTGCACGCAAAGGCTCAGGACACCATCCTGTCGCTGGCCGCTGAAGCCGGTTCG
GTGGAAGACCTCGAGATCGATGATGTGATGAAGGTGGGCTATCGCGACATCCGTTGCGTGGAGTCCGGTG
GTCCTGAGCCCGGCGTGGGCTGTGCCGGCCGTGGCGTGATCACCTCGATCAACTTCCTGGAAGAAAACGG
TGCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGC

>turkkubo2011_EU916669.1 Uncultured bacterium clone balntcK19 nitorgenase reductase (nifH) gene, partial cds
TCGACCCGCCTGATCCTGCACGCGAAGGCTCAGGACACCATCTTGTCGCTGGCCACTGAAGCCGGTTCGG
TGGAGGACCTCGAACTGGAAGACGTGATGAAGGTCGGGTACCGCGACATCCGTTGCGTGGAATCCGGCGG
CCCTGAGCCTGGGGTTGGCTGCGCCGGCCGCGGCGTGATCACTTCGATCAACTTCCTGGAAGAAAACGGC
GCCTACGAAGGCGTGGACTATGTGTCCTACGACGTGCTGGGCGACGTGGTGTGCGGTGGCTTTGCCATGC
CCATCCGTGAGAACAAGGCACAGGAAATCTACATCGTCATGTCC

>turkkubo2011_AB198373.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conI2-4
CTCCACGCGCCTCATCCTGCACGCCAAGGCTCAGGACACCATCCTCAGCCTCGCCGCCGAGCAGGGCAGC
GTCGAGGACCTCGAACTCGAAGACGTAATGAAGATCGGCTACCAAAACATCCGTTGTGTGGAATCCGGCG
GTCCGGAGCCGGGCGTCGGCTGCGCTGGCCGCGGCGTCATCACCTCGATCAACTTCCTCGAGGAAAACGG
CGCCTACGAGGACATCGACTACGTCTCCTACGACGTGCTGGGC

##>turkkubo2011_AY225107.1 Uncultured bacterium 2580A2 nitrogenase Fe protein (nifH) gene, partial cds
##TCGACCCGTCTTATCCTTCACTCGAAGGCCCAGGACACCATCCTCAGCCTCGCCGCTGCCGCTGGTTCCG
##TTGAGGACCTCGAAATCGAAGACGTCATGAAGGTCGGTTACCTCGACATCCGTTGCGTCGAGTCGGGTGG
##TCCGGAGCCGGGCGTTGGCTGCGCGGGTCGTGGTGTTATCACCTCGATCAACTTCCTCGAGGAAAACGGC
##GCTTACGAAGACGTTGACTACGTTTCCTACGACGTTCTCGGCGACGTGGTCTGCGGCGGTTTCGCCATGC
##CGATCCGTGAGAACAAGGCTCAGGAAATTTACATCGTCATGTCC

>turkkubo2011_AB198377.1 Uncultured bacterium nifH gene for nitrogenase Fe protein, partial cds, clone: conI4-2
TTCTACAAGGCTCCTGCTTGGCGGCCTTGCCCAAAAAACAGTATTAGACACACTTCGCGAAGAAGGAGAA
GATGTAGAATTAGACGATGTTGTAAAACTCGGCTTCATGGGTTCACGTTGTGTAGAATCGGGTGGGCCTG
AACCCGGTGTTGGTTGTGCAGGCCGTGGTATCATTACTTCCATTAACCTGCTTGAACAATTAGGTGCTTA
TAACGAAGAACAACTTGAAACAGATTTTGTGTTTTATGATGTTTTGGGT
_known_contaminants_


## Translate all contaminants and make the DB. Strip orf bounds from
## the contaminant id's.
echo
echo "Translating contaminant sequences using FragGeneScan..."
run_FragGeneScan.pl -genome="$CONTAMFASTA" -out=contam -complete=0 -train=complete
cat contam.faa | sed 's/_[0-9]*_[0-9]*_[+]$//' > tmp.faa
mv tmp.faa contam.faa
echo
echo "Making the contaminants protein DB..."
makeblastdb -in contam.faa -dbtype prot -out nifHContamDB \
      -title "nifH amplicon contaminants known from literature and to user"
rm contam.{out,gff,ffn,faa} "$CONTAMFASTA"


## Translate ASVs
echo
echo "Translating ASVs..."
run_FragGeneScan.pl -genome="$asvFasta" -out=asvs -complete=0 -train=complete
cat asvs.faa | sed 's/_[0-9]*_[0-9]*_[+]$//' > tmp.faa
mv tmp.faa asvs.faa
rm asvs.{out,gff,ffn}


## BLASTp ASVS against the contaminants DB.
echo
echo "Searching NifH amino ASVs against the NifH contaminants database"
echo "using default parameters (expect=10) except will keep only the top"
echo "contaminant hit for each ASV..."
blastp -db nifHContamDB -query asvs.faa -outfmt 6 \
  -max_target_seqs 1 -max_hsps 1 -out blast.out
echo
echo "Now filtering for hits at >${aminoId}..."
awk -v pid="$aminoId" '$3 > pid' blast.out > "blast.${aminoId}id.out"
echo "Done!"
