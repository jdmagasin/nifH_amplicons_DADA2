#!/bin/bash

## Stream me a fastq.  I will use FragGeneScan (old version) to predict ORFs.
## (No quality info in the FGS output of course.)

## Unfortunately FGS (my old copy at least) cannot use stdin as 'genome'
## so I have to save the fasta to a temp file.
TMP=$(mktemp --tmpdir Dada2pipe.fastq2orf.fasta.XXXXXXXXXXXX)

## Key off of the '+' line to grab the sequence ID and the sequence.
cat - | grep -B 2 '^+$' | egrep -v '^(--|\+)' | sed 's/^@/>/' > $TMP

## Error rate of 1% selected because the Guam reads often have phred 38.
#FragGeneScan -genome=$TMP -out=orfs -complete=0 -train=illumina_10
run_FragGeneScan.pl -genome=$TMP -out=orfs -complete=0 -train=illumina_10 -thread=6
rm $TMP
