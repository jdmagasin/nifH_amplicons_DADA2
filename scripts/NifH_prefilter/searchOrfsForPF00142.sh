#!/bin/bash

## Hmmer3 search for NifH-like proteins in the amino acid fasta on stdin.
## Output only tabular results for NifH domains in hmmsearch.Fer_NifH.domtab.
## (Domain table because it has alignment bounds so I can get alignment
## length. The sequence tabular output doesn't have bounds.)
##
## Use --cut_tc ("use profile's TC trusted cutoffs to set all thresholding"
## [per-sequence and per-domain cut offs]) as I did in the pooled assembly BINF
## paper. From the HMMER3 userguide:
##   TC thresholds are generally considered to be the score of the
##   lowest-scoring known true positive that is above all known false positives.
##
## Don't output alignments (noali) and also throw out usual text output (>
## /dev/null) since the table is all we need.

## Expect the HMM to be in the same dir as this script.
NIFH_HMM="$(dirname $0)/Fer4_NifH.hmm"
if [ ! -f "$NIFH_HMM" ] ; then
    echo "Cannot find the NifH HMM. It should be alongside this script."
    exit -1
fi
hmmsearch  --cut_tc  --noali --cpu 6  \
	   --domtblout hmmsearch.Fer_NifH.domtab \
	   $NIFH_HMM  -  > /dev/null

