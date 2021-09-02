#!/bin/bash

##
## For all the paired-end fastqs (full file paths) in the passed list, run
## cutadapt on them to trim nifH primers (and any adapters that are more 5' or
## 3' to the primers), saving the results in an output directory (OutDir, with
## same structure as in the listed fastqs).
##
## A few motes:
##  -- Careful, will clobber previous output in the OutDir.
##
##  -- Presumes fwd and rev nifH primers I got from Kendra.  These include some
##     wildcards. By default cutadapt does match IUPAC wildcards in "adapters"
##     (but you can tell it not to with --no-match-adapter-wildcards).
##
##  -- Note that cutadapt docs use "adapters" but we use it for primers, which
##     are between the Illumina adapters.  Because Illumina adapters are
##     "outside" of the primers, I give unanchored primer patterns to cutadapt.
##
##  -- Use parameters suggested in user documentation (search for http below)
##     for removing primers from paired-end reads.  Reads must have _both_
##     primers or they will be discarded. Because primers are interior to
##     Illumina adapters, I don't know or check for adapter sequences.  This
##     approach worked in my toy example with 5 data sets: > 97% of the reads
##     were retained for each set.
##
##  -- Also drop empty reads by using -m 1.  This is only to avoid an error in
##     plotting quality profiles from DADA2 -- chokes on empty reads.  Defer
##     length-based filtering until after DADA2 merges paired reads (e.g.Kendra
##     uses PEAR to restrict to reads 300-400nt).
##
##  -- Most parameters to cutadapt are left as default.  Including -e which
##     defaults to 0.1, so up to 10% of the "adapter" (the primer) can
##     mismatch. (Note that cutadapt doesn't distinguish between primers and
##     adapters.)  And see first note about wildcards.  Also, any default
##     quality filtering apparently throws out very few reads.  That's fine,
##     just leave it for the DADA2 filter and trim step.
##

## By default use Kendra's favorite primer. Note that they include IUPAC wildcards.
fwd="TGYGAYCCNAARGCNGA"  #  nifH_up_inner_F:   5’ – TGY GAY CCN AAR GCN GA – 3’
rev="ADNGCCATCATYTCNCC"  #  nifH_down_inner_R: 5’ – ADN GCC ATC ATY TCN CC – 3’

usage="
Use cutadapt to trim nifH paired-end Illumina reads their primers.
Usage:
\trunCutadapt.sh  FastqList  OutDir  [fwdPrimer] [revPrimer]

Paired files in FastqList of this form:
   path/to/prefix_R1_suffix.fastq.gz
   path/to/prefix_R2_suffix.fastq.gz
will be trimmed of their primers and stored in:
   OutDir/to/prefix_R1_suffix.trimmed.fastq.gz
   OutDir/to/prefix_R2_suffix.trimmed.fastq.gz
   OutDir/to/cutadapt.log
Note that the first part of the file path is replaced with OutDir.

By default these nifH primers are used:
\tforward:  5' - $fwd - 3'  for the R1 fastqs
\treverse:  5' - $rev - 3'  for the R2 fastqs
You may specify your own primers 5' to 3' with IUPAC codes (upper or lower case).

This script will discard read pairs that are missing one or both of the primers.
"

FastqList=$1
OutDir=$2
if [ ! -z "$3" ] ; then
    fwd=`echo "$3" | tr [:lower:] [:upper:]`
    rev=`echo "$4" | tr [:lower:] [:upper:]`
    ## Check both defined below.
fi

## Note on the cutadapt parameters used in the loop below:
## Pass params that will trim paired reads of their primers and discard pairs
## that do not contain primers.  Described at:
##   https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads
## The options below (in particular --discard-untrimmed) will cause a pair of
## reads to be discardedd if the R1 read lacks the forward primer or the R2 read
## lacks the reverse primer (i.e. both reads must have their respective 5'
## primers).  You can change this e.g. to require that both primers were missing
## by also passing --pair-filter=both (rather than 'any' which is the default).

if [ ! -f "$FastqList" ] || [ -z "$OutDir" ] ; then echo -e "$usage"; exit -1; fi
if [ -z `which cutadapt 2> /dev/null` ] ; then
    echo "Error: Cannot find cutadapt. Are you in a conda environment that includes it?"
    exit -1
fi
if [ -z "$rev" ] ; then echo "Error: You may not specify only a forward primer."; exit -1; fi

## Calculate reverse complements
legalIUPAC="ACGTURYSWKMBDHVN"  # From https://www.bioinformatics.org/sms/iupac.html
complIUPAC="TGCAAYRSWMKVHDBN"  # By hand and verified at www.reverse-complement.com
bad=`echo $fwd | tr -d "$legalIUPAC"`
if [ -n "$bad" ] ; then echo "Non IUPAC codes in your forward primer: $bad"; exit -1; fi
bad=`echo $rev | tr -d "$legalIUPAC"`
if [ -n "$bad" ] ; then echo "Non IUPAC codes in your reverse primer: $bad"; exit -1; fi

rcfwd=`echo $fwd | rev | tr "$legalIUPAC" "$complIUPAC"`
rcrev=`echo $rev | rev | tr "$legalIUPAC" "$complIUPAC"`

while read r1 ; do
    if [ ! -z `echo "$r1" | grep _R2_` ] ; then
        ## We will handle the R2 when we see the R1.
        continue
    fi
    ## Replace the first part of the path (which may or may not have
    ## a leading /) with OutDir.
    odir=`echo $(dirname "$r1") | sed -e 's:^\/::' -e 's:[^\/]*\/::'`
    odir="$OutDir/$odir"
    echo "Working on $r1 and its R2 file..."
    r2=`echo $r1 | sed 's/_R1_/_R2_/'`  # Not robust but ~okay.
    if [ ! -f "$r2" ] ; then
	echo "No R2 file for $r1"
	exit -1
    fi
    if [ ! -d "$odir" ] ; then mkdir -p "$odir"; fi
    out1="$odir/$(basename $r1 .fastq.gz).trimmed.fastq.gz"
    out2="$odir/$(basename $r2 .fastq.gz).trimmed.fastq.gz"
    if [ -f "$out1" ] ; then echo "$out1 already exists. Will clobber."; fi
    if [ -f "$out2" ] ; then echo "$out2 already exists. Will clobber."; fi
    ## As noted above, the primers are not anchored because there is Illumina
    ## adapter sequence up- and downstream of the primers.
    echo "### Working on $r1 and its R2 file." >> $odir/cutadapt.log
    echo "Will use these primers with cutadapt:" >> $odir/cutadapt.log
    echo -e "\tforward:  5'- ${fwd} -3'  revcomp 5'- ${rcfwd} -3'" >> $odir/cutadapt.log
    echo -e "\treverse:  5'- ${rev} -3'  revcomp 5'- ${rcrev} -3'" >> $odir/cutadapt.log
    echo >> $odir/cutadapt.log
    cutadapt -a "${fwd}...${rcrev}"  -A "${rev}...${rcfwd}"  -m 1  --discard-untrimmed \
             -o $out1  -p $out2  $r1  $r2 >> $odir/cutadapt.log
    echo -e "### Finished working on $r1 and its R2 file.\n\n" >> $odir/cutadapt.log
done < $FastqList
echo "Done!"
exit -1
