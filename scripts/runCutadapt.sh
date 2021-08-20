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
##  -- Hack! For some reason I run into failed checksums (gzip.py) if I run this
##     script on Kendra's fastq.gz's, but not on my own.  Perhaps we use
##     different compression libraries. To work around, uncompress her files.
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
##     plotting quality profiles from DADA2 -- chokes on empy reads.  Defer
##     length-based filtering until after DADA2 merges paired reads (e.g.Kendra
##     uses PEAR to restrict to reads 300-400nt).
##
##  -- Most parameters to cutadapt are left as default.  Including -e which
##     defaults to 0.1, so up to 10% of the "adapter" (the primer) can
##     mismatch. (Note that cutadapt doesn't distinguish between primers and
##     adapters.)  And see first note about wildcards.  Also, any default
##     quality filtering apparently throws out very few reads.  That's fine,
##     just leave it for DADA2 filter and trim step.
##

FastqList=$1
OutDir=$2

## Primers Kendra provided are as follows, and note that they include IUPAC wildcards.
fwd="TGYGAYCCNAARGCNGA"  #  nifH_up_inner_F:   5’ – TGY GAY CCN AAR GCN GA – 3’
rev="ADNGCCATCATYTCNCC"  #  nifH_down_inner_R: 5’ – ADN GCC ATC ATY TCN CC – 3’

## Their reverse complements (from dada2::rc() and checkd with www.reverse-complement.com).
rcfwd="TCNGCYTTNGGRTCRCA"
rcrev="GGNGARATGATGGCNHT"

usage="
Use cutadapt to trim nifH paired-end Illumina reads their primers.
Usage:
\trunCutadapt.sh  FastqList  OutDir

Paired files in FastqList of this form:
   path/to/prefix_R1_suffix.fastq.gz
   path/to/prefix_R2_suffix.fastq.gz
will be trimmed of their primers and stored in:
   OutDir/to/prefix_R1_suffix.trimmed.fastq.gz
   OutDir/to/prefix_R2_suffix.trimmed.fastq.gz
   OutDir/to/cutadapt.log
Note that the first part of the file path is replaced with OutDir.

The R1 and R2 fastqs should be for nifH sequences that have these primers:
\tforward:  5' - $fwd - 3'
\treverse:  5' - $rev - 3'
This script will discard read pairs that are missing one or both of the primers.
You should be in a conda environment that includes cutadapt.
"

## Note on the cutadapt parameters used in the loop below:
## Pass params that will trim paired reads of their primers and discard pairs
## that do not contain primers.  Described at:
##   https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads
## The options below (in particular --discard-untrimmed) will cause a pair of
## reads to be discardedd if the R1 read lacks the forward primer or the R2 read
## lacks the reverse primer (i.e. both reads must have their respective 5'
## primers).  You can change this e.g. to require that both primers were missing
## by also passing --pair-filter=both (rather than 'any' which is the default).

if [ ! -x `which cutadapt` ] ; then
    echo "Error: Cannot find cutadapt. Are you in a conda environment that includes it?"
    echo -e "$usage"
    exit -1
fi
if [ ! -f "$FastqList" ] || [ -z "$OutDir" ] ; then
    echo -e "$usage"
    exit -1
fi

chackdir=`mktemp -d --tmpdir "${USER}_cutadapthack_XXXXX"`
if [ "$?" -ne 0 ] ; then
    echo "Error:  Failed to create temporary directory for runCutadapt.sh."
    exit -1
fi

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
    ## Hack! Uncompress Kendra's fastq.gz to work around failed checksum
    ## (gzip.py) probably due to she and I using different zlibs.
    tr1="$chackdir/$(basename $r1 .gz)"
    tr2="$chackdir/$(basename $r2 .gz)"
    zcat $r1 > $tr1
    zcat $r2 > $tr2
    ## As noted above, the primers are not anchored because there is Illumina
    ## adapter sequence up- and downstream of the primers.
    echo "### Working on $r1 and its R2 file." >> $odir/cutadapt.log
    cutadapt -a "${fwd}...${rcrev}"  -A "${rev}...${rcfwd}"  -m 1  --discard-untrimmed \
             -o $out1  -p $out2  $tr1  $tr2 >> $odir/cutadapt.log
    echo -e "### Finished working on $r1 and its R2 file.\n\n" >> $odir/cutadapt.log
    rm $tr1 $tr2
done < $FastqList
echo "Done!"
exit -1
