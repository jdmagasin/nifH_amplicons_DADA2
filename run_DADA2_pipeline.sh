#!/bin/bash

OUTDIR="Data.nobackup"
PARAMS=$1
## If missing params file or if -h,--help, -help, etc., then print usage.
if [ -z "$PARAMS" ] || [ ! -z `echo $PARAMS | grep '\-h'` ] ; then
    cat <<-EOUSAGE
	 This script uses DADA2 to identify nifH amplicon sequence
	 variants (ASVs) in paired MiSeq data sets.

	 Usage:
            run_DADA2_pipeline.sh params.tsv > log.date.txt

	Will run the pipeline using the specified parameters, most of which specify
	how DADA2 should filter reads for quality. Saving output in the log file is
	recommended but not required.

	Before you run this script you must: (1) Install the required tools [see
	note at end]; (2) Use organizeFastqs.R to create a directory structure (in
	LinksToFastqs) that organizes the FASTQs into "processing groups":
	
	    Processing groups: The FASTQ files are paritioned into non-overlapping
	    sets that are run through DADA2 independently (during stage 4).  They 
	    also have their own error models (stage 3).  Although you can define
	    processing groups however you like, probably you want them to partition
	    the FASTQs by sequencing type (DNA vs. mRNA) and size fraction, e.g. so
	    that nifH genes vs. transcripts are processed separately for each size
	    fraction.  I prefer this approach because amplicon relative abundances
	    differ greatly with respect to sequence type and size fraction, and I
	    do not want that to influence DADA2's ASV inference. Second, we gain
	    confidence in the ASVs if we see them in different processing groups.
	    Third, it is easier to manage and write scripts for the outputs when
	    they are organized by directories that follow the processing groups,
	    e.g. path/to/output/DNA/Filter0.2_3.

	 The main stages of this script are:
	 1. Primer trimming using cutadapt.  Primers must be:
	      forward:  5’ – TGY GAY CCN AAR GCN GA – 3’
	      reverse:  5’ – ADN GCC ATC ATY TCN CC – 3’
	    which are hard coded in the runCutadapt.sh script that does this stage.
	    Trimmed FASTQs are placed in $OUTDIR/Data.trimmed

	 2. Identification of reads that likely encode NifH. These are used in the
	    next stage.  Output is stored in $OUTDIR/NifH_prefilter

	 3. Build error models for each processing group. Only use the NifH-like
	    reads found in stage 2 because I have found that excluding the non-NifH-
	    like reads (presumably off-target PCR products) results in better error
	    models.  The intuition is that non-nifH amplicons are highly variable
	    (compared to nifH) garbage that will make it harder for DADA2 to learn
	    the sequencing error models.  Error models are stored in 
	    $OUTDIR/ErrorModels

	 4. Run DADA2 on each processing group. The script that runs DADA2,
	    dada2_nifH_amplicons.R, has stages: Read quality filtering; ASV
	    inference ("denoising", but use the stage 3 error models); paired read
	    merging to create full length ASVs; removal of chimeric ASVs.  Output
	    is stored in $OUTDIR/Dada2PipeOutput/Out.<datestamp>

	 If you rerun this script, previous outputs will be reused. Usually this is
	 desired for stages 1-3 (e.g. there is no point in trimming primers again)
	 and allows you to experiment with DADA2 (stage 4) parameters. However, if
	 you want to rerun DADA2 but the datestamp already exists, just rename or
	 delete the previous Out.<datestamp> directory.

	 Required tools: This script depends on many external tools (R packages
	 cutadapt, HMMER3, ...) nearly all of which can be installed using
	 miniconda3.  The document INSTALL.txt describes how to get set up to run
	 the pipeline.

EOUSAGE
    exit -1
elif [ ! -f "$PARAMS" ] ; then
    echo "Must specify a parameters file. This can be the example params.csv file in"
    echo "which every parameter has been commented out (in which case defaults will be"
    echo "used)."
    exit -1
fi


## Absolute path of the scripts dir. 'scripts' must be next to run_DADA2_pipeline.sh.
SDIR="$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)/scripts"
CWD=`pwd`

echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
echo "Started DADA2 analysis scripts on" `date`
echo "Working directory is $CWD"
echo "Will use parameters file $PARAMS"
echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
echo

echo "All outputs will be in the directory $OUTDIR.  Previous intermediate"
echo "results will generally be reused, e.g. primer-trimmed reads and error models."
echo "The DADA2 analyses for each processing group will be redone and saved in"
echo "date-stamped directories."
mkdir -p "$OUTDIR"

## Use LinksToFastqs which is hard coded in organizeFastqs.R.  Could instead
## pass a name to this script but seems ~error prone.  One should not run the
## pipeline on different data paritions in the same working directory.
FQDIR=LinksToFastqs
if [ ! -d "$FQDIR" ] ; then
    echo "There should be a directory called $FQDIR that you created with the"
    echo "organizeFastqs.R script."
    exit -1
fi


echo
echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
echo "Building FASTQ lists that define which samples will be processed together."

## Note that this code is simplified compared to the Baltic data where there are
## separate lists for {DNA,RNA}x{sizeFracA,sizeFracB}.  See the Baltic script
## for future projects that require separate processing groups.
rflist="RawFastq.list"
if [ ! -f "$rflist" ] ; then
    find -L "$FQDIR" -name "*.fastq.gz" | sort > "$rflist"
else
    echo "$rflist already exists."
fi
echo "There are "`cat "$rflist" | wc -l`" FASTQ files in $rflist."

## Processing groups are defined as the distinct pathways up to but excluding
## the folder that holds the fastq's, and exclude the LinksToFastqs.  E.g.,
##   LinksToFastqs/foo/bar/s1/s1_R1.fastq.gz
##   LinksToFastqs/foo/bar/s1/s1_R2.fastq.gz
##   LinksToFastqs/foo/bif/s2/s2_R1.fastq.gz
##   LinksToFastqs/foo/bif/s2/s2_R2.fastq.gz
## has two processing groups, foo/bar and foo/bif.  For this script the key is
## that the processing groups encode partial paths so we can store and find
## outputs as we go.
while read line; do
   x="$(dirname $line)"
   echo "$(dirname $x)"
done < "$rflist" | sed 's/^'"$FQDIR"'\///' | sort | uniq > processingGroups.txt


echo
echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
echo "Running cutadapt. Presuming paired reads and the standard nifH primers."
echo "Results will be in Data.trimmed and each sample will have its own detailed"
echo "log from cutadapt."
if [ -d "$OUTDIR/Data.trimmed" ] ; then
    echo "Directory Data.trimmmed already exists so skipping cutadapt."
else
    ## Not necessary to use processing groups for cutadapt.
    ## Note that runCutadapt output fastq paths have the same directory
    ## structure as their inputs execpt the first folder is renamed as directed.
    $SDIR/runCutadapt.sh "$rflist" "$OUTDIR/Data.trimmed" 2>&1
    ## One summary across all processing groups.
    casum="$OUTDIR/Data.trimmed/summary.cutadapt.txt"
    echo "Summarizing all of the cutadapt logs (all processing groups) in $casum"
    echo "# Summary of reads that passed cutadapt filters." > "$casum"
    for x in `find "$OUTDIR/Data.trimmed" -name cutadapt.log`; do
        echo -n -e "$x:\t\t" >> "$casum"
        echo `grep '^Pairs written (passing filters):' $x | cut -d: -f2 | tr -d "[:space:]"` \
          >> "$casum"
    done
fi


echo
echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
echo "Identifying reads that look like NifH because those will lead to better DADA2"
echo "error models.  Will predict ORFs and then search them for Fer4_NifH domains"
echo "using HMMER and trusted cutoffs for the domain (defined in the hmm), and then"
echo "additionally filter for predicted domains >33 residues and bit score > 150."
echo "Will only use the R1 reads."
if [ -d "$OUTDIR/Data.NifH_prefilter" ] ; then
    echo "Looks like you already did this step. Skipping."
else
    cd "$OUTDIR"
    ## Not necessary to use processing groups.
    $SDIR/NifH_prefilter/run_NifH_prefilter.sh  Data.trimmed
    echo
    echo "Done. The NifH-like search results are in $OUTDIR/Data.NifH_prefilter."
    echo "Now summarizing how many R1 reads were processed, how many ORFs were predicted,"
    echo "and how many of the reads had NifH identified at cutoffs mentioned above."
    pfsum="Data.NifH_prefilter/summary.NifH_prefilter.txt"
    printf "%-30s%15s%15s%15s%15s\n" Sample R1.reads Orfs OrfsTC R1.reads.NifH > "$pfsum"
    for lf in `find Data.NifH_prefilter -name log.nifScan.txt`; do
        ## Print just the matching text (only on any lines matched).
        samp=`echo $(basename $(dirname $lf))`
        nseq=`cat $lf  | sed -n 's/^no\. of seqs: \([0-9][0-9]*\)/\1/p'`
        norf=`cat $lf  | sed -n 's/^Predicted \([0-9][0-9]*\) orfs\./\1/p'`
        ## Approximately how many ORFs beat the trusted cut offs for the NifH model.
        ## Slow so just count the domtab lines -- assume an ORF only gets 1 table entry.
        hmmtc=`gunzip -c "$(dirname $lf)/hmmsearch.Fer_NifH.domtab.gz" | grep -vc '^#'`
        nread=`cat $lf | sed -n 's/^Identified \([0-9][0-9]*\) reads with NifH domains\./\1/p'`
        printf "%-30s%15s%15s%15s%15s\n" $samp $nseq $norf $hmmtc $nread >> "$pfsum"
    done
    cd "$CWD"
fi


echo
echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
echo "Building error models from the forward reads that appeared to have NifH domains."
echo "Separate models will be built for each processing group."
while read pgrp; do
    pdesc=`echo "$pgrp" | tr '/' '.'`
    logfile="log.learnErrors.${pdesc}.txt"
    if [ -f "$OUTDIR/ErrorModels/$pgrp/${logfile}" ] ; then
        echo "Error models already exist for $pgrp ($logfile found). Skipping."
    else
        cd "$OUTDIR"
        echo " -- Working on error models for processing group $pgrp"
        ## pgrp encodes the path for the processing group.
        Rscript "${SDIR}/learnErrorsFromFastq.R"  "Data.NifH_prefilter/$pgrp" > "$logfile"
        mkdir -p "ErrorModels/$pgrp"
        mv errorModel.rds errorModel.pdf "$logfile" "ErrorModels/$pgrp"
        cd "$CWD"
    fi
done < processingGroups.txt


echo
echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
echo "The main event, DADA2!  Will run DADA2 separately for each processing group,"
echo "with independently inferred ASVs across DADA2 runs. A later stage of the"
echo "pipeline will relabel ASV id's to AUID's (ASV unique id's) for identifying"
echo "sequences found across DADA2 runs."
## One date stamp is used for all pgroups (even if DADA2 finishes the last group
## after midnight).  The stamp is so that we can have different DADA2 runs (with
## parameter tweaks) stored alongside each other.  (Previous stages of the
## pipeline change less frequently and their params should be finalized before
## playing with DADA2 params.)
STAMP=$(echo $(date) | awk -F' ' '{print $6 $2 $3}')
while read pgrp; do
    pdesc=`echo "$pgrp" | tr '/' '.'`
    logfile="log.dada2.${pdesc}.${STAMP}.txt"
    if [ -f "$OUTDIR/Dada2PipeOutput/$pgrp/${logfile}" ] ; then
        echo "DADA2 output already exists for $pgrp ($logfile found). Skipping."
    else
        cd "$OUTDIR"
        echo ""
        echo " #######################################################"
        echo " ## Running DADA2 on processing group $pgrp"
        echo " ##"
        echo " ##"
        ## Shared results directory for this processing group.  Add a time stamp for
        ## the DADA2 outdir so different runs (e.g. parameter tweaks) can be stored
        ## alongside each other.
        mkdir -p "Dada2PipeOutput/$pgrp"
        trimdir="Data.trimmed/$pgrp"
        dada2dir="Dada2PipeOutput/$pgrp/Out.${STAMP}"
        echo "For this DADA2 run the output will be in $dada2dir"
        Rscript "$SDIR/dada2_nifH_amplicons.R" "$trimdir" "$dada2dir"  \
                "ErrorModels/$pgrp/errorModel.rds" "../$PARAMS" \
          2> "$logfile"
        mv "$logfile" "${dada2dir}/"
        ## Tally reads retained at stage.  The first arg is the raw reads dir
        ## for this pgroup and takes awhile to load. Although the script caches
        ## the counts for subsequent runs, it uses a hard-coded cache name which
        ## would have collisions across pgroups. So delete the cache.
        $SDIR/countReadsAtStages.R "$CWD/$FQDIR/$pgrp" "$trimdir" "$dada2dir"
        rm -f readCountsCacheForScript_countReadsAtStages.RData
        rcdp="readCountsDuringProcessing.${pdesc}.${STAMP}.csv"
        rcdpPng="readCountsDuringProcessing.${pdesc}.${STAMP}.png"
        mv readCountsDuringProcessing.csv "$rcdp"
        mv readCountsDuringProcessing.png "$rcdpPng"
        mv "$rcdp" "$rcdpPng" "${dada2dir}/"
        echo " ##"
        echo " ## Finished DADA2 processing group $pgrp"
        echo " #######################################################"
        echo ""
        cd "$CWD"
    fi
done < processingGroups.txt


echo
echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
echo "Finished DADA2 analysis scripts on" `date`
echo "Working directory is $CWD"
echo "##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####"
exit 0
### END END END ####







