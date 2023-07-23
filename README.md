# DADA2 pipeline for processing _nifH_ amplicon sequencing data

This repository contains our DADA2 pipeline for processing nitrogenase (_nifH_) amplicon data that were sequenced using paired-end MiSeq.

All scripts run from the command-line in a Unix/Linux shell (BASH recommended) and provide documentation when run with no parameters.  For example, the documentation from the main pipeline script, run_DADA2_pipeline.sh, is included below.


## Repository contents

- **[Installation](Installation/INSTALL.txt):**  Documentation for installing (mini)conda and required conda and R packages.
- **[Example](Example/EXAMPLE.txt):** A small example data set for testing your installation and learning how to create the parameter file and table of input FASTQ files.
- **run_DADA2_pipeline.sh:**  Main script that runs the whole pipeline.
  
- **scripts:** Helper scripts used by run_DADA2_pipeline.sh.
- **scripts.ancillary:**  Additional tools that are not part of the pipeline, organized within four subdirectories.  Most tools includes an Example (subdirectory).
   - ASVs_to_AUIDS:  For combining results for different data sets (run separately through the pipeline) into one abundance table and FASTA with new sequence identifiers (AUIDs).
  - Annotation:  Several tools for annotating _nifH_ ASVs.
  - Pre_pipeline:  Several tools for evaluating data sets before running them through the pipeline.
  - Post_pipeline:  For identifying ASVs that are not likely _nifH_.
  
- **bin:** Symbolic links to main scripts so they can be run easily from your unix-like shell.


## Documentation for run_DADA2_pipeline.sh

```
This script uses DADA2 to identify nifH amplicon sequence
 variants (ASVs) in paired MiSeq data sets.

 Usage:
            run_DADA2_pipeline.sh params.csv > log.date.txt

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
 1. Primer trimming using cutadapt.  The default primers are:
      forward:  5’ – TGY GAY CCN AAR GCN GA – 3’
      reverse:  5’ – ADN GCC ATC ATY TCN CC – 3’
    but you can specify other primers by setting 'forward' and 'reverse'
    in the parameters file.  The script runCutadapt.sh uses the specified
    primers and puts the trimmed FASTQs in Data.nobackup/Data.trimmed

 2. Identification of reads that likely encode NifH. These are used in the
    next stage.  Output is stored in Data.nobackup/NifH_prefilter

 3. Build error models for each processing group. Only use the NifH-like
    reads found in stage 2 because I have found that excluding the non-NifH-
    like reads (presumably off-target PCR products) results in better error
    models.  The intuition is that non-nifH amplicons are highly variable
    (compared to nifH) garbage that will make it harder for DADA2 to learn
    the sequencing error models.  Error models are stored in 
    Data.nobackup/ErrorModels

 4. Run DADA2 on each processing group. The script that runs DADA2,
    dada2_nifH_amplicons.R, has stages: Read quality filtering; ASV
    inference ("denoising", but use the stage 3 error models); paired read
    merging to create full length ASVs; removal of chimeric ASVs.  Output
    is stored in Data.nobackup/Dada2PipeOutput/Out.<datestamp>

 If you rerun this script, previous outputs will be reused. Usually this is
 desired for stages 1-3 (e.g. there is no point in trimming primers again)
 and allows you to experiment with DADA2 (stage 4) parameters. However, if
 you want to rerun DADA2 but the datestamp already exists, just rename or
 delete the previous Out.<datestamp> directory.

 Required tools: This script depends on many external tools (R packages
 cutadapt, HMMER3, ...) nearly all of which can be installed using
 miniconda3.  The document INSTALL.txt describes how to get set up to run
 the pipeline.
```

***

_Copyright (C) 2023 Jonathan D. Magasin_
