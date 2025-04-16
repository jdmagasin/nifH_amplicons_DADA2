# DADA2 pipeline for processing _nifH_ amplicon sequencing data

This repository contains our DADA2 pipeline for processing nitrogenase (_nifH_) amplicon data that
were sequenced using paired-end MiSeq.  Features of our pipeline include:

  - A text file provides key DADA2 [parameters](Example/params.example.csv) (e.g., for filterAndTrim()).  No coding is required to
    run the pipeline.

  - Samples are input using a [FASTQ map](Example/fastqMap.tsv) (text table) that includes structured descriptions,
    e.g. whether the sample was DNA or RNA, the collection site, and the size fraction.  The
    pipeline uses the descriptions to partition the samples into "processing groups" which are run
    separately through DADA2.  Advantages of this approach include:
    
       1. Results organized hierarchically e.g., DNA/Station41/SizeFrac_0.2.
       2. Faster processing by DADA2.
       3. Separate error models for each processing group.

  - Error models can be precalculated using only the reads that appear to be _nifH_, not PCR
    artifacts.  On average, this results in up to a few thousand more reads retained in each sample
    and fewer ASVs (by as much as a few K).  [These plots](Readmore/nifH_error_models.md) compare
    results with vs. without using _nifH_ error models.

  - Automatically runs cutadapt.

  - Simple set up by creating a [conda environment](https://docs.conda.io/en/latest/index.html) called DADA2_nifH.

Once the [parameters file](Example/params.example.csv) and [FASTQ map](Example/fastqMap.tsv) are created, two commands launch the pipeline:

```bash session
(DADA2_nifH) [jmagasin@thalassa]$ organizeFastqs.R fastqMap.tsv
(DADA2_nifH) [jmagasin@thalassa]$ run_DADA2_pipeline.sh params.csv &> log.pipeline.txt &
```

All scripts run from the command-line in a Unix/Linux shell (BASH recommended) and provide
documentation when run with no parameters.  For example, the documentation from the main pipeline
script, run_DADA2_pipeline.sh, is included below.


## Can I use the pipeline with other types of amplicon data?

Yes.

1. Set _forward_ and _reverse_ in your parameters file to your PCR primers. Otherwise
   the pipeline will default to primers _nifH1_ and _nifH4_.

2. Disable precalculated error models by setting _skipNifHErrorModels_ to _true_.

We have used our pipeline with 16S rRNA MiSeq data using these steps.  No coding
changes were required.


## Repository contents

- **[Installation](Installation/INSTALL.txt):**  Documentation for installing (mini)conda and required conda and R packages.
- **[Example](Example/EXAMPLE.txt):** A small example data set for testing your installation and learning how to create the parameter file and the table of input FASTQ files.  A second example shows how to compare pipeline outputs, e.g., ASVs and total abundances obtained when using different pipeline parameters. This is an important step before committing to the ASVs that you will use in your analysis.
- **run_DADA2_pipeline.sh:**  Main script that runs the whole pipeline.
  
- **scripts:** Helper scripts used by run_DADA2_pipeline.sh.
- **scripts.ancillary:**  Additional tools you might find useful, mainly for quality-filtering and annotating ASVs after running the pipeline.  Most tools include an Example subdirectory.
   - ASVs_to_AUIDS:  For combining results across runs of the pipeline: Merge abundance tables and assign new sequence identifiers (AUIDs).
  - Annotation:  Several tools for annotating _nifH_ ASVs.
  - Pre_pipeline:  Several tools for evaluating data sets before running them through the pipeline.
  - Post_pipeline:  Quality filters for identifying ASVs that are probably not _nifH_.
  
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

 If you rerun this script, previous outputs will be reused. This behavior
 is usually desired for stages 1-3 (e.g., there is no reason to trim
 primers again) and allows you to experiment with DADA2 parameters (stages
 4-9).  However, if you want to rerun DADA2 but the datestamp already
 exists, you can (1) rename previous Out.<datestamp> directory, or (2) add
 to your parameters file a Dada2OutdirTag.  For example, if you include
 "Dada2OutdirTag, truncLen187" then output will go in directory
 Out.truncLen187.<datestamp>.

 Required tools: This script depends on many external tools (R packages
 cutadapt, HMMER3, ...) nearly all of which can be installed using
 miniconda3.  The document INSTALL.txt describes how to get set up to run
 the pipeline.
```

***

_Copyright (C) 2023 Jonathan D. Magasin_
