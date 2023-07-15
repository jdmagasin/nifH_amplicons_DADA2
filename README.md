# DADA2 pipeline for processing _nifH_ amplicon data sets

This repository contains our DADA2 pipeline for processing nitrogenase (_nifH_) amplicon data that was sequenced using paired-end MiSeq.

- **Installation:**  Documentation for installing (mini)conda and required conda and R packages.
- **Example:** A small example data set for testing your installation and learning how to create parameters files.
- **run_DADA2_pipeline.sh:**  Main script that runs the whole pipeline.
  
- **scripts:** Helper scripts used by run_DADA2_pipeline.sh
- **scripts.ancillary:**  Additional tools that are not part of the pipeline, organized by type:
   - ASVs_to_AUIDS:  For combining results for different data sets (run separately through the pipeline) into a since abundance table and FASTA with new sequence identifiers (AUIDs)
  - Annotation:  Several tools for annotating _nifH_ ASVs
  - Pre_pipeline:  Several tools for evaluating data sets before running them through the pipeline
  - Post_pipeline:  For identifying ASVs that are not likely _nifH_
  
- **bin:** Symbolic links to main scripts so they can be run easily from your unix-like shell
