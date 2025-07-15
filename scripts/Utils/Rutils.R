## Copyright (C) 2025 Jonathan D. Magasin

## Wrapper for dada2::filterAndTrim() to handle crashes that sometimes occur when running with
## multithread = TRUE. I suspect that the crash happens within one of mcmapply()'s parallelized
## calls to:
##     fastqPairedFilter() ---> FastqStreamer::yield()
## possibly because one of the cores fails to allocate memory.  This crash was observed when running
## on a 14-sample processing group.  For 13 samples the output .trimmed.fastq.gz all had rounghly
## the same last modification time stamp, but one sample (the largest) had a time stamp 7 hours
## later.  The dada2 log showed this warning:
##    In mclapply(seq_len(n), do_one, mc.preschedule = mc.preschedule,  :
##      scheduled cores 5, 11 did not deliver results, all values of the jobs will be affected
## and also a crash when mcmapply() tried to assign 14 names/paths to just 12 results received
## from the paralleized fastqPairedFilter() calls:
##    Error in names(answer) <- names1 : 
##      'names' attribute [14] must be the same length as the vector [12]
##    Calls: filterAndTrim -> mcmapply
## I suspect the failed allocation (or whatever caused cores 5 and 11 to not respond) prevented
## the last line of fastqPairedFilter():
##        return(invisible(c(reads.in = inseqs, reads.out = outseqs)))
## which caused the "did not deliver results" crash.
##
## This function wraps filterAndTrim() in a tryCatch() so that the user can be given some helpful
## context if a crash occurs.
##
## Relevant material:
##  * Closed dada2 issue:   https://github.com/benjjneb/dada2/issues/1915
##  * dada2::filterAndTrim() docs:  See Details and description on 'n' parameter
##  * ShortRead::yield() docs:  See 'n' and 'readerBlockSize' parameters.
##
Graceful_filterAndTrim <- function(...)
{
    ## These are errors in the source of filterAndTrim() in dada2 1.30.0.  Most of them are for
    ## problems with parameters or similarly simple problems.  If any of these appear, then do not
    ## report 'msg' below which would add confusion.
    IsKnownError <- function(e) {
        knownErrors <- c("File paths must be provided as character vectors.",
                         "Some input files do not exist.",
                         "Every input file must have a corresponding output file.",
                         "All output files must be distinct.",
                         "Output files must be distinct from the input files.",
                         "Output files for the reverse reads are required.",
                         "File paths (rev/filt.rev) must be provided as character vectors.",
                         "Some input files (rev) do not exist.",
                         "Paired forward and reverse input files must correspond.",
                         "Every input file (rev) must have a corresponding output file (filt.rev).",
                         "All output files must be distinct.",
                         "Output files must be distinct from the input files.",
                         ## This error is not within filterAndTrim() but can occur.  The solution is to
                         ## set id.field in the parameters file.
                         "Couldn't automatically detect the sequence identifier field in the fastq id string.")
                    # These errors deal with multithreaded processing.  DO report 'msg' if they occur.
                    #    "These are the errors (up to 5) encountered in individual cores..."
                    #    "Some input files were not processed, perhaps due to memory issues. Consider lowering ncores."
        return( conditionMessage(e) %in% knownErrors )
    }

    msg <- paste0("\nThe pipeline tried to run dada2::filterAndTrim() with multithread = TRUE but an\n",
                  "error occurred that was not handled by filterAndTrim(). There should be an error message\n",
                  "just above that explains the problem.  Additionally you can check the log for a warning\n",
                  "that identifies cores that \"did not deliver results\".  This might indicate that you\n",
                  "ran out of memory.  Also check the \"filtered\" directory for the processing group with\n",
                  "the crash:\n",
                  "  1. Does any sample have .trimmed.fastq.gz files that are too small (bytes or number\n",
                  "     of reads)?  This could indicate partial results.\n",
                  "  2. Does any sample have a modification time that is much later than all the other\n",
                  "     samples?  Possibly hours were spent waiting for the core that failed on this sample.\n",
                  "Probably #1 or #2 will be true for one sample.  Splitting the processing group into\n",
                  "two ~equally sized groups should avoid the crash.  If not, then the crash probably has\n",
                  "a different cause. You should start by checking that the sample has valid FASTQs and\n",
                  "try running them through filterAndTrim() alone.\n\n")
    track <- tryCatch({ filterAndTrim(..., multithread = TRUE) },
                      error = function (e) {
                          cat("\nThe following error occurred during filterAndTrim():\n",
                              conditionMessage(e), "\n")
                          if (!IsKnownError(e)) { cat(msg) }
                          NULL
                      })
    return(track)
}
