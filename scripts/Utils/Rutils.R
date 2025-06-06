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
    msg <- paste0("\nThe pipeline tried to run dada2::filterAndTrim() with multithread = TRUE but an\n",
                  "error occurred that was not handled by filterAndTrim(). If you see a warning that\n",
                  "identifies cores that \"did not deliver results\", then you might have run out of\n",
                  "memory.  Check the \"filtered\" directory for the processing group with the crash:\n",
                  "  1. Does any sample have .trimmed.fastq.gz files that are too small (bytes or number\n",
                  "     of reads)?  This could indicate partial results.\n",
                  "  2. Does any sample have a modification time that is much later than all the other\n",
                  "     samples?  Possibly hours were spent waiting for the that failed on this sample.\n",
                  "Probably #1 or #2 will be true for one sample.  Splitting the processing group into\n",
                  "two ~equally sized groups should avoid the crash.  If not, then the crash probably has\n",
                  "a different cause. You should start by checking that the sample has valid FASTQs and\n",
                  "try running it through filterAndTrim() alone.\n\n")
    track <- tryCatch({ filterAndTrim(..., multithread = TRUE) },
                      error = function (e) { cat(msg); e})
    return(track)
}
