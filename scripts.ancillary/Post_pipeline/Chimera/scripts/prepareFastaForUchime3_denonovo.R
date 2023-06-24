#!/bin/env Rscript

##
## Make a FASTA file with definition lines that are expected by
## vsearch's uchime3_denovo.  It is okay to pass a FASTA that
## contains more IDs than are in the abundance table. Only the
## sequences in the abundance table will be in the output FASTA.
##

args <- commandArgs(trailingOnly=T)
if (!file.exists(args[1:2]) || dir.exists(args[3])) {
    stop("Need an ASV abundance table, a FASTA file, and a name for the output ",
         "directory for creating temporary FASTA files for vsearch, one per sample. ",
         "All ASVs in the abundance table must be in the input FASTA.  Extra ASVs ",
         "may be in the input FASTA.")
}

cat("Reading abundances...")
## Rows are ASV ID's, columns are samples.
abundTab <- read.table(args[1], row.names=1, header=T, sep="\t", check.names=T)
cat("got", paste(dim(abundTab),c('ASVs','samples')),"\n")

## read.table's check.names=T (which is default) ensures no duplicated sample
## names (cols) and should also check for duplicates.  To be sure though...
stopifnot( length(unique(colnames(abundTab))) == ncol(abundTab) )

## Fetch ASV sequences.
cat("Loading libraries...\n")
suppressMessages(library(ShortRead))
cat("Loading ASV sequences...\n")
fasta <- readFasta(args[2])

## Ensure that every sequence named in the table is in the FASTA.
stopifnot(rownames(abundTab) %in% as.character(id(fasta)))
## Now restrict the FASTA to the ASVs in the table.
fasta <- fasta[ as.character(id(fasta)) %in% rownames(abundTab) ]

## Make a separate FASTA file for every sample so that uchime3 denovo can be run
## on each sample individually -- because chimerism should be tested separately
## for each PCR reaction (sample).
## - Each FASTA will have "size" information (ASV counts for the sample) at the
##   start of each definition line as required by uchime3_denovo.  Apparently
##   there cannot be a space after the ';' following the size (unlike for Bob
##   Edgar's uchime3_denovo).
##
outdir <- args[3]
if (!dir.create(outdir)) { stop("Failed to create ", outdir) }
for (samp in colnames(abundTab)) {
    asvTots <- abundTab[,samp]
    names(asvTots) <- rownames(abundTab)
    if (sum(asvTots) == 0) next;  # No reads in this sample
    asvTots <- asvTots[asvTots > 0]
    ## Fasta that has just the ASVs in this sample
    fas4samp <- fasta[ as.character(id(fasta)) %in% names(asvTots) ]
    idx <- match(as.character(id(fas4samp)), names(asvTots))
    deflines <- paste0('size=',asvTots[idx],';',names(asvTots)[idx])
    newfasta <- ShortRead(sread = sread(fas4samp), id = BStringSet(deflines))
    ## Check.
    stopifnot( sub('^size=[0-9]+;','',id(newfasta)) == id(fas4samp) )
    stopifnot( sread(newfasta) == sread(fas4samp) )
    fnam <- file.path(outdir,paste0(samp,".fasta"))
    writeFasta(newfasta, file=fnam)
    ##cat("Wrote out FASTA file",fnam,"\n")
}
