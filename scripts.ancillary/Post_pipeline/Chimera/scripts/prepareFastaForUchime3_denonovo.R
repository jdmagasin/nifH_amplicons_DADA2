#!/bin/env Rscript

##
## Make a FASTA file with definition lines that are expected by
## vsearch's uchime3_denovo.  It is okay to pass a FASTA that
## contains more IDs than are in the abundance table. Only the
## sequences in the abundance table will be in the output FASTA.
##

args <- commandArgs(trailingOnly=T)
if (!file.exists(args[1:2]) || file.exists(args[3])) {
    stop("Need an ASV abundance table, a FASTA file, and a name ",
         "for the output FASTA file which does not yet exist.",
         "The FASTA can have more sequences than are in the table",
         "but IDs in the table will be used for the output FASTA.")
}

cat("Reading abundances...\n")
## Rows are ASV ID's so next line sums the counts for each ASV.
asvTots <- rowSums(read.table(args[1], row.names=1, header=T, sep="\t"))
stopifnot(table(names(asvTots))==1)  # no duplicate rows/ASVs

## Fetch ASV sequences.
cat("Loading libraries...\n")
suppressMessages(library(ShortRead))
cat("Loading ASV sequences...\n")
fasta <- readFasta(args[2])

## Ensure that every sequence named in the table is in the FASTA.
stopifnot(names(asvTots) %in% as.character(id(fasta)))
## Now restrict the FASTA to what is needed.
fasta <- fasta[ as.character(id(fasta)) %in% names(asvTots) ]

## Make new FASTA in which "size" info (ASV counts) is at the start of each
## definition line as required by uchime3_denovo.  Apparently there cannot be a
## space after the ';' following the size (unlike for Bob Edgar's
## uchime3_denovo).
idx <- match(as.character(id(fasta)), names(asvTots))
deflines <- paste0('size=',asvTots[idx],';',names(asvTots)[idx])
newfasta <- ShortRead(sread = sread(fasta), id = BStringSet(deflines))
## Check.
stopifnot( sub('^size=[0-9]+;','',id(newfasta)) == id(fasta) )
stopifnot( sread(newfasta) == sread(fasta) )

writeFasta(newfasta, file=args[3])
cat("Wrote out FASTA file",args[3],"\n")
