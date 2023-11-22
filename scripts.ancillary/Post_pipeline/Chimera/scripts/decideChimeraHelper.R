#!/bin/env Rscript

## Copyright (C) 2023 Jonathan D. Magasin

##
## Called by check_chimera_denovo.sh after it runs uchime3_denovo on each of the
## samples (columns in the abundance table).
##
## An ASV is defined as a chimera if uchime3_denovo called in a chimera in a
## majority of samples.  Ties --> not a chimera.
##

## Returns a character vector (possibly empty) of ASV IDs in FASTA file 'fas'.
ReadAsvIds <- function(fas) {
    lns <- readLines(fas)
    sub('^>([^\\s]+)', '\\1', lns[grep("^>",lns)])
}

TallyAsvs <- function(pat,desc) {
    flist <- list.files(pattern=pat)
    df <- data.frame(V1='dummyASV',V2=0)
    if (length(flist) > 0) {
        df2 <- data.frame(table(unlist(lapply(flist, ReadAsvIds))))
        if (nrow(df2) > 0) df <- df2
    }
    colnames(df) <- c('id',desc)
    df
}

df.chim    <- TallyAsvs('\\.chimera\\.fasta',    'samps.chimera')
df.nonchim <- TallyAsvs('\\.nonchimera\\.fasta', 'samps.nonchimera')

voteTab <- merge(df.chim, df.nonchim, all.x=T, all.y=F)  # need ASVs with chimera calls
stopifnot(nrow(voteTab) > 0)
stopifnot(!is.na(voteTab$samps.chimera))
voteTab[is.na(voteTab)] <- 0
voteTab$VotedChimera <- (voteTab$samps.nonchimera < voteTab$samps.chimera)
write.csv(voteTab, 'chimeraVotes.csv', quote=F, row.names=F)

## It is helpful to have ID files of the chimera and non-chimera that
## check_chimera_denovo.sh can use to extractFasta.pl.
chimIds <- voteTab$id[voteTab$VotedChimera]
writeLines(paste(chimIds), 'chimera.ids')
writeLines(paste(setdiff(union(df.chim$id,df.nonchim$id), chimIds)), 'nonchimera.ids')
