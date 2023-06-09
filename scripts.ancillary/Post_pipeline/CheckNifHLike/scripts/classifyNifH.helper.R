#!/usr/bin/env Rscript

## Usage:
##   classifyNifH.R  positveExampleBlastx.tab  negativeExampleBlastx.tab
##

args <- commandArgs(trailingOnly=T)
blastxPosTab <- args[1]
blastxNegTab <- args[2]
cat("Loading positive blastx results",blastxPosTab,"and negative results", blastxNegTab,"\n")
cat("These are the tabular blastx results from searching a set of sequences against\n",
    "positive NifH examples and negative NifH examples.\n\n")
stopifnot(file.exists(blastxPosTab))
stopifnot(file.exists(blastxNegTab))

cnams <- c('qry','sbj','pid','alen','mm','gap','qs','qe','ss','se','eval','bits')
ptab <- read.table(blastxPosTab)
colnames(ptab) <- cnams
ntab <- read.table(blastxNegTab)
colnames(ntab) <- cnams

## Only use the best hits (by e-value).  match() finds the first matches for each qry,
## which is the best hit b/c that's how tabular blast output comes.
ptab <- ptab[match(unique(ptab$qry), ptab$qry),]
ntab <- ntab[match(unique(ntab$qry), ntab$qry),]
## Since each qry now appears just once, use as the row names.
rownames(ptab) <- ptab$qry
rownames(ntab) <- ntab$qry
origIds <- union(rownames(ptab),rownames(ntab))


onlyPosHits <- setdiff(ptab$qry,ntab$qry)
onlyNegHits <- setdiff(ntab$qry,ptab$qry)
bothHits    <- intersect(ntab$qry,ptab$qry)
cat("Num sequences with hits only to positive examples:    ",length(onlyPosHits),"\n")
cat("Num sequences with hits only to negative examples:    ",length(onlyNegHits),"\n")
cat("Num sequences with hits to both pos and neg examples: ",length(bothHits),"\n")
## Total universe of ASVs is unknown since I have only the blastx tables.
stopifnot(length(intersect(onlyPosHits,onlyNegHits))==0)


cat("\nClassifying sequences based on whether their best positive hits were 10x better\n",
    "(based on e-value) than their best negative hits, or vise versa.\n",
    "Hits with e-value >=1 will be ignored.\n",
    "Sequences without a clearly better positive or negative hit will be classified as 'unsure'.\n\n")

cat("Positive sequence best hits: Are e-values <1?:")
print(table(ptab$eval<1))
cat("\nNegative sequence best hits: Are e-values <1?:")
print(table(ntab$eval<1))


## Initialize the lists of accepted and rejected.  Heuristic: Sequence appears
## only among the positives(negatives) and has e-value < 1.
results <- list()
results$Positives <- rownames(subset(ptab, qry %in% onlyPosHits & eval < 1))
results$Negatives <- rownames(subset(ntab, qry %in% onlyNegHits & eval < 1))
###results$Unsure    <- setdiff(origIds, union(results$Positives, results$Negatives))
stopifnot(length(intersect(results$Positives,results$Negatives))==0)

## Add to the lists based on whether the hits in the positive are much better
## than in the negative, or vise versa.
## First restrict to sequences with hits to pos and neg (regardless of whether e-value was < 1).
ptab <- subset(ptab, qry %in% bothHits)
ntab <- subset(ntab, qry %in% bothHits)

## For each sequence the 'superiority' of the positive hit over the negative hit
## is defined as in the ARBitrator paper during the specificity phase (Heller et
## al. 2014): -log10(posEval) - -log10(negEval).  So superiority >1 is
## interpreted as the best positive hit has an e-value 10X better than the best
## negative hit.  And if superiority is < -1 then the negative hit is the much
## better hit.  (The blastx's were done assuming equal sizes of pos and neg DBs
## so e-values are comparable.)
stopifnot(rownames(ptab) == rownames(ntab))
superiorities <- -log10(ptab$eval) - -log10(ntab$eval)
ids <- rownames(subset(ptab[superiorities > +1,], eval < 1))  # 10X better pos hit has e-value < 1
results$Positives <- union(results$Positives, ids)

superiorities <- -log10(ntab$eval) - -log10(ptab$eval)
ids <- rownames(subset(ntab[superiorities > +1,], eval < 1))  # 10X better neg hit has e-value < 1
results$Negatives <- union(results$Negatives, ids)
stopifnot(length(intersect(results$Positives,results$Negatives))==0)

## What's left goes into unsure. This includes sequences for which abs(superiorities) <= 1.
results$Unsure    <- setdiff(origIds, union(results$Positives, results$Negatives))

cat("\nHere are the classification totals:\n")
print(sapply(results,length))

stopifnot(setequal(unlist(results), origIds))
writeLines(results$Positives, 'positives.ids')
writeLines(results$Negatives, 'negatives.ids')
writeLines(results$Unsure,    'unsure.ids')
cat("\nWrote out files {positive,negative,unsure}.ids.\n")
cat("Done.\n")
