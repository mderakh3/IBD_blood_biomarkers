### Differential Expression Analysis

library(limma)
library(GEOquery)
library(reshape2)
library(plyr)
library(Biobase)
library(vsn)

ibd_ex <- read.delim(file = "Metadataoutrem.txt", header = TRUE) # IBD data
IBDColDataOutRem <- read.delim(file = "IBDColDataOutRem.txt", header = T)
ibd_ex <- ibd_ex[,IBDColDataOutRem$accision]
ibd_gr <- IBDColDataOutRem$Class

ra_ex <- read.delim(file = "data6outremGenefil.txt", header = TRUE) # RA data
RAColDataoutrem <- read.delim(file = "RAColDataoutrem.txt", header = T)
ra_ex <- ra_ex[,RAColDataoutrem$accision]
ra_gr <- RAColDataoutrem$Class

dif.ex <- function(ex, gr, con, ann=NULL, number=Inf) {
  if (! is.null(ann)) rownames(ex) <- 1:nrow(ex)
  gr <- factor(gr)
  mat  <- data.frame(G=gr)
  design <- model.matrix(~G+0,mat)
  colnames(design) <- levels(gr)
  fit <- lmFit(ex, design)
  cont.matrix <- makeContrasts(contrasts=con, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01,trend=TRUE,robust=TRUE)
  tT=topTable(fit2, adjust="fdr", number=number)
  if (! is.null(ann)) tT <- cbind(tT, ann[as.numeric(rownames(tT)),,drop=F])
  tT
}

IBD_Con <- dif.ex(ex = ibd_ex, gr = ibd_gr, con = "IBD-Con") # IBD vs Con
IBD_Con <- IBD_Con[IBD_Con$adj.P.Val < 0.05,]
write.table(IBD_Con, file = "IBD_Con.txt", sep = "\t")

RA_Con <- dif.ex(ex = ra_ex, gr = ra_gr, con = "RA-Con") # RA vs Con
RA_Con <- RA_Con[RA_Con$adj.P.Val < 0.05,]
write.table(RA_Con, file = "RA_Con.txt", sep = "\t")

shared_genes <- intersect(IBD_Con$ID, RA_Con$ID)
IBD_specific_genes <- IBD_Con[!IBD_Con$ID %in% shared_genes,]
write.table(IBD_specific_genes, file = "Biomarkers.txt", sep = "\t")






