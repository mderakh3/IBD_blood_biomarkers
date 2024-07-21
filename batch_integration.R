setwd("D:/Projects/IBD Dignosis Project/Data Collection and Preprocessing/Batch Effect and Gene Filtering")

### Preparing IBD Metadata 
data1outrem <- read.delim(file = "data1outrem.txt", header = TRUE)
data4outrem <- read.delim(file = "data4outrem.txt", header = TRUE)

mdata <- merge(x = data1outrem, y = data4outrem, by = "row.names")
rownames(mdata) <- mdata$Row.names
mdata$Row.names <- NULL

min(mdata)
max(mdata)

gr1 <- c(rep("UCa", 11), rep("UCi", 8), rep("Con", 16), 
         rep("CDa", 30), rep("CDi", 4))
gr4 <- c(rep("Con", 11), rep("CDa", 15), rep("UCa", 15), 
         rep("Con", 10), rep("CDa", 16), rep("UCa", 13), rep("Con", 11),
         rep("CDa", 8), rep("UCa", 8), rep("Con", 12), rep("CDa", 8), 
         rep("UCa", 8))

gr <- c(gr1, gr4)
rm(gr1,gr2,gr3,gr4,gr5,gr6)
gc()

Datasets <- factor(c(rep(c("GSE94648","GSE119600"), 
                         times=c(ncol(data1outrem),ncol(data4outrem)))))

pc1 <- prcomp(mdata)
pcr1 <- data.frame(pc1$r[,1:6],Datasets, gr)

library(ggplot2)
ggplot(pcr1,aes(x = PC1, y = PC2, colour=Datasets, shape = gr))+ geom_point(size=3)+
  theme_bw()+theme(legend.text=element_text(size=13,family = "serif"),
                   legend.title = element_text(size=15,family = "serif"),
                   axis.title =element_text(size=15,family = "serif"))

write.table(mdata, "MetadataBeforeBatchoutrem.txt", sep = "\t")

### Batch Effect Removal
library(sva)

mdata <- read.delim("MetadataBeforeBatchoutrem.txt", header = T)
mdata.c <- as.data.frame(ComBat(as.matrix(mdata),Datasets))

min(mdata.c)
max(mdata.c)

pc1 <- prcomp(mdata.c)
pcr1 <- data.frame(pc1$r[,1:6],Datasets, gr)

ggplot(pcr1,aes(x = PC1, y = PC2, colour=Datasets, shape = gr))+ geom_point(size=3)+
  theme_bw()+theme(legend.text=element_text(size=13,family = "serif"),
                   legend.title = element_text(size=15,family = "serif"),
                   axis.title =element_text(size=15,family = "serif"))

write.table(mdata.c, "mdata.c.txt", sep = "\t")

### Gene Filtering
Q1 <- as.numeric(quantile(as.numeric(as.matrix(mdata.c)), probs = c(0.25)))
col_UCa <- which(gr == "UCa")
col_CDa <- which(gr == "CDa")
col_UCi <- which(gr == "UCi")
col_CDi <- which(gr == "CDi")
col_RA <- which(gr == "RA")
col_Con <- which(gr == "Con")

for (j in nrow(mdata.c):1) {
  if(sum(mdata.c[j,col_UCa] > Q1) < round(length(col_UCa)*0.7)){
    if(sum(mdata.c[j,col_CDa] > Q1) < round(length(col_CDa)*0.7 )){
      if(sum(mdata.c[j,col_RA] > Q1) < round(length(col_RA)*0.7)){
        if(sum(mdata.c[j,col_Con] > Q1) < round(length(col_Con)*0.7)){
          mdata.c <- mdata.c[-j,]
        }
      }
    }
  }
}
write.table(mdata.c, "Metadataoutrem.txt", sep = "\t")

### RA Data Gene Filtering
data6outrem <- read.delim(file = "data6outrem.txt", header = TRUE)

Q1 <- as.numeric(quantile(as.numeric(as.matrix(data6outrem)), probs = c(0.25)))
gr <- c(rep("Con", 30), rep("RA", 29), rep("Con", 13), 
        rep("RA", 89))

col_RA <- which(gr == "RA")
col_Con <- which(gr == "Con")

for (j in nrow(data6outrem):1) {
      if(sum(data6outrem[j,col_RA] > Q1) < round(length(col_RA)*0.7)){
        if(sum(data6outrem[j,col_Con] > Q1) < round(length(col_Con)*0.7)){
          data6outrem <- data6outrem[-j,]
    }
  }
}

write.table(data6outrem, "data6outremGenefil.txt", sep = "\t")
