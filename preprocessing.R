setwd("D:/Projects/IBD Diagnosis Using Liquid Biopsies/Data Collection and Preprocessing/Dataset 1/")

### Converting Entrez Gene ID to Gene Symbol
library(org.Hs.eg.db)
library(annotate)

a <- read.table(file = "Entrez.txt", header = TRUE)
a <- as.vector(a$ENTREZ_GENE_ID)
c <- as.character(a)
b <- getSYMBOL(c(c), data = 'org.Hs.eg.db')
d <- as.data.frame(b)
colnames(d) <- "Gene Symbol"
d[is.na(d)]<- 0
d <- subset(d, d$`Gene Symbol` != 0)

write.table(x = d, file = "Entrez to Gene Symbol.txt", sep = "\t")

rm(d,a,b,c)
gc()

d <- read.delim(file = "Entrez to Gene Symbol.txt", header = T)
annot <- read.delim("GPL19109.txt", header = TRUE)
annot <- merge(d, annot, by = "Gene.Symbol")
write.table(x = annot, file = "GPL19109.txt", sep = "\t")

rm(d,annot)
gc()

### Creating Data File
data <- read.delim("GSE94648.txt", header = TRUE)

annot <- read.delim("GPL19109.txt", header = TRUE)

colnames(annot)[1] <- "ID_REF"

merged_data <- merge(annot, data, by = "ID_REF")

for(j in 2:ncol(merged_data)){
  if(typeof(merged_data[,j]) != "double"){print(paste("j=",j))}
}

merged_data[is.na(merged_data)]<- 0

data <- aggregate(by = list(unique.values = merged_data$Gene.Symbol)
                    , x = merged_data, FUN = mean)

data$ID_REF <- NULL
data$Gene.Symbol <- NULL
rownames(data) <- data$unique.values
data$unique.values <- NULL

if(max(data > 100)){
  data <- log2(data+1)
}

min(data)
max(data)

### Finishing Pre Processing
library(limma)

boxplot(data, las = 2, col = c(rep("blue", 11), rep("red", 24)),
       main = "Samples' Value Distributions")

pdf("Boxplot1.pdf")
boxplot(data, las = 2, col = "red",
        main = "Samples' Value Distributions")
dev.off()

dataset1 <- normalizeQuantiles(data)
min(dataset1)
max(dataset1)

write.table(x = data, file = "dataset1.txt", sep = "\t")
gr1 <- c(rep("UCa", 15), rep("UCi", 8), rep("Con", 22), 
         rep("CDa", 32), rep("CDi", 4))

### PCA
library(ggplot2)
dataset1 <- read.delim("dataset1.txt", header = TRUE)
gr <- c(rep("UCa", 15), rep("UCi", 8), rep("Con", 22), 
         rep("CDa", 32), rep("CDi", 4))

pc1 <- prcomp(dataset1)
pcr1 <- data.frame(pc1$r[,1:6], gr)

ggplot(pcr1,aes(x = PC2, y = PC3, colour=gr, shape = gr))+ geom_point(size=3)+
  theme_bw()+theme(legend.text=element_text(size=13,family = "serif"),
                   legend.title = element_text(size=15,family = "serif"),
                   axis.title =element_text(size=15,family = "serif"))

### Removing Outlier
data1outrem <- dataset1[,-c(3,5,8,14,24,28,33,34,44,45,65,66)]
gr <- c(rep("UCa", 11), rep("UCi", 8), rep("Con", 16), 
        rep("CDa", 30), rep("CDi", 4))

pc1 <- prcomp(data1outrem)
pcr1 <- data.frame(pc1$r[,1:6], gr)

library(ggplot2)
ggplot(pcr1,aes(x = PC1, y = PC2, colour=gr, shape = gr))+ geom_point(size=3)+
  theme_bw()+theme(legend.text=element_text(size=13,family = "serif"),
                   legend.title = element_text(size=15,family = "serif"),
                   axis.title =element_text(size=15,family = "serif"))

write.table(data1outrem, file = "data1outrem.txt", sep = "\t")


### Boxplot for Screening Intensity
gene_expression_data <- read.delim("dataset1.txt", header = TRUE)

mean_expression <- apply(gene_expression_data, 2, mean)

boxplot(mean_expression, main = "Boxplot of Mean Gene Expression", ylab = "Mean Expression")

stripchart(gene_expression_data, method = "jitter", add = TRUE, pch = 32, col = "blue")

for (i in seq_along(mean_expression)) {
  points(jitter(rep(i, nrow(gene_expression_data))), gene_expression_data[, i], col = "blue", pch = 16)
}

