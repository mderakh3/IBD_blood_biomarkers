### Functional Annotations

library(foreign)
library(igraph)
library(dplyr)

# Pearson Correlation
ex <- read.delim(file = "Metadataoutrem.txt", header = TRUE)
degs <- read.delim(file = "IBD_Con.txt", header = TRUE)
core_matrix <- ex[which(row.names(ex) %in% row.names(degs)),]
core_matrix <- t(core_matrix)
dim(core_matrix)
cor_mat <- cor(core_matrix, method = "pearson")
dim(cor_mat)
AA <- graph.adjacency(cor_mat, mode = "undirected", weighted = TRUE, diag = TRUE)
cor_edge_list <- get.data.frame(AA, "edges")
net <- cor_edge_list[cor_edge_list$weight != 1,]
for (i in nrow(net):1) {
  if (net[i,1] == net[i,2]) {net <- net[-i,]}
}

only_sig <- net[abs(net$weight) > 0.8,]
nodes <- unique(c(only_sig$from, only_sig$to))
write.table(only_sig, file = "CorrelationEdge.txt", sep = "\t", row.names = FALSE)

# Network Construction
core <- read.delim("CorrelationEdge.txt", header = T)
ppi <- read.delim("PPIEdges.txt", header = T)
netedges <- merge(ppi, core, by = "node1.to.node2")
netedges$node1.to.node2 <- NULL
write.table(netedges, "netedges.txt", sep = "\t")

# Pathview with Absolute Expressions
netgenes <- read.delim("Big Model default node.txt", header = T)
exdata <-  read.delim("Metadataoutrem.txt", header = T)
IBDColDataOutRem <- read.delim(file = "IBDColDataOutRem.txt", header = T)
exdata <- exdata[,IBDColDataOutRem$accision]
gr <- IBDColDataOutRem$Class
colnames(exdata) <- gr
pathviewInput <- exdata[netgenes$Genes,]
write.table(pathviewInput, "Gene Data.txt", sep = "\t")

# Pathview with Relative Expressions
netgenes <- read.delim("Big Model default node.txt", header = T)
exdata <-  read.delim("IBD_Con.txt", header = T)
pathviewInput <- exdata[netgenes$Genes,]
write.table(pathviewInput, "Gene Data (relative).txt", sep = "\t")

# Pathways bubble plot
data <- read.table("GO Network.txt", header = TRUE)
data <- data[,c(1,3,2)]
names(data)[2] <- paste("adj.P.value")

p = ggplot(data,aes(adj.P.value, Gene.Ontology))
p = p + geom_point()  
p = p + geom_point(aes(size = Gene.Ratio))
pbubble = p + geom_point(aes(size = Gene.Ratio, color = adj.P.value))
pr = pbubble + scale_color_gradient(low = "#EAA563", high = "#7F97B6")
pr = pr + labs(color = expression ("-Log10(adj.P.value)"),size="Gene.Ratio",x="-Log10(adj.P.value)",y="Gene Ontology",title="Gene Ontology Analysis of Network")+theme(plot.title=element_text(hjust=0.5))
pr = pr + theme_bw()
pr

