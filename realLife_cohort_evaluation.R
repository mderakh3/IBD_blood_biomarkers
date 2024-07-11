### Real Life Cohort Panel Evaluation

library("ggplot2")
library("plotly")
library("viridis")
library("hrbrthemes")
library("foreign")
library("pheatmap")
library("caret")
library("pROC")
library("dplyr")
library("doParallel")
library("e1071")
library("RANN")
library("class")
library("MASS")
library("sva")

###10-fold Cross Validation
data <- read.delim("qPCR 4 samples removed.txt", header = T)
data <- tibble::column_to_rownames(data,"Patient")

# PCA for Checking SVM
gr <- data$gr
data_pca <- data
data_pca$gr <- NULL
data_pca$fold <- NULL
data_pca <- t(data_pca)
pc1 <- prcomp(data_pca)
pcr1 <- data.frame(pc1$r[,1:2],gr, gr)
ggplot(pcr1,aes(x = PC1, y = PC2, colour=gr,shape = gr))+ geom_point(size=3)+
  theme_bw()+theme(legend.text=element_text(size=13,family = "serif"),
                   legend.title = element_text(size=15,family = "serif"),
                   axis.title =element_text(size=15,family = "serif"))

# SVM Cros fold Validation
CM_table <- data.frame("Sensitivity"=0,"Specificity"=0,"Accuracy"=0,
                       "fold"=0, "NO.features"=0, "method" = "a")

D1 <- data.frame("pred" = NA, "Ref" = NA)
for (i in 1:6) {
  set.seed(5783125)
  y_train = data[data$fold != i,]
  y_test = data[data$fold == i,]
  data_train <- data[row.names(y_train),]
  data_test <- data[row.names(y_test),]
  set.seed(100)
  
  new_train3 <- data_train[,c(2:4)]
  new_test3 <- data_test[,c(2:4)]

  svm.model2 <- svm(as.factor(y_train$gr) ~ ., data = new_train3)
  svm.pred3 <- predict(svm.model2, new_test3)
  D1 <- rbind(D1,data.frame("pred" = svm.pred3, "Ref" = as.factor(y_test$gr)))
  
  
}
D1 <- D1[-1,]
CM3 <- confusionMatrix(data =  factor(D1$pred),
                       reference = factor(D1$Ref), 
                       positive = "IBD")

CM_table <- c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"svm")
CM_table
write.table(CM_table,"10-fold CV.txt")

# ROC Curve 
ROC_table <- data.frame("response"= "a","predictor"= 0, 
                        "fold" = 0,"method" = "a")

for (i in 1:6) {
  set.seed(5783125)
  y_train = data[data$fold != i,]
  y_test = data[data$fold == i,]
  data_train <- data[row.names(y_train),]
  data_test <- data[row.names(y_test),]
  cl <- makeCluster(detectCores(), type ='PSOCK')
  registerDoParallel(cl)
  set.seed(100)
  
  #y_train$gr <- as.factor(y_train$gr)
  #y_test$gr <- as.factor(y_test$gr)
  
  train_gr <- as.factor(y_train$gr)
  test_gr <- as.factor(y_test$gr)
  
  new_train3 <- data_train[,c(1:3)]
  new_test3 <- data_test[,c(1:3)]
  
  cl <- makeCluster(detectCores(), type ='PSOCK')
  registerDoParallel(cl)
  
  #SVM
  
  svm.model <- svm(train_gr ~ ., data = new_train3, probability = TRUE)
  
  svm.pred <- predict(svm.model, new_test3, probability = TRUE)
  
  svm_pr  <- as.data.frame(attr(svm.pred ,"probabilities"))
  
  ROC <- data.frame("response"=as.character(test_gr),"predictor"= svm_pr$IBD,
                    "fold"= rep(i,length(svm.pred)),
                    "method" = rep("svm",length(svm.pred)))
  
  ROC_table <- rbind(ROC_table,ROC)

}

ROC_table <- ROC_table[-1,]

ME <- unique(ROC_table$method)
FO <- as.numeric(unique(ROC_table$fold))

ROC_svm <- ROC_table[ROC_table$method == "svm",]

library(pROC)
roc.svm <- roc(response =as.character(ROC_svm$response), 
               predictor = as.numeric(ROC_svm$predictor),
               smooth = FALSE,legacy.axes = TRUE, 
               levels = as.factor(c("Con","IBD")))

LAB <- c(paste("3-featured Panel: AUC =",round(roc.svm$auc,2)))

roc.list <- list(roc.svm)

g.list <- ggroc(roc.list, alpha = 0.5, size = 0.85, legacy.axes = TRUE, print.auc = TRUE)
g.list + scale_colour_manual(values = c("#7B94B4"),
                             labels= LAB)+theme_bw()+
  theme(legend.text=element_text(size=19,family = "sans"),
        legend.title = element_text(size=0,family = "sans"),
        axis.title =element_text(size=26,family = "sans"),
        axis.text=element_text(size=12))+
  theme(legend.position = c(0.7, 0.4),
        legend.direction = "vertical")
