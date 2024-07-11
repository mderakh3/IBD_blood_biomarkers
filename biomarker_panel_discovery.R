### In Silico Panel Discovery

library(reshape2)
library(plyr)
library(Biobase)
library(randomForest)
library(caret)
library(data.table)
library(dplyr)
library(doParallel)
library(rpart)
library(e1071)
library(RANN)
library(class)
library(MASS)
library(foreign)
library(glmnet)

ibd_ex <- read.delim(file = "Metadataoutrem.txt", header = TRUE)
IBDColDataOutRem <- read.delim(file = "IBDColDataOutRem.txt", header = T)
ex <- ex[,IBDColDataOutRem$accision]
gr <- IBDColDataOutRem$Class

Biomarkers <- read.delim(file = "Biomarkers.txt", header = TRUE)
IBDBioData <- ex[Biomarkers$Genes,]
FSI <- as.data.frame(t(IBDBioData)) #FSI: Feature Selection Input
FSI <- as.matrix(FSI[IBDColDataOutRem$accision,])
gr <- as.factor(IBDColDataOutRem$Class)

cl <- makeCluster(detectCores(), type="PSOCK")
registerDoParallel(cl)

set.seed(100478)
cv.lasso_cont <- cv.glmnet(FSI, gr, family = "binomial", nfolds = 10 ,
                           alpha = 1,  standardize = TRUE, type.measure = "mse")
Coef_cont <- coef(cv.lasso_cont, s = cv.lasso_cont$lambda.min)
registerDoSEQ()

Coef_cont <- as.matrix(Coef_cont)
Coef1_cont <- data.frame("Gene" = row.names(Coef_cont),"coef" = Coef_cont[,1],
                         "abs" = abs(Coef_cont[,1]))
Coef1_cont <- Coef1_cont[Coef1_cont$coef !=0,]
Coef1_cont <- Coef1_cont[-1,]
write.table(Coef1_cont, "IBDFeatures.txt", sep = "\t")

data <- FSI[,Coef1_cont$Gene] 
data$fold <- c(rep(1:5, 38), rep(1:2, 1))
mark <- colnames(data)

library(gtools)
COM3 <- combinations(n = 15,r = 3, v = mark)
COM4 <- combinations(n = 15,r = 4, v = mark)
COM5 <- combinations(n = 15,r = 5, v = mark)
COM6 <- combinations(n = 15,r = 6, v = mark)
COM7 <- combinations(n = 15,r = 7, v = mark)
COM8 <- combinations(n = 15,r = 8, v = mark)
COM9 <- combinations(n = 15,r = 9, v = mark)
COM10 <- combinations(n = 15,r = 10, v = mark)
COM11 <- combinations(n = 15,r = 11, v = mark)
COM12 <- combinations(n = 15,r = 12, v = mark)
COM13 <- combinations(n = 15,r = 13, v = mark)
COM14 <- combinations(n = 15,r = 14, v = mark)

set.seed(83345786)
test_fold <- sample(1:5,1)
ColData_train <- ColData[ColData$fold1 != test_fold,]
ColData_test <- ColData[ColData$fold1 == test_fold,]
test_data <- data[ColData_test$accision,]
train_data <- data[ColData_train$accision,]

row.names(train_data)[27] == ColData_train$accision[27]
row.names(train_data)[nrow(ColData_train)] == ColData_train$accision[nrow(ColData_train)]

row.names(test_data)[27] == ColData_test$accision[27]
row.names(test_data)[nrow(ColData_test)] == ColData_test$accision[nrow(ColData_test)]

ColData_train$gr <- as.factor(ColData_train$gr)
ColData_test$gr <- as.factor(ColData_test$gr)

list_mark <- list(COM3, COM4, COM5, COM6, COM7, COM8, COM9, COM10, COM11, COM12,
                  COM13, COM14)

cl <- detectCores()
registerDoParallel(25)  

CM_table1 <- data.frame("Sensitivity"= 0, # 80-20 Combination Machine Learning
                        "Specificity"= 0,
                        "Accuracy"= 0,
                        "NO.features"= 0, "method" = "a",
                        C1 = 0, C2 = 0)

test_data <- as.matrix(test_data)
train_data <- as.matrix(train_data)

for (c1 in seq(length(list_mark))) {
  M1 <- list()
  for(t in seq(nrow(list_mark[[c1]]))){
    M1[[t]] <- list_mark[[c1]][t,]
  }
  length(M1) == nrow(list_mark[[c1]])
  
  foreach (M = M1) %dopar% {
    library(randomForest)
    library(caret)
    library(rpart)
    library(MASS)
    library(e1071)
    train_data3 <- train_data[,M]
    test_data3 <- test_data[,M]
    set.seed(358632)
    RF3 <- randomForest(x = train_data3, 
                        y = factor(ColData_train$gr),
                        ntree=1000, mtry = length(M))
    pred_test3 <- predict(RF3 ,newdata = test_data3 ,type="class")
    CM3 <- confusionMatrix(data = factor(pred_test3),
                           reference = factor(ColData_test$gr),
                           positive = "IBD")
    CM_table <- data.frame("Sensitivity"= CM3$byClass[1],
                           "Specificity"= CM3$byClass[2],
                           "Accuracy"= CM3$overall[1],
                           "NO.features"= length(M),
                           "method" = "RF",C1 = c1)
    svm.model3 <- svm(ColData_train$gr ~ ., data = train_data3)
    svm.pred3 <- predict(svm.model3, test_data3)
    CM3 <- confusionMatrix(data =  factor(svm.pred3),
                           reference = as.factor(ColData_test$gr), 
                           positive = "IBD")
    CM_table <- rbind(CM_table, data.frame("Sensitivity"= CM3$byClass[1],
                                           "Specificity"= CM3$byClass[2],
                                           "Accuracy"= CM3$overall[1],
                                           "NO.features"= length(M),
                                           "method" = "SVM",C1 = c1))
    
  } -> CM_T
  
  results <- CM_T
  
  for (z in seq(length(results))) {
    results[[z]]$C2 <- z
    CM_table1 <- rbind(CM_table1,results[[z]])
  }
  CM_table1 <- CM_table1[-1,]
  write.table(CM_table1,"CM_table1_COM8.txt",sep = "\t",
              row.names = F)
  rm(results,CM_T, CM_table1)
  gc()
  
}

registerDoSEQ() 


# SVM Cross Fold Validation
data <- read.delim("MLI.txt",header = TRUE)
data <- data[,-16]
markers <- as.vector(colnames(data))

ppanel <- list()
count <- 1
for (l in seq(3,8)) {
  comb_opt <- combn(markers, l, simplify = FALSE)
  for (m in seq(length(comb_opt))) {
    ppanel[count] <- data.frame(comb_opt[[m]])
    count <- count + 1
  }
}

metData <- read.delim("mtData.txt", header = T)
ColData <- read.delim("IBDColDataoutrem.txt",header = TRUE)
Class <- ColData$gr
metData <- cbind(metData, Class)
data <- cbind(data, Class)

foreach(M = ppanel) %dopar% {
  library(caret)
  library(MASS)
  library(stringi)
  tdata2 <- data[,c(M,"Class")]
  AllPred <- c()
  ActualClass <- c()
  FoldNumber <- c()
  CoefList <- list()
  SampleID <- c()
  PredClass <- c()
  for (i in 1:5) { 
    set.seed(1234567)
    y_train = metData[metData$Fold != i,]
    y_test = metData[metData$Fold == i,]
    data_train <- tdata2[as.numeric(row.names(y_train)),]
    data_test <- tdata2[as.numeric(row.names(y_test)),]
    y_train$Class <- as.factor(y_train$Class)
    y_test$Class <- as.factor(y_test$Class)
    fit <- svm(y_train$Class ~ ., data = data_train)
    pred <- predict(fit, newdata = data_test[,-ncol(data_test)], type="response")
    #Coefs <- coefficients(fit)
    #CoefList[[i]] <- Coefs
    AllPred <- rbind(AllPred,data.frame(pred))
    SampleID <- rownames(AllPred)
    ActualClass <- c(ActualClass, as.character(data_test$Class))
    #FoldNumber <- c(FoldNumber, rep(i,nrow(pred$posterior)))
  }
  
  All_preds_complete <- data.frame(pred = AllPred, ActualClass = ActualClass, row.names = SampleID)
  
  confResult <- confusionMatrix(data = as.factor(All_preds_complete$pred),
                                reference = as.factor(All_preds_complete$ActualClass), positive = "IBD")
  CMtable <- data.frame(
    "Sensitivity" = round(confResult$byClass[1],3),
    "Specificity" = round(confResult$byClass[2],3),
    "Accuracy" = round(confResult$overall[1],3),
    "Features_Number"= length(M),
    "Panel" = toString(M))
} -> CMT

Ttable <- data.frame()
for (z in seq(length(CMT))) {
  Ttable <- rbind(Ttable,CMT[[z]])
}
write.table(Ttable, "SVM_Combination.txt", sep = "\t")

# ROC curve for best gene sets
library("pROC")
library("caret")
library("MASS")
library("stringi")
data <- read.delim("MLI.txt",header = TRUE)
markers <- as.vector(colnames(data[,-16]))
IBDColDataOutRem <- read.delim(file = "IBDColDataOutRem.txt", header = T)
data <- data[IBDColDataOutRem$accision,]
gr <- IBDColDataOutRem$gr
data <- cbind(data, gr)

ppanel <- list()
count <- 1
for (l in seq(3,8)) {
  comb_opt <- combn(markers, l, simplify = FALSE)
  for (m in seq(length(comb_opt))) {
    ppanel[count] <- data.frame(comb_opt[[m]])
    count <- count + 1
  }
}

ROC_table <- data.frame("response"= "a","predictor"= 0, "fold" = 0,
                        "NO.features"= 0,
                        "method" = "a")

for (i in 1:5) {
  set.seed(5783125)
  y_train = data[data$fold != i,]
  y_test = data[data$fold == i,]
  data_train <- data[row.names(y_train),]
  data_test <- data[row.names(y_test),]
  cl <- makeCluster(detectCores(), type ='PSOCK')
  registerDoParallel(cl)
  set.seed(100)
  y_train$gr <- as.factor(y_train$gr)
  y_test$gr <- as.factor(y_test$gr)
  
  #data_train1 <- data_train[,ppanel[[1]]]
  data_train1 <- data_train[,ppanel[[2]]]
  data_train2 <- data_train[,ppanel[[593]]]
  data_train3 <- data_train[,ppanel[[2065]]]
  data_train4 <- data_train[,ppanel[[5002]]]
  data_train5 <- data_train[,ppanel[[10272]]]
  data_train6 <- data_train[,ppanel[[16648]]]
  
  #data_test1 <- data_test[,ppanel[[1]]]
  data_test1 <- data_test[,ppanel[[2]]]
  data_test2 <- data_test[,ppanel[[593]]]
  data_test3 <- data_test[,ppanel[[2065]]]
  data_test4 <- data_test[,ppanel[[5002]]]
  data_test5 <- data_test[,ppanel[[10272]]]
  data_test6 <- data_test[,ppanel[[16648]]]
  
  cl <- makeCluster(detectCores(), type ='PSOCK')
  registerDoParallel(cl)
  
  #SVM
  
  svm.model1 <- svm(y_train$gr ~ ., data = data_train1, probability = TRUE)
  svm.model2 <- svm(y_train$gr ~ ., data = data_train2, probability = TRUE)
  svm.model3 <- svm(y_train$gr ~ ., data = data_train3, probability = TRUE)
  svm.model4 <- svm(y_train$gr ~ ., data = data_train4, probability = TRUE)
  svm.model5 <- svm(y_train$gr ~ ., data = data_train5, probability = TRUE)
  svm.model6 <- svm(y_train$gr ~ ., data = data_train6, probability = TRUE)
  
  svm.pred1 <- predict(svm.model1, data_test1, probability = TRUE)
  svm.pred2 <- predict(svm.model2, data_test2, probability = TRUE)
  svm.pred3 <- predict(svm.model3, data_test3, probability = TRUE)
  svm.pred4 <- predict(svm.model4, data_test4, probability = TRUE)
  svm.pred5 <- predict(svm.model5, data_test5, probability = TRUE)
  svm.pred6 <- predict(svm.model6, data_test6, probability = TRUE)
  
  svm_pr1  <- as.data.frame(attr(svm.pred1 ,"probabilities"))
  svm_pr2  <- as.data.frame(attr(svm.pred2 ,"probabilities"))
  svm_pr3  <- as.data.frame(attr(svm.pred3 ,"probabilities"))
  svm_pr4  <- as.data.frame(attr(svm.pred4 ,"probabilities"))
  svm_pr5  <- as.data.frame(attr(svm.pred5 ,"probabilities"))
  svm_pr6  <- as.data.frame(attr(svm.pred6 ,"probabilities"))
  
  ROC1 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr1$IBD,
                     "fold"= rep(i,length(svm.pred1)), 
                     "NO.features"= rep(3,length(svm.pred1)),
                     "method" = rep("svm",length(svm.pred1)))
  ROC2 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr2$IBD,
                     "fold"= rep(i,length(svm.pred2)), 
                     "NO.features"= rep(4,length(svm.pred2)),
                     "method" = rep("svm",length(svm.pred2)))
  ROC3 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr3$IBD,
                     "fold"= rep(i,length(svm.pred3)), 
                     "NO.features"= rep(5,length(svm.pred3)),
                     "method" = rep("svm",length(svm.pred3)))
  ROC4 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr4$IBD,
                     "fold"= rep(i,length(svm.pred4)), 
                     "NO.features"= rep(6,length(svm.pred4)),
                     "method" = rep("svm",length(svm.pred4)))
  ROC5 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr5$IBD,
                     "fold"= rep(i,length(svm.pred5)), 
                     "NO.features"= rep(7,length(svm.pred5)),
                     "method" = rep("svm",length(svm.pred5)))
  ROC6 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr6$IBD,
                     "fold"= rep(i,length(svm.pred6)), 
                     "NO.features"= rep(8,length(svm.pred6)),
                     "method" = rep("svm",length(svm.pred6)))
  
  ROC_table <- rbind(ROC_table,ROC1)
  ROC_table <- rbind(ROC_table,ROC2)
  ROC_table <- rbind(ROC_table,ROC3)
  ROC_table <- rbind(ROC_table,ROC4)
  ROC_table <- rbind(ROC_table,ROC5)
  ROC_table <- rbind(ROC_table,ROC6)
  
  rm(ROC1,ROC4,ROC7,ROC10,svm.pred1,svm.pred4,svm.pred8,svm.pred9,svm.model1,
     ROC2,ROC5,ROC8,ROC11,svm.pred2,svm.pred5,svm.pred7,svm.pred10,svm.model2,
     ROC3,ROC6,ROC9,svm.pred3,svm.pred6,svm.pred6,svm.pred11,svm.model3,
     svm.model4,svm.model5,svm.model6,svm.model7,svm.model8,svm.model3,
     svm.model9)
  
  gc()

}

ROC_table <- ROC_table[-1,]

ME <- unique(ROC_table$method)
NF <- as.numeric(unique(ROC_table$NO.features))
FO <- as.numeric(unique(ROC_table$fold))

ROC_svm <- ROC_table[ROC_table$method == "svm",]
ROC_lda <- ROC_table[ROC_table$method == "LDA",]
ROC_logit <- ROC_table[ROC_table$method == "logit",]

ROC_svm_1 <- ROC_svm[ROC_svm$NO.features == 3,]
ROC_svm_2 <- ROC_svm[ROC_svm$NO.features == 4,]
ROC_svm_3 <- ROC_svm[ROC_svm$NO.features == 5,]
ROC_svm_4 <- ROC_svm[ROC_svm$NO.features == 6,]
ROC_svm_5 <- ROC_svm[ROC_svm$NO.features == 7,]
ROC_svm_6 <- ROC_svm[ROC_svm$NO.features == 8,]
                     
library(pROC)
roc.svm1 <- roc(response =as.character(ROC_svm_1$response), 
                predictor = as.numeric(ROC_svm_1$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("Con","IBD")))
roc.svm2 <- roc(response =as.character(ROC_svm_2$response), 
                predictor = as.numeric(ROC_svm_2$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("Con","IBD")))
roc.svm3 <- roc(response =as.character(ROC_svm_3$response), 
                predictor = as.numeric(ROC_svm_3$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("Con","IBD")))
roc.svm4 <- roc(response =as.character(ROC_svm_4$response), 
                predictor = as.numeric(ROC_svm_4$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("Con","IBD")))
roc.svm5 <- roc(response =as.character(ROC_svm_5$response), 
                predictor = as.numeric(ROC_svm_5$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("Con","IBD")))
roc.svm6 <- roc(response =as.character(ROC_svm_6$response), 
                predictor = as.numeric(ROC_svm_6$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("Con","IBD")))

LAB <- c(paste("best 3 features: AUC =",round(roc.svm1$auc,2)), 
         paste("best 4 features: AUC =",round(roc.svm2$auc,2)),
         paste("best 5 features: AUC =",round(roc.svm3$auc,2)),
         paste("best 6 features: AUC =",round(roc.svm4$auc,2)),
         paste("best 7 features: AUC =",round(roc.svm3$auc,2)),
         paste("best 8 features: AUC =",round(roc.svm4$auc,2)))

roc.list <- list(roc.svm1,roc.svm2,roc.svm3,roc.svm4,roc.svm5,roc.svm6)

g.list <- ggroc(roc.list, alpha = 0.5, size = 0.85, legacy.axes = TRUE, print.auc = TRUE)
g.list + scale_colour_manual(values = c("red","#DF8B00","blue","green", "yellow",
                                        "purple","pink", "brown", "black"),
                             labels= LAB)+theme_bw()+
  theme(legend.text=element_text(size=19,family = "sans"),
        legend.title = element_text(size=0,family = "sans"),
        axis.title =element_text(size=26,family = "sans"),
        axis.text=element_text(size=12))+
  theme(legend.position = c(0.7, 0.4),
        legend.direction = "vertical")




