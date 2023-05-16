

install.packages('splitstackshape', repos='http://cran.us.r-project.org')
library(splitstackshape)
library(tidyverse)
library(biomformat) 
library(vegan)
library(glmnet)
library(torch)
library(microDecon)
# install_torch()

#### CODE FROM KNIGHT LAB'S GITUHB ####


# Load dependencies
require(devtools)
require(doMC)
require(tibble)
require(gbm)
require(splitstackshape)
require(reshape2)
require(ggpubr)
require(caret) # for model building
require(pROC) # for AUC calculations
require(purrr) # for functional programming using map()
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(tibble) # for df operations
require(cowplot) # for plotting
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning
require(caret) # for machine learning


sessionInfo()


seed = 2525


#### unused 1  ####

# moved to be closer to function that we use that uses this.
# defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
#                                n.trees = floor((1:3) * 50),
#                                shrinkage = 0.1,
#                                n.minobsinnode = 5)

## unused
# customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
#                               n.trees = floor((1:3) * 50),
#                               shrinkage = 0.1,
#                               n.minobsinnode = 1)

# this values are only used in unused code. The same values are defined internally in the function that is used (loocvDTs).
# numKFold <- 4
# numResampleIter <- 1

#### unused 2 ####

## This function is defined twice but never used.

# ml2DTs <- function(snmData, 
#                    classOfInterest = "Lung Adenocarcinoma", 
#                    cutPoint = 0.5, 
#                    samplingSize = 20, 
#                    caretTuneGrid = defaultGBMGrid){
#   
#   metaTmp1 <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% c("PRAD",
#                                                                                                                  "SKCM",
#                                                                                                                  "NSCLC")),])
#   tmp <- metaTmp1
#   tmp$disease_type_consol <- factor(ifelse(metaTmp1$disease_type_consol == classOfInterest, yes = classOfInterest, no = "Other"))
#   metadataSimSampled <- as.data.frame(stratified(tmp,
#                                                  group = "disease_type_consol",
#                                                  size = samplingSize,
#                                                  keep.rownames = TRUE,
#                                                  replace = FALSE,
#                                                  bothSets = FALSE))
#   rownames(metadataSimSampled) <- metadataSimSampled$rn
#   mlDataY <- metadataSimSampled
#   mlDataX <- snmData[rownames(mlDataY),]
#   
#   set.seed(seed)
#   index <- createDataPartition(mlDataY$disease_type_consol, p = 0.7, list = FALSE)
#   trainX <- mlDataX[index,]
#   trainY <- mlDataY[index,]$disease_type_consol
#   testX <- mlDataX[-index,]
#   testY <- mlDataY[-index,]$disease_type_consol
#   # print(testY)
#   
#   refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
#   refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
#   
#   set.seed(seed)
#   ctrl <- trainControl(method = "repeatedcv",
#                        number = numKFold,
#                        repeats = numResampleIter,
#                        sampling = "up",
#                        summaryFunction = twoClassSummary,
#                        classProbs = TRUE,
#                        verboseIter = TRUE,
#                        savePredictions = TRUE,
#                        allowParallel=TRUE)
#   
#   mlModel <- train(x = trainX,
#                    y = refactoredTrainY,
#                    method = "gbm",
#                    preProcess = c("scale","center"),
#                    trControl = ctrl,
#                    verbose = TRUE,
#                    metric = "ROC",
#                    tuneGrid = customGBMGrid)
#   
#   positiveClass <- gsub(" ","", classOfInterest)
#   negativeClass <- "Other"
#   
#   predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
#   fg <- predProbs[refactoredTestY == positiveClass]
#   bg <- predProbs[refactoredTestY == negativeClass]
#   
#   prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
#   prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
#   
#   # par(mfrow = c(1,2))
#   plot(prroc_roc)
#   plot(prroc_pr)
#   # dev.off()
#   
#   
#   predClass <- predict(mlModel, newdata = testX)
#   
#   confusionMatrix(table(predict(mlModel, newdata = testX, type="prob")[,positiveClass] >= cutPoint,
#                         refactoredTestY == positiveClass))
# }

#-----------------------------------------#
# Machine learning
#-----------------------------------------#


# mlHvsC <- function(snmData){
# Load dependencies


## not used
# 
# mlHvsC <- function(snmData){
#   
#   numCores <- detectCores()
#   registerDoMC(cores=numCores)
#   
#   defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
#                                  n.trees = floor((1:3) * 50),
#                                  shrinkage = 0.1,
#                                  n.minobsinnode = 5)
#   customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
#                                 n.trees = floor((1:3) * 50),
#                                 shrinkage = 0.1,
#                                 n.minobsinnode = 1)
#   
#   caretTuneGrid <- defaultGBMGrid
#   numKFold <- 4
#   numResampleIter <- 1
#   
#   mlDataY <- metadataPSMatchedDPQCFiltered
#   mlDataX <- snmData[rownames(mlDataY),]
#   
#   set.seed(seed)
#   index <- createDataPartition(mlDataY$HvsC, p = 0.7, list = FALSE)
#   trainX <- mlDataX[index,]
#   trainY <- mlDataY[index,]$HvsC
#   testX <- mlDataX[-index,]
#   testY <- mlDataY[-index,]$HvsC
#   
#   refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
#   refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
#   
#   set.seed(seed)
#   ctrl <- trainControl(method = "repeatedcv",
#                        number = numKFold,
#                        repeats = numResampleIter,
#                        sampling = "up",
#                        summaryFunction = twoClassSummary,
#                        classProbs = TRUE,
#                        verboseIter = TRUE,
#                        savePredictions = TRUE,
#                        allowParallel=TRUE)
#   
#   mlModel <- train(x = trainX,
#                    y = refactoredTrainY,
#                    method = "gbm",
#                    preProcess = c("scale","center"),
#                    trControl = ctrl,
#                    verbose = TRUE,
#                    metric = "ROC",
#                    tuneGrid = defaultGBMGrid)
#   
#   positiveClass <- "Cancer"
#   negativeClass <- "Control"
#   predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
#   fg <- predProbs[refactoredTestY == positiveClass]
#   bg <- predProbs[refactoredTestY == negativeClass]
#   
#   prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
#   prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
#   
#   plot(prroc_roc)
#   plot(prroc_pr)
#   
#   predClass <- predict(mlModel, newdata = testX)
#   print(confusionMatrix(data = predClass, reference = refactoredTestY, positive = positiveClass))
# }


#### loocvDTs - the magic ####

metadataPSMatchedDPQCFiltered <- read.csv('../data/Fig4_plasma/Metadata-Plasma-Filtered-For-Analysis.csv', row.names=1)

defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                               n.trees = floor((1:3) * 50),
                               shrinkage = 0.1,
                               n.minobsinnode = 5)

## This loocvDTs function is the main event.

loocvDTs <- function(snmData, samplingSize = 15, DTs, caretTuneGrid = defaultGBMGrid,
                     filenameString = paste(DTs,collapse = "__"), HvsCFlag = FALSE){
  
  if(HvsCFlag){
    metaTmpX <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% DTs),])
    metaTmpX$disease_type_consol <- metaTmpX$HvsC
    classes <- gsub(" ","",levels(metaTmpX$disease_type_consol))
  } else{
    metaTmpX <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% DTs),])
    classes <- gsub(" ","",DTs)
  }
  
  # Do LOOCV model building and testing
  
  multiClassSummaryStats <- list()
  multiClassSummaryStatsDist <- list()
  numKFold <- 4
  numResampleIter <- 1
  metaData <- metaTmpX
  snmData <- snmData # dataPSUniqueDecontamQC # 
  iterSize <- 1
  for(jj in 1:iterSize){
    metadataSimSampled <- as.data.frame(stratified(metaData,
                                                   group = "disease_type_consol",
                                                   size = samplingSize,
                                                   keep.rownames = TRUE,
                                                   replace = FALSE,
                                                   bothSets = FALSE))
    rownames(metadataSimSampled) <- metadataSimSampled$rn
    mlDataY <- metadataSimSampled
    mlDataX <- snmData[rownames(mlDataY),]
    dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
    
    # Create data partitions
    # set.seed(42)
    indexSuper <- 1:dim(mlDataY)[1]
    predProbs <- list()
    obsClass <- vector()
    predClass <- vector()
    varImpBestModelDF2OrderedNonzeroList <- list()
    
    for(ii in 1:length(indexSuper)){
      print(sprintf("Iteration: %d/%d", ii, length(indexSuper)))
      index <- indexSuper[ii]
      trainX <- mlDataX[-index,]
      trainY <- mlDataY[-index,]$disease_type_consol
      testX <- mlDataX[index,,drop=FALSE]
      testY <- mlDataY[index,,drop=FALSE]$disease_type_consol
      
      refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
      refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
      
      obsClass[ii] <- as.character(refactoredTestY)
      
      set.seed(seed)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = multiClassSummary,
                           classProbs = TRUE,
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      mlModel <- train(x = trainX,
                       y = refactoredTrainY,
                       method = "gbm",
                       preProcess = c("scale","center"),
                       trControl = ctrl,
                       verbose = FALSE,
                       metric = "ROC",
                       tuneGrid = caretTuneGrid)
      
      predProbs[ii] <- list(predict(mlModel, newdata = testX, type = "prob"))
      predClass[ii] <- as.character(predict(mlModel, newdata = testX, type = "raw"))
      
      varImpBestModelDF <- as.data.frame(varImp( mlModel$finalModel, scale = FALSE ))
      varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
      varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
      colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
      varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
      varImpBestModelDF2OrderedNonzeroList[[ii]] <- varImpBestModelDF2OrderedNonzero
      
      rm(mlModel)
    }
    
    loocvPreds <- cbind(obs = factor(obsClass,
                                     levels = classes),
                        pred = factor(predClass,
                                      levels = classes),
                        do.call(rbind,predProbs))

    
    multiClassSummaryStats[[jj]] <- multiClassSummary(loocvPreds, lev = classes)
    print(multiClassSummaryStats[[jj]])
    
    loocvPreds %>% write.csv( paste0(filenameString, "__Preds.csv"))
    
    filenameROC <- paste0(filenameString,"__ROC.png")
    filenamePR <- paste0(filenameString,"__PR.png")
    filenameROCData <- paste0(filenameString,"__Data__ROC.csv")
    filenamePRData <- paste0(filenameString,"__Data__PR.csv")
    filenameSink <- paste0(filenameString,"__CM.txt")
    
    predProbs <- loocvPreds[,DTs[1]]
    fg <- predProbs[loocvPreds$obs == DTs[1]]
    bg <- predProbs[loocvPreds$obs == DTs[2]]
    
    prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
    
    png(filename=filenameROC, width = 6, height = 4, units = 'in', res = 300)
    plot(prroc_roc)
    dev.off()
    
    png(filename=filenamePR, width = 6, height = 4, units = 'in', res = 300)
    plot(prroc_pr)
    dev.off()
    
    rocCurveData <- cbind(as.data.frame(prroc_roc$curve), DT1 = DTs[1], DT2 = DTs[2])
    prCurveData <- cbind(as.data.frame(prroc_pr$curve), DT1 = DTs[1], DT2 = DTs[2])
    
    write.table(prCurveData, sep=",", file = filenamePRData, col.names = FALSE)
    write.table(rocCurveData, sep=",", file = filenameROCData, col.names = FALSE)
  }
  
  print(confusionMatrix(loocvPreds$obs, loocvPreds$pred))
  multiClassSummaryStatsDist <- data.frame(do.call(rbind, multiClassSummaryStats))
  
  sink(filenameSink)
  print(print(confusionMatrix(loocvPreds$obs, loocvPreds$pred)))
  sink()
  
  return(multiClassSummaryStats[jj])
}


#### unused 3 ####
## This function is defined twice but never used.
# ml2DTs <- function(snmData, 
#                    classOfInterest = "Lung Adenocarcinoma", 
#                    cutPoint = 0.5, 
#                    samplingSize = 20, 
#                    caretTuneGrid = defaultGBMGrid){
#   
#   metaTmp1 <- droplevels(metadataPSMatchedDPQCFiltered[(metadataPSMatchedDPQCFiltered$disease_type_consol %in% c("PRAD",
#                                                                                                                  "SKCM",
#                                                                                                                  "NSCLC")),])
#   tmp <- metaTmp1
#   tmp$disease_type_consol <- factor(ifelse(metaTmp1$disease_type_consol == classOfInterest, yes = classOfInterest, no = "Other"))
#   metadataSimSampled <- as.data.frame(stratified(tmp,
#                                                  group = "disease_type_consol",
#                                                  size = samplingSize,
#                                                  keep.rownames = TRUE,
#                                                  replace = FALSE,
#                                                  bothSets = FALSE))
#   rownames(metadataSimSampled) <- metadataSimSampled$rn
#   mlDataY <- metadataSimSampled
#   mlDataX <- snmData[rownames(mlDataY),]
#   
#   set.seed(seed)
#   index <- createDataPartition(mlDataY$disease_type_consol, p = 0.7, list = FALSE)
#   trainX <- mlDataX[index,]
#   trainY <- mlDataY[index,]$disease_type_consol
#   testX <- mlDataX[-index,]
#   testY <- mlDataY[-index,]$disease_type_consol
# 
#   
#   refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
#   refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
#   
#   set.seed(seed)
#   ctrl <- trainControl(method = "repeatedcv",
#                        number = numKFold,
#                        repeats = numResampleIter,
#                        sampling = "up",
#                        summaryFunction = twoClassSummary,
#                        classProbs = TRUE,
#                        verboseIter = TRUE,
#                        savePredictions = TRUE,
#                        allowParallel=TRUE)
#   
#   mlModel <- train(x = trainX,
#                    y = refactoredTrainY,
#                    method = "gbm",
#                    preProcess = c("scale","center"),
#                    trControl = ctrl,
#                    verbose = TRUE,
#                    metric = "ROC",
#                    tuneGrid = customGBMGrid)
#   
#   positiveClass <- gsub(" ","", classOfInterest)
#   negativeClass <- "Other"
#   
#   predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
#   fg <- predProbs[refactoredTestY == positiveClass]
#   bg <- predProbs[refactoredTestY == negativeClass]
#   
#   prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
#   prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
#   
#   # par(mfrow = c(1,2))
#   plot(prroc_roc)
#   plot(prroc_pr)
#   # dev.off()
#   
#   
#   predClass <- predict(mlModel, newdata = testX)
#   
#   confusionMatrix(table(predict(mlModel, newdata = testX, type="prob")[,positiveClass] >= cutPoint,
#                         refactoredTestY == positiveClass))
# }

#### Poore et al (unused) ####
## following pipeline implemented by Poore et al

# xy <-  metadataPSMatchedDPQCFiltered %>% 
#   count(disease_type_consol) %>% 
#   arrange(n) %>%
#   pull(disease_type_consol) %>% as.character() %>%
#   combn(2) 

#copy the sampling size's from the knight lab's notebook
# (https://github.com/biocore/tcga/blob/master/jupyter_notebooks/Plasma%20Kraken%20Machine%20Learning%20LOO%20Analysis.ipynb)
# sampling_sizes <- c(25, 59, 69, 59, 69, 69)

# 
# ## loop through the 6 idifferent prediction tasks
# for(idx in 1:ncol(xy)){
#   scrubbed_hVsC_tmp <- loocvDTs(snmData = scrubbed_normalized,
#                                 samplingSize = sampling_sizes[idx], 
#                                 DTs = xy[ ,idx],
#                                 filenameString = paste0('../results/data/Fig4_plasma/trial_', seed, '/scrubbed_preds/', xy[,idx][1], '_', xy[,idx][2] ),
#                                 caretTuneGrid = defaultGBMGrid)
# }
# 
# ## loop through the 6 idifferent prediction tasks
# for(idx in 1:ncol(xy)){
#   raw_hVsC_tmp <- loocvDTs(snmData = raw_normalized,
#                            samplingSize = sampling_sizes[idx], 
#                            DTs = xy[ ,idx],
#                            filenameString = paste0('../results/data/Fig4_plasma/trial_', seed, '/raw_preds/', xy[,idx][1], '_', xy[,idx][2] ),
#                            caretTuneGrid = defaultGBMGrid)
# }
# 
# 
# ## loop through the 6 idifferent prediction tasks
# for(idx in 1:ncol(xy)){
#   raw_hVsC_tmp <- loocvDTs(snmData = dec_normalized,
#                            samplingSize = sampling_sizes[idx],
#                            DTs = xy[ ,idx],
#                            filenameString = paste0('../results/data/Fig4_plasma/trial_', seed, '/dec_preds/', xy[,idx][1], '_', xy[,idx][2] ),
#                            caretTuneGrid = defaultGBMGrid)
# }


#### find data files ####

decontaminationFolder = "../results/data/ivory/decontamination"
infiles = dir(decontaminationFolder, pattern="vsnm.csv", full.names = T, recursive = T)

predictionFolder = "../results/data/ivory/prediction"
suppressWarnings(dir.create(predictionFolder))

#### prediction tasks ####

xy <-  metadataPSMatchedDPQCFiltered %>% 
  count(disease_type_consol) %>% 
  arrange(n) %>%
  pull(disease_type_consol) %>% as.character() %>%
  combn(2) 

# for time, limit this set to just two tasks.
xy = xy[,c(3,5)]
print("Prediction tasks:")
print(xy)

#### read and process data ####

totalCycles = length(infiles) * ncol(xy)
ithCycle = 0

for (infile in infiles){
  message("Training model based on data in file: ", infile)
  data = read.csv(infile, row.names = 1)
  
  outfile = sub(decontaminationFolder, predictionFolder, infile)
  outfile = sub("vsnm.csv", "preds.csv", outfile)
  
  # loop through tasks
  for(idx in 1:ncol(xy)){
    ithCycle = ithCycle + 1
    message("Starting cycle ", ithCycle, " of ", totalCycles, ".") # keep this in inner-most loop
    
    catA = xy[ ,idx][[1]]
    catB = xy[ ,idx][[2]]
    message("Training to distinguish between ", catA, " and ", catB, ".")
    
    fileString = paste0("preds_", catA, "vs", catB)
    outfileString = sub("preds.csv", fileString, outfile)
    
    scrubbed_hVsC_tmp <- loocvDTs(snmData = data,
                                  samplingSize = 69, #sampling_sizes[idx], 
                                  DTs = xy[ ,idx],
                                  filenameString = outfileString,
                                  caretTuneGrid = defaultGBMGrid)
    
    message("Finished cycle ", ithCycle, " of ", totalCycles, ".") # keep this in inner-most loop
  }
  
}




message("Done!)")
