################################################################################
# Build a machine learning model (xgbTree) predicting the Tm/pHt of a limited 
# sub-universe of CT-based i-motifs.
# Athena, R/3.6.1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/SahakyanLab/iMotif_dev"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/iMotif_dev"
  } else if(whorunsit == "LiezelLinuxDesk"){
    wk.dir = "/home/ltamon/SahakyanLab/iMotif_dev"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
} 
target = "pHt" # "Tm" | "pHt"
data.dir = paste0(wk.dir, "/out_featureTable")
out.dir = paste0(wk.dir, "/3_xgbTree/", target, "/model_2-noLength/retrain")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(compiler)))
suppressWarnings(suppressPackageStartupMessages(library(caret)))
suppressWarnings(suppressPackageStartupMessages(library(xgboost)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(doParallel)))
suppressWarnings(suppressPackageStartupMessages(library(corrplot)))
### OTHER SETTINGS #############################################################
train.FILEPATH = paste0(data.dir, "/DATA_feat.csv")
features.v = c("C", "T1", "T2", "T3")
seed          = 825
# Fraction of data for training
trainDataFr   = 0.8
CenterScaling = TRUE
cpu           = 2
#-------------------
# model no. 
ind      = 2
# K-fold validation
CVnum    = 5
# Repeat of K=fold validation
CVrep    = 3
# pHt - with Length
tuneGrid = expand.grid(nrounds=c(750),
                       # Same as in GBM, typically [3,10]
                       max_depth=c(8),
                       # Learning rate
                       eta=c(0.01),
                       # Minimum number of values assigned to a leaf
                       min_child_weight=c(3),
                       # Same as the subsample of GBM. Denotes the fraction of 
                       # observations to be randomly sampled for each tree, typically [0.5,1]
                       # To avoid trees becoming highly correlated
                       subsample=c(0.6),
                       # A node is split only when the resulting split gives a positive 
                       # reduction in the loss function. Gamma specifies the minimum loss 
                       # reduction required to make a split.
                       # gamma=0 means No regularisation, [0,infinity]
                       gamma=c(0),
                       # Similar to max_features in GBM. Denotes the fraction of columns/features 
                       # to be randomly sampled for each tree.
                       # colsample_bytree=1 means no subsampling
                       colsample_bytree = c(1))
# pHt - no Length
tuneGrid = expand.grid(nrounds=c(1500),
                       # Same as in GBM, typically [3,10]
                       max_depth=c(6),
                       # Learning rate
                       eta=c(0.01),
                       # Minimum number of values assigned to a leaf
                       min_child_weight=c(10),
                       # Same as the subsample of GBM. Denotes the fraction of 
                       # observations to be randomly sampled for each tree, typically [0.5,1]
                       # To avoid trees becoming highly correlated
                       subsample=c(0.6),
                       # A node is split only when the resulting split gives a positive 
                       # reduction in the loss function. Gamma specifies the minimum loss 
                       # reduction required to make a split.
                       # gamma=0 means No regularisation, [0,infinity]
                       gamma=c(0),
                       # Similar to max_features in GBM. Denotes the fraction of columns/features 
                       # to be randomly sampled for each tree.
                       # colsample_bytree=1 means no subsampling
                       colsample_bytree = c(1))
#-------------------
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
fitname <- paste0("FIT_GBM_",ind)
registerDoParallel(makeCluster(cpu))

print("##################################################", quote=FALSE)
print(paste0("NOTE: Model tuning for the data ", train.FILEPATH), quote=FALSE)
print(paste0("NOTE: fitname - ", fitname), quote=FALSE)
print(paste0("NOTE: SEED - ", seed), quote=FALSE)
print(paste0("NOTE: Number of CPUs used  - ",        cpu), quote=FALSE)
print("NOTE: Machine learning method used - GBM (xgbTree)", quote=FALSE)
print(paste0("NOTE: Validation scheme - ",CVnum,"-fold CV"), quote=FALSE)
print(paste0("NOTE:                      repeated ",CVrep," times."), quote=FALSE)
print(paste0("NOTE: Training data fraction - ", trainDataFr), quote=FALSE)
print(paste0("NOTE: CenterScaling - ",CenterScaling), quote=FALSE)
print("##################################################", quote=FALSE)
#-------------------------------------------------------------------------------
data <- read.csv(train.FILEPATH, header=TRUE, stringsAsFactors=FALSE)
#data <- data[, colnames(data)%in%c(features.v, target)]
data.len <- nrow(data)
#-------------------------------------------------------------------------------
# Define seeds for the training
if(!is.null(seed)){
  set.seed(seed) 
} else {
  stop("Set seed!")
}
seeds <- vector(mode="list", length=((CVnum*CVrep)+1))
for(si in 1:(CVnum*CVrep)){
  seeds[[si]] <- sample.int(n=10000, length(tuneGrid[,1]))
}
seeds[[((CVnum*CVrep)+1)]] <- sample.int(10000, 1)
#-------------------------------------------------------------------------------
FEATQCIMP <- list()
if(trainDataFr < 1 & trainDataFr > 0){
  set.seed(seed)
  # Initially shuffle once the whole dataset
  data <- data[sample(x=1:data.len, size=data.len, replace=FALSE),]
  ind <- sample(x=1:data.len, size=data.len*trainDataFr, replace=FALSE)
  FEATQCIMP[["testingData"]] <- data[-ind, c(features.v, target)]
} else if(trainDataFr==1){
  print("Training with all data.", quote=FALSE)
  ind <- 1:data.len
  FEATQCIMP[["testingData"]] <- NULL
}
trainData <- data[ind,features.v]
targetData <- data[ind,target]
# Save split of data to train and test to be used later (for plottting and for eureqa)
SPLIT <- list(trainData=data[ind,], testData=data[-ind,])
save(SPLIT, file=paste0(data.dir, "/split_seed", seed, "_trainDataFr", trainDataFr, ".RData"))
rm(ind, data.len, data); gc()
#-------------------------------------------------------------------------------
if(CenterScaling){
  PREPROC <- rbind(MEAN=colMeans(x=trainData),
                   SD=apply(X=trainData, MARGIN=2, FUN=sd))
  trainData <- sapply(X=colnames(trainData), simplify=FALSE,FUN=function(x){
    (trainData[[x]]-as.numeric(PREPROC["MEAN",x]))/as.numeric(PREPROC["SD",x])
  })
  trainData <- do.call(cbind.data.frame, trainData)
  
  FEATQCIMP[["PREPROC"]] <- PREPROC; rm(PREPROC)
  
  print("Predictors centered and scaled w.r.t. mean and sd.", quote=FALSE)
}
#write.csv(trainData, file=paste0(wk.dir, "/out_featureTable/DATA_feat_scaled.csv"),
#          row.names=FALSE)
#-------------------------------------------------------------------------------
# Assess features
FEATQCIMP[["NZV"]] <- nearZeroVar(trainData, saveMetrics=TRUE)
if(any(FEATQCIMP$NZV$zeroVar)==TRUE){stop("Checkpoint: ZeroVar.")}

FEATQCIMP[["COR.MX"]] <- cor(trainData)
FEATQCIMP[["COR.SIG.TEST"]] <- cor.mtest(trainData, conf.level = .95)

pdf(file=paste0(out.dir, "/", fitname, "_trainDataFr", trainDataFr,
                "_corrplot.pdf"), width=10, height=10)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(FEATQCIMP$COR.MX, method="color", col=col(200), order="hclust", mar=c(5, 5, 1.8, 2) + 0.1,
         number.cex=1, number.font=1, number.digits=3,
         # Add coefficient of correlation
         addCoef.col="black", 
         # Text label color and rotation
         tl.cex=1, tl.col="black", tl.srt=90,
         addgrid.col=TRUE, title=paste0(fitname, "_trainDataFr", trainDataFr)
         # Combine with significance
         #p.mat = sig.test$p, sig.level = 0.01, insig = "blank"
         )
dev.off()
rm(col)
#-------------------------------------------------------------------------------
# Doing the architecture tuning across multiple parameters and with 
# CVnum-fold CV repeated CVrep times
print("NOTE: Model fitting...", quote=FALSE)
FIT <- train(x=trainData,
             y=targetData,
             weights = NULL,
             method = "xgbTree",
             metric = "RMSE",
             trControl=trainControl(method = "repeatedcv",
                                    number = CVnum,
                                    repeats = CVrep,
                                    seeds = seeds),
             tuneGrid=tuneGrid,
             verbose=TRUE
)
#-------------------------------------------------------------------------------
# Feature importance
FEATQCIMP[["FeatureIMP"]] <- varImp(FIT, scale=FALSE)
save(FEATQCIMP, file=paste0(out.dir, "/", fitname, "_trainDataFr", trainDataFr,
                        "_FeatQCIMP.RData"))
rm(FEATQCIMP, seeds); gc()
#-------------------------------------------------------------------------------
# Saving the outcome of the architecture tuning
print("NOTE: Saving the results for the first level.", quote=FALSE)
eval(parse(text=paste(fitname," <- FIT",sep="")))
eval(parse(text=paste0(
  "save(",fitname,", file=paste0(out.dir, '/', fitname, '_trainDataFr', trainDataFr, 
  '.RData')); rm(",fitname,"); gc()"
)
))
# rm(list=ls())
