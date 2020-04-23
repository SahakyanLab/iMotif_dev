################################################################################
# Generate plots for the models; use same test data as in GBMs
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/iMotif_dev"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/iMotif_dev"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
} 
target = "pHt"
out.dir = paste0(wk.dir, "/2_eureqa/out/", target)
### OTHER SETTINGS #############################################################
train.FILEPATH = paste0(wk.dir, "/out_featureTable/DATA_feat.csv")
# Tm
out.name = "iMotif_chosen_BE_unscaled"
# pHt
#out.name = "Motif_pHt_BE_unscaled"

# Should be the same as in GBM's
trainDataFr = 0.8
seed = 825 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/lmEqn_string.R"))

eumodel <- list()
# Tm; from ./out/Tm/iMotif_chosen_B.fxp; complexity=14; no scaling of features
eumodel$Tm <- function(data=data){
  outcome <- apply(X=data, MARGIN=1, FUN=function(x){
    val1 <- 137-(x["T2"]*x["T3"])+x["T1"]
    return( 102-x["T3"]-(val1/x["C"]) )
  }) 
  return(outcome)
}
# pHt; from ./out/pHt/iMotif_pHt_BE.fxp; complexity=13; no scaling of features
eumodel$pHt <- function(data=data){
  outcome <- apply(X=data, MARGIN=1, FUN=function(x){
    val1 <- 7.38-(3.70/x["C"])
    val2 <- (0.00565*x["length"])/x["T2"]
    return(val1-val2)
  }) 
  return(outcome)
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
data <- read.csv(train.FILEPATH, header=TRUE, as.is=TRUE)
data.len <- nrow(data)
#-------------------------------------------------------------------------------
# Split into test and train data exactly the same as GBM's
set.seed(seed)
# Initially shuffle once the whole dataset
data <- data[sample(x=1:data.len, size=data.len, replace=FALSE),]
ind <- sample(x=1:data.len, size=data.len*trainDataFr, replace=FALSE)
TF <- colnames(data)==target

trainData <- data[ind,!TF]
trainTarget <- data[ind, TF]
testData <- data[-ind,!TF]
testTarget <- data[-ind, TF]
rm(data, data.len, ind, TF); gc()
# Save so I can compare with GBM to make sure train and test data are identical
DTA <- list(train=trainData, test=testData)
save(DTA, file=paste0(out.dir, "/split_eu_", target, "_seed", seed, "_trainDataFr", trainDataFr, ".RData"))
#-------------------------------------------------------------------------------
trainPred <- eumodel[[target]](data=trainData)
testPred <- eumodel[[target]](data=testData)

PRED <- rbind(
  cbind.data.frame(obs=trainTarget, pred=trainPred, dataType="Training"),
  cbind.data.frame(obs=testTarget, pred=testPred, dataType="Test")
)
PRED$dataType <- factor(as.factor(PRED$dataType), levels=c("Test", "Training"))

#p.param <- NULL
#p.param$Tm <- list(c(40,85), 74, c(42,40), 45, 85)
#p.param$pHt <- list(c(5.5,7), 6.6, c(5.6,5.5), 5.7, 7)

lim <- list(Tm=c(40,85), pHt=c(5.5,7))

ggplot(data=PRED, aes(x=pred, y=obs)) +
  geom_abline(intercept=0, slope=1, colour="gray70") + 
  geom_point(aes_string(colour="dataType")) + 
  scale_x_continuous(limits=lim[[target]]) + 
  scale_y_continuous(limits=lim[[target]]) + 
  scale_colour_manual(values=c("#b3e2e3", "#37aaad")) +
  labs(title=out.name, 
       y=expression(bold("Predicted")), 
       x=expression(bold("Observed")), 
       colour=NULL
  ) +
  annotate(geom="text", x=lim[[target]][2], y=lim[[target]][1], 
           size=4, parse=TRUE, hjust=1, vjust=0,
           label=lmEqn_string(x="obs", y="pred", 
                              data=PRED[PRED$dataType=="Test",c("obs", "pred")])
  ) + 
  annotate(geom="text", x=lim[[target]][1], y=lim[[target]][2],
           size=6, label="Eureqa", fontface=2, hjust=0, vjust=1
  ) + 
  bgr1 +
  theme(legend.position="top")

ggsave(filename=paste0(out.dir, "/", out.name, "_", target, "_OvsP_plot.pdf"), 
       units="in", width=5.5, height=5.5)

# rm(list=ls())
