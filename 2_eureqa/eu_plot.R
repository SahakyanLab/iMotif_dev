################################################################################
# Generate plots for the models; use test data on retrained model
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/Collaboration/iMotif_Opt"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/Collaboration/iMotif_Opt"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
### OTHER SETTINGS #############################################################
train.FILEPATH  = paste0(wk.dir, "/out_featureTable/DATA_feat.csv")
out.name = "iMotif_pHt_B_CT2length"
#out.name = "iMotif_pHt_B_CT1T2T3"
features.v = c("length", "C", "T1", "T2", "T3")
target = "pHt"
CenterScaling = FALSE
trainDataFr = 0.8
seed = 825 

out.dir = paste0(wk.dir, "/2_eureqa/out/", target)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/lmEqn_string.R"))

# Eureqa sample model (without scaling)
# Tm
eumodel <- function(data=data){
  
  len <- length(data[,1])
  
  outcome <- apply(X=data, MARGIN=1, FUN=function(x){
    
    val1 <- 137-(x["T2"]*x["T3"])+x["T1"]
    return( 102-x["T3"]-(val1/x["C"]) )
    
  }) 
  return(outcome)
}

# pHt all features except length
eumodel <- function(data=data){
  
  len <- length(data[,1])
  
  outcome <- apply(X=data, MARGIN=1, FUN=function(x){
    
    val1 <- 7.37+( 0.00469*(x["T1"])*(x["T2"])*(x["T3"]) )
    val2 <- 0.0272*(x["T3"])
    val3 <- ( 3.2+(0.124*x["T1"]) )/x["C"]
    return(val1-val2-val3)
    
  }) 
  return(outcome)
}

# pHt
eumodel <- function(data=data){
  
  len <- length(data[,1])
  
  outcome <- apply(X=data, MARGIN=1, FUN=function(x){
    
    val1 <- 7.37-(3.69/x["C"])
    val2 <- (0.00549*x["length"])/x["T2"]
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
#-------------------------------------------------------------------------------
trainPred <- eumodel(data=trainData)
testPred <- eumodel(data=testData)

PRED <- rbind(
  cbind.data.frame(obs=trainTarget, pred=trainPred, dataType="Training"),
  cbind.data.frame(obs=testTarget, pred=testPred, dataType="Test")
)
PRED$dataType <- factor(as.factor(PRED$dataType), levels=c("Test", "Training"))

p.param <- NULL
p.param$Tm <- list(c(40,85), 74, c(42,40), 45, 85)
p.param$pHt <- list(c(5.5,7), 6.6, c(5.6,5.5), 5.7, 7)

ggplot(data=PRED, aes(x=pred, y=obs)) +
  geom_abline(intercept=0, slope=1, colour="gray70") + 
  geom_point(aes(colour=PRED$dataType)) + 
  scale_x_continuous(limits=p.param[[target]][[1]]) + 
  scale_y_continuous(limits=p.param[[target]][[1]]) + 
  scale_colour_manual(values=c("#7CAE00", "#425c03")) +
  labs(title=out.name, 
       y=expression(bold("Predicted")), 
       x=expression(bold("Observed")), 
       colour=NULL
  ) +
  annotate(geom="text", x=p.param[[target]][[2]], y=p.param[[target]][[3]], 
           size=4, parse=TRUE,
           label=lmEqn_string(x="obs", y="pred", 
                              data=PRED[PRED$dataType=="Test",c("obs", "pred")])
  ) + 
  annotate(geom="text", x=p.param[[target]][[4]], y=p.param[[target]][[5]],
           size=6, label="Eureqa", fontface=2
  ) + 
  bgr1 +
  theme(legend.position="top")

ggsave(filename=paste0(out.dir, "/", out.name, "_OvsP_plot.pdf"), 
       units="in", width=5.5, height=5.5)

# rm(list=ls())
