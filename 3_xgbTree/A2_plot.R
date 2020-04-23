################################################################################
# Generate plots measuring performance of models; apply retrained model on test data
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
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
target = "pHt" # "Tm" | "pHt"
model.dir = out.dir = paste0(wk.dir, "/3_xgbTree/", target, "/model_2-noLength/retrain")
### OTHER SETTINGS #############################################################
#out.name = "FIT_GBM_2_trainDataFr1"
out.name = "FIT_GBM_2_trainDataFr0.8"
facet.v = c("eta", "min_child_weight", "subsample")
  group = "max_depth"
  group.lab = "Max Tree Depth"
  x = "nrounds"
  x.lab = "No. of boosting iterations"
  y = "RMSE"
  y.lab = y
train.FILEPATH  = paste0(wk.dir, "/out_featureTable/DATA_feat.csv")
facetPlot = FALSE
obsVSpredPlot = TRUE
varImpPlot = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(caret)))
library(reshape)
library(ggplot2)
library(ggpubr)
library(ggrepel)
source(paste0(wk.dir, "/lib/GG_bgr.R"))
source(paste0(wk.dir, "/lib/lmEqn_string.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
model.nme <- load(file=paste0(model.dir, "/", out.name, ".Rdata"))
eval(parse(text=paste0("FIT <- ", model.nme, "; rm(", model.nme, "); gc() ")
           ))
#-------------------------------------------------------------------------------
if(facetPlot){
  # Metric vs. No. of boosting iterations (number of trees)
  df <- FIT$results[,c(x, y, group, facet.v)]
  df <- melt(df, id=c(x, y, facet.v))
  ## Scatter plot
  ggplot( data=df, aes_string(x=x, y=y), aes(group=value) ) +
    geom_point(size=1.5, aes(colour=factor(value)) ) + 
    geom_line(size=1, aes(colour=factor(value)) ) + 
    labs(title=paste0(out.name, "plot"), 
         x=bquote(bold( .(x.lab) )),
         y=bquote(bold( .(y.lab) )),
         colour=bquote(bold( .(group.lab) ))
    ) +
    bgr1 +
    theme(legend.text=element_text(size=20, face="bold"),
          legend.title=element_text(size=25, face="bold"),
          legend.position="top",
          strip.text.x=element_text(size=12),
          strip.background=element_rect(color="black", fill="gray80", size=1, 
                                        linetype="solid")
    ) +
    facet_wrap(~subsample+eta+min_child_weight, labeller=label_both)
  ggsave(filename=paste0(out.dir, "/", out.name, "_facetplot.pdf"), units="in", width=15, height=20)
}
#-------------------------------------------------------------------------------
# Observed vs. Predicted plot
if(obsVSpredPlot){
  testX <- NULL
  testY <- NULL
  dataType="Training"
  
  # Load FEATQCIMP
  load(file=paste0(model.dir, "/", out.name, "_FEATQCIMP.Rdata"))
  
  if("testingData"%in%names(FEATQCIMP) & !is.null(FEATQCIMP$testingData) ){
    print("Validating on test data.", quote=FALSE)
    # Scale test data using mean and sd of training data
    TF <- colnames(FEATQCIMP$testingData)==target
    testY <- FEATQCIMP$testingData[,TF]
    testX <- FEATQCIMP$testingData[,!TF]; rm(TF)
    testX <- sapply(X=colnames(testX), simplify=FALSE, FUN=function(x){
      (testX[[x]]-as.numeric(FEATQCIMP$PREPROC["MEAN",x]))/as.numeric(FEATQCIMP$PREPROC["SD",x])
      })
    testX <- do.call(cbind.data.frame, testX)
    # Calculate metric on?
    dataType.met <- "Test"
  }
  lim <- list(Tm=c(40,85), pHt=c(5.5,7))
  # extractPrediction to also get outcome from training data
  PRED <- extractPrediction(models=list(FIT), testX=testX, testY=testY)
  # Calculate RMSE for test set
  test.TF <- PRED$dataType=="Test"
  rmse.val <- RMSE(obs=PRED[test.TF,"obs"], pred=PRED[test.TF,"pred"]); rm(test.TF)
  
  PRED$dataType <- factor(as.factor(PRED$dataType), levels=c("Test", "Training"))
  
  ggplot(data=PRED, aes(x=pred, y=obs)) +
    geom_abline(intercept=0, slope=1, colour="gray70") + 
    geom_point(size=3, stroke=0.6, shape=1, aes_string(colour="dataType")) + 
    scale_x_continuous(limits=lim[[target]]) + 
    scale_y_continuous(limits=lim[[target]]) + 
    scale_colour_manual(values=c("#654f8f", "paleturquoise3")) +
    labs(title=paste0(out.name, "_", target, "_RMSE=", format(rmse.val,digits=7)), 
         y=expression(bold("Predicted")), 
         x=expression(bold("Observed")), 
         colour=NULL
    ) +
    annotate(geom="text", size=4, parse=TRUE, x=lim[[target]][2], y=lim[[target]][1], 
             hjust=1, vjust=0,
             label=lmEqn_string(x="obs", y="pred", data=PRED[PRED$dataType==dataType.met,
                                                             c("obs", "pred")])
    ) + 
    #annotate(geom="text", size=6, label="xgbTree", fontface=2, hjust=0, vjust=1,
    #         x=lim[[target]][1], y=lim[[target]][2]
    #) + 
    bgr1 +
    theme(legend.position="top",
          plot.title=element_text(size=8, face="bold"))
  ggsave(filename=paste0(out.dir, "/", out.name, "_OvsP_plot.pdf"), 
         units="in", width=5.5, height=5.5)
}
#-------------------------------------------------------------------------------
# Plot feature/variable importance 

if(varImpPlot){
  imp.df <- varImp(FIT, scale=FALSE)
  imp.df <- imp.df$importance
  # Scale importance scores to maximum value
  imp.df[,1] <- imp.df[,1]/imp.df[1,1]*100
  imp.df <- data.frame(feature=rownames(imp.df), Overall=imp.df[,1])
  imp.df <- imp.df[order(imp.df$Overall, decreasing=FALSE),]
  imp.df$feature <- factor(as.factor(imp.df$feature), levels=as.character(imp.df$feature))
  
  ggplot(data=imp.df, aes(x=feature, y=Overall)
  ) +
    geom_segment(aes(x=feature, xend=feature, y=0, yend=Overall),
                 size=1) +
    geom_point(size=4, color="#C77CFF") + 
    #fill=alpha("black", 0.3), alpha=0.7, shape=21, stroke=2) +
    geom_text_repel(aes(label=format(Overall, digits=7)),
                    #box.padding=0.35, 
                    point.padding=0.5,
                    #segment.color="grey50",
                    direction="y",
                    nudge_x=0.1
    ) + 
    labs(title=out.name, 
         x=expression(bold("Importance")), 
         y=expression(bold("Predictor")) 
    ) +
    bgr1 +
    coord_flip() 
  ggsave(filename=paste0(out.dir, "/", out.name, "_varImp_plot.pdf"), 
         units="in", width=5.5, height=5.5)
}
# rm(list=ls())


