################################################################################
# Compare G4Hunter coefficients with Optimus-optimised ones
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
iMotif.dir = wk.dir
out.dir = paste0(wk.dir, "/1_Optimus/out/pHt")
# Optimus saves output to current directory
setwd(out.dir)
### OTHER SETTINGS #############################################################
out.name = "SA_iMotif_Hunter1_model"
k.initial = c(k1=1, k2=2, k3=3, k4=4, k5=5, k6=6)
col.v = c(G4Hunter="#F8766D", Optimised="#00BFC4")
target = "pHt" # "Tm" | "pHt" 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(S4Vectors)
library(reshape)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/lmEqn_string.R"))

HunterScore <- function(seq="CCTCCTCCCT", K=K){
  
  MAX = 6
  pos.base = "C"
  neg.base = "T"
  
  x <- Rle(strsplit(as.character(seq),NULL)[[1]])
  runValue.x <- runValue(x)
  runLength.x <- runLength(x)
  
  xres <- as.character(runValue.x) 
  
  #-------------
  for(i in 1:(MAX-1)){
    xres[runValue.x==pos.base & runLength.x==i]  <-  K[i]  
    xres[runValue.x==neg.base & runLength.x==i]  <- -K[i]  
  }
  #-------------
  xres[runValue.x==pos.base & runLength.x>=MAX]  <-  K[MAX]  
  xres[runValue.x==neg.base & runLength.x>=MAX]  <- -K[MAX]
  #-------------  
  xres[runValue.x!=pos.base & runValue.x!=neg.base] <- 0
  #-------------
  
  return( sum( (as.numeric(xres)*runLength.x) )/sum(runLength.x) )     
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(file=paste0(out.dir, "/", out.name, "_K.Rdata"))
data <- read.table(paste0(iMotif.dir, "/DATA.txt"), header=TRUE, as.is=TRUE)

scores <- list()
scores[["G4Hunter"]] <- sapply(X=data$Sequence, FUN=HunterScore, K=k.initial,
                               USE.NAMES=FALSE, simplify=TRUE)
scores[["Optimised"]] <- sapply(X=data$Sequence, FUN=HunterScore, K=K.stored,
                                USE.NAMES=FALSE, simplify=TRUE)
scores[["Actual"]] <- data[[target]]
rm(data, K.stored); gc()

scores <- data.frame(do.call(cbind, scores), stringsAsFactors=FALSE)

#---------------------------------------

ylb <- list(expression(bold("T"["m"])), expression(bold("pH"["t"])))
names(ylb) <- c("Tm", "pHt")

p.param <- NULL
p.param$Tm <- list(c(40,85), 14, c(42,40), -8, 85)
p.param$pHt <- list(c(5.5,7), 10, c(5.6,5.5), -7.7, 7)

p.lst <- list()
for( i in c("G4Hunter", "Optimised") ){
  
  ## Scatter plot
  p.lst[[i]] <- ggplot( data=scores, aes_string(x=i, y="Actual") ) +
    geom_point(colour=col.v[i]) +
    scale_y_continuous( limits=p.param[[target]][[1]] ) + 
    labs(title=paste0(i, "_", out.name), 
         y=ylb[[target]], 
         x=expression(bold("Score")), 
         colour=NULL
    ) +
    annotate(geom="text", 
             x=ifelse(i=="G4Hunter", 4, p.param[[target]][[2]]),
             y=p.param[[target]][[3]],
             size=4, parse=TRUE,
             label=lmEqn_string(x=i, y="Actual", data=scores[,c(i, "Actual")])
    ) + 
    annotate(geom="text", 
             x=ifelse( i=="G4Hunter", 0.1, p.param[[target]][[4]] ), 
             y=p.param[[target]][[5]], 
             size=6, 
             label=i, fontface=2
    ) + 
    bgr1 +
    theme(legend.text=element_text(size=20, face="bold"),
          legend.title=element_text(size=25, face="bold"),
          plot.title=element_text(size=15, face="bold"))
  
}

p.lst <- ggarrange(plotlist=p.lst, nrow=1, ncol=2, 
                   legend=NULL, common.legend=FALSE)
ggexport(p.lst, width=10, height=8,
         filename=paste0(out.dir, "/", out.name, "_plot.pdf"))
#---------------------------------------

# rm(list=ls())


