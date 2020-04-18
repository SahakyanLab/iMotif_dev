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
out0.dir = paste0(wk.dir, "/1_Optimus/out_zeroT_opti")
### OTHER SETTINGS #############################################################
out.name = "SA_iMotif_zeroT_model"
model = list(Tm=c(model.trad="SA_iMotif_Hunter_opti_zeroT1_model", 
                  model.ext="SA_iMotif_Hunter_optiExt_zeroT1_model"),
             pHt=c(model.trad="SA_iMotif_Hunter_opti_zeroT1_model", 
                  model.ext="SA_iMotif_Hunter_optiExt_zeroT1_model")
             )
col.v = c(Traditional="#C77CFF", Extended="#7CAE00")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(S4Vectors)
library(reshape)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/lmEqn_string.R"))

HunterScore <- function(seq="CCCTCCCTCCCTCCC", K=K, 
                        pos.base = "C", neg.base = NULL){
  
  MAX <- length(K)
  x <- Rle(strsplit(as.character(seq),NULL)[[1]])
  runValue.x <- runValue(x)
  runLength.x <- runLength(x)
  rm(x)
  
  xres <- as.character(runValue.x) 
  
  #-------------
  
  for(i in 1:(MAX-1)){
    xres[runValue.x==pos.base & runLength.x==i]  <-  K[i]  
    if(!is.null(neg.base)){
      xres[runValue.x==neg.base & runLength.x==i]  <- -K[i]
    }
  }
  
  #-------------
  
  xres[runValue.x==pos.base & runLength.x>=MAX]  <-  K[MAX]
  
  #------------- 
  
  if(!is.null(neg.base)){
    xres[runValue.x!=pos.base & runValue.x!=neg.base] <- 0
    xres[runValue.x==neg.base & runLength.x>=MAX]  <- -K[MAX]
  } else {
    xres[runValue.x!=pos.base] <- 0
  }
  
  #-------------
  
  return( sum( (as.numeric(xres)*runLength.x) )/sum(runLength.x) )  
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
dta <- read.table(paste0(iMotif.dir, "/DATA.txt"), header=TRUE, as.is=TRUE)

p.lst <- list()
p.param <- list()
p.param$Tm <- c(42, 40)
p.param$pHt <- c(5.95, 5.9)
p.param$gen <- list(
  bgr1, 
  theme(legend.text=element_text(size=20, face="bold"),
        legend.title=element_text(size=25, face="bold"),
        plot.title=element_text(size=15, face="bold"))
)

ylb <- list(expression(bold("T"["m"])), expression(bold("pH"["t"])))
names(ylb) <- c("Tm", "pHt")

for(target in c("pHt", "Tm")){
  
  out.dir <- paste0(out0.dir, "/", target)
  
  # Traditional
  load(file=paste0(out.dir, "/", model[[target]]["model.trad"], "_K.Rdata"))
  trad <- sapply(X=dta$Sequence, FUN=HunterScore, K=K.stored, pos.base="C", 
                 neg.base=NULL, USE.NAMES=FALSE, simplify=TRUE)
  
  # Extended
  load(file=paste0(out.dir, "/", model[[target]]["model.ext"], "_K.Rdata"))
  ext <- sapply(X=dta$Sequence, FUN=HunterScore, K=K.stored, pos.base="C", 
                neg.base=NULL, USE.NAMES=FALSE, simplify=TRUE)
  
  dta.sub <- cbind.data.frame(trad, ext, dta[,target])
  colnames(dta.sub) <- c("Traditional", "Extended", target)
  
  rm(K.stored); gc()
  
  #---------------------------------------
  # Scatter plot
  
  for(scoring in c("Traditional", "Extended") ){
    
    p.lst[[paste0(target, "_", scoring)]] <- ggplot(data=dta.sub,
                                                    aes_string(x=scoring, y=target)) +
      geom_point(colour=col.v[scoring]) +
      labs(title=paste0(target, "_Optimus_", scoring), 
           y=ylb[[target]], 
           x=expression(bold("Score")), 
           colour=NULL
      ) +
      annotate(geom="text", 
               x=mean(dta.sub[,scoring]),
               y=p.param[[target]],
               size=4, parse=TRUE,
               label=lmEqn_string(x=scoring, y=target, data=dta.sub)
      ) +
      p.param$gen
    
  }
  
}

p.lst <- ggarrange(plotlist=p.lst, nrow=2, ncol=2, 
                   legend=NULL, common.legend=FALSE)
ggexport(p.lst, width=12, height=9.6,
         filename=paste0(out0.dir, "/G4Hunter_opti_", out.name, "_plots.pdf"))

# rm(list=ls())


