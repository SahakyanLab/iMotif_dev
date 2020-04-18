################################################################################
# G4Hunter scoring applied to i-motifs
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/Collaboration/iMotif_Opt/1_Optimus"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
datafile = paste0(wk.dir, "/DATA.txt")
#out.dir = paste0(wk.dir, "/1_Optimus/out/Tm")
out.dir = paste0(wk.dir, "/out_zeroT_G4Hunter")
# Optimus saves output to current directory
### OTHER SETTINGS #############################################################
# Change K and MAX in the MAIN CODE, if necessary. 
neg.base = "T" #NULL
out.name = "negT" #"zeroT"
col.v = c(Traditional="#F8766D", Extended="#00BFC4")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(S4Vectors)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/lmEqn_string.R"))
### FUNCTION ###################################################################
HunterScore <- function(seq="CCTCCTCCCT", K=K, 
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
dta <- read.table(file=datafile, header=TRUE, as.is=TRUE)
len <- length(dta[,1])
trad <- ext <- rep(NA, times=len)

# Traditional G4Hunter score
K = 1:4
for(n in 1:len){
  trad[n] <- HunterScore(seq=dta[n,"Sequence"], K = K,
                         pos.base = "C", neg.base = neg.base)
}

# Extended G4Hunter score
K = 1:6
for(n in 1:len){
  ext[n] <- HunterScore(seq=dta[n,"Sequence"], K = K, 
                        pos.base = "C", neg.base = neg.base)
}

dta <- cbind(Traditional=trad, Extended=ext, dta[,c("pHt", "Tm")])

#---------------------------------------
# Scatter plot

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
  for(scoring in c("Traditional", "Extended") ){
    
    dta.sub <- cbind(dta[, c(scoring, target)])
    p.lst[[paste0(target, "_", scoring)]] <- ggplot(data=dta.sub,
                               aes_string(x=scoring, y=target)) +
      geom_point(colour=col.v[scoring]) +
      labs(title=paste0(target, "_", scoring), 
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
         filename=paste0(out.dir, "/G4Hunter_", out.name, "_plots.pdf"))
  
# rm(list=ls())

