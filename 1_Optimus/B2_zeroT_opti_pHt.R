################################################################################
# Optimise G4Hunter for i-motif
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/Collaboration/iMotif_Opt/1_Optimus"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/Collaboration/iMotif_Opt/1_Optimus"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
datafile = paste0(wk.dir, "/DATA.txt")
out.dir = paste0(wk.dir, "/out_zeroT_opti/pHt")
# Optimus saves output to current directory
setwd(out.dir)
### OTHER SETTINGS #############################################################
#CHANGE
# Change Tm/pHt in m() and out.dir above
# 6 elements for MAX=6, 5 elements for MAX=5
K.lst <- list(
  `iMotif_Hunter_optiExt_zeroT` = c(k1=1, k2=2, k3=3, k4=4, k5=5, k6=6),
  `iMotif_Hunter_opti_zeroT` = c(k1=1, k2=2, k3=3, k4=4)
)
opt.type = "SA" # "SA" | "RE"
replica = 3 # 3 | # 12

seed = 840
# If opt.type = "RE" 
accratio = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2)
# If opt.type = "RE", replica should be equal to length of ACCRATIO
numiter = 1e6
cycles = 4
# For opt.type="SA", statwindow=70 (default); For opt.type="RE", statwindow=50
statwindow = 50
dump.freq = 1e4 #5e4

DATA <- NULL
DATA$data <- read.table(file=datafile, header=TRUE, as.is=TRUE)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
# BiocManager::install("S4Vectors")
library(S4Vectors)
source(paste0(wk.dir, "/lib/Optimus.R"))
source(paste0(wk.dir, "/lib/OptimusSA.R"))
source(paste0(wk.dir, "/lib/OptimusRE.R"))
source(paste0(wk.dir, "/lib/TemperatureControlUnit.R"))
### FUNCTION ###################################################################
DATA$HunterScore <- function(seq="CCTCCTCCCT", K=K, 
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
# Define interfacing functions for Optimus
m <- function(K, DATA){
  
  scores <- sapply(DATA$data$Sequence, FUN=DATA$HunterScore, K=K,
                   USE.NAMES=FALSE, simplify=TRUE)

  #return( abs(cor(x=DATA$data$Tm, y=scores)) )
  return( abs(cor(x=DATA$data$pHt, y=scores)) )

}

u <- function(O, DATA=NULL){
  result <- NULL
  result$Q <- O
  result$E <- -O
  return(result)
}

r <- function(K){
  
  move.step = 0.1
  
  K.new <- K
  # Randomly selecting a coefficient to alter:
  K.ind.toalter <- sample(size = 1, x = 1:length(K.new))
  # Creating a potentially new set of coefficients where one randomly selected 
  # entry is altered by either +move.step or -move.step:
  K.new[K.ind.toalter] <- K.new[K.ind.toalter] + 
    sample(size = 1, x = c(-move.step, move.step))
  ## Setting the negative coefficients to 0 (not necessary in this example,
  ## useful for optimising rate constants):
  K.new[K.new < 0] <- 0
  K.new <- sort(K.new, decreasing=FALSE)
  return(K.new)
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(K.name in names(K.lst)){
  
  if(opt.type=="SA"){
    
    Optimus(NCPU=replica,
            SEED=seed,
            DATA=DATA,
            K.INITIAL=K.lst[[K.name]],
            rDEF=r,
            mDEF=m,
            uDEF=u,
            OPT.TYPE=opt.type,
            OPTNAME=paste0(opt.type, "_", K.name),
            NUMITER=numiter,                                                   
            CYCLES=cycles,
            DUMP.FREQ=dump.freq,
            LONG=TRUE)
    
  } else if(opt.type=="RE") {
    
    Optimus(NCPU=replica,
            SEED=seed,
            DATA=DATA,
            K.INITIAL=K.lst[[K.name]],
            rDEF=r,
            mDEF=m,
            uDEF=u,
            ACCRATIO=accratio,
            OPT.TYPE=opt.type,
            OPTNAME=paste0(opt.type, "_", K.name),
            NUMITER=numiter,                                                   
            CYCLES=cycles,
            DUMP.FREQ=dump.freq,
            LONG=TRUE)
    
  } else {
    stop("opt.type can only be either SA or RE")
  }
  
}

# rm(list=ls())

