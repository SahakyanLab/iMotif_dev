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
    wk.dir = "/Users/ltamon/DPhil/Collaboration/iMotif_Opt"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/Collaboration/iMotif_Opt"
  } else if(whorunsit == "LiezelLinuxDesk"){
    wk.dir = "/home/ltamon/DPhil/Collaboration/iMotif_Opt"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
datafile = paste0(wk.dir, "/1_Optimus/DATA.txt")
#out.dir = paste0(wk.dir, "/1_Optimus/out/Tm")
out.dir = paste0(wk.dir, "/1_Optimus/out/pHt")
# Optimus saves output to current directory
setwd(out.dir)
### OTHER SETTINGS #############################################################
#CHANGE
# Change Tm/pHt in m() and out.dir above
out.name = "iMotif_Hunter"
# 6 elements for MAX=6, 5 elements for MAX=5
k.initial <- c(k1=1, k2=2, k3=3, k4=4, k5=5, k6=6)
# MAX in HunterScore
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
dump.freq = 5e4

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
source(paste0(wk.dir, "/1_Optimus/lib/Optimus.R"))
source(paste0(wk.dir, "/1_Optimus/lib/OptimusSA.R"))
source(paste0(wk.dir, "/1_Optimus/lib/OptimusRE.R"))
source(paste0(wk.dir, "/1_Optimus/lib/TemperatureControlUnit.R"))
### FUNCTION ###################################################################
DATA$HunterScore <- function(seq="CCTCCTCCCT", K=K){
  
  MAX = 6
  pos.base = "C"
  neg.base = "T"
  
  x <- Rle(strsplit(as.character(seq),NULL)[[1]])
  runValue.x <- runValue(x)
  runLength.x <- runLength(x)
  rm(x)
  
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
  return(K.new)
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(opt.type=="SA"){
  
  Optimus(NCPU=replica,
          SEED=seed,
          DATA=DATA,
          K.INITIAL=k.initial,
          rDEF=r,
          mDEF=m,
          uDEF=u,
          OPT.TYPE=opt.type,
          OPTNAME=paste0(opt.type, "_", out.name),
          NUMITER=numiter,                                                   
          CYCLES=cycles,
          DUMP.FREQ=dump.freq,
          LONG=TRUE)
  
} else if(opt.type=="RE") {
  
  Optimus(NCPU=replica,
          SEED=seed,
          DATA=DATA,
          K.INITIAL=k.initial,
          rDEF=r,
          mDEF=m,
          uDEF=u,
          ACCRATIO=accratio,
          OPT.TYPE=opt.type,
          OPTNAME=paste0(opt.type, "_", out.name),
          NUMITER=numiter,                                                   
          CYCLES=cycles,
          DUMP.FREQ=dump.freq,
          LONG=TRUE)

} else {
  stop("opt.type can only be either SA or RE")
}

# rm(list=ls())

