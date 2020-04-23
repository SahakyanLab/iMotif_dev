################################################################################
# Function to calculate G4Hunter-like score for i-motifs
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(S4Vectors)
################################################################################
HunterScore <- function(seq = "CCTCCTCCCT", K=K, 
                        pos.base = "C", neg.base = DATA$neg.base){
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