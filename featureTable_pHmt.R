################################################################################
# Make table of features of iMotif Data for eureqa and GBM
# Mac, R/3.6.1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/Collaboration/iMotif_Opt"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
datafile = paste0(wk.dir, "/DATA.txt")
#datafile = paste0(wk.dir, "/Data_i-motif_pHt_Tm.csv")                  
out.dir = paste0(wk.dir, "/out_featureTable")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(S4Vectors)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
src.df <- read.table(file=datafile, header=TRUE, as.is=TRUE)
#src.df <- read.csv(file=datafile, header=TRUE, as.is=TRUE)

##which(grepl(colnames(src.df), pattern="Name")) = 2 10 18 27
#src.df <- cbind(src.df[,1:8], src.df[,11:16], src.df[,19:25], src.df[,26:49])

seq.Rle <- sapply(X=src.df$Sequence, simplify=FALSE, FUN=function(seq){
  return(
    Rle(strsplit(as.character(seq),NULL)[[1]])
         )
})

#-------------------
# Identify maximum number of runs
numRuns.v <- sapply(X=src.df$Sequence, simplify=TRUE, FUN=function(seq){
  length( seq.Rle[[seq]]@lengths )
})
# All sequences have run=7
numRuns.v <- unique(numRuns.v)
#-------------------

#-------------------
#  Confirmed that the order of C-T is the same
orderCT <- c("C", "T", "C", "T", "C", "T", "C")
for(seq in src.df$Sequence){
  TF <- identical(seq.Rle[[seq]]@values, orderCT)
  if(TF!=TRUE){
    stop("Checkpoint")
  }
}
#-------------------

#-------------------
# Build table
feat.df <- sapply(X=src.df$Sequence, simplify=FALSE, FUN=function(seq){
  seq.Rle[[seq]]@lengths 
})
feat.df <- data.frame(do.call("rbind", feat.df), row.names=NULL)
colnames(feat.df) <- c("C", "T1", "C2", "T2", "C3", "T3", "C4")
feat.df <- cbind(length=nchar(src.df$Sequence), feat.df)
feat.df <- cbind(feat.df, Tm=src.df$Tm, pHt=src.df$pHt)
#-------------------

# Max length of run is 6
max( unlist(feat.df[,-c(1,9)], use.names=FALSE) )

#-------------------
# C lengths identical
identical(feat.df$C1, feat.df$C2, feat.df$C3, feat.df$C4)
# TRUE

# Since C(C1)=C2=C3=C4, denote length as C
feat.df <- feat.df[, -(which(colnames(feat.df)%in%c("C2", "C3", "C4")))]

#-------------------

write.csv(feat.df, file=paste0(out.dir, "/DATA_feat.csv"),
          row.names=FALSE)

# rm(list=ls())
