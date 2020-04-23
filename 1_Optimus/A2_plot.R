################################################################################
# Plot performance of G4Hunter coefficients/Optimus-optimised ones
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
target = "pHt" # "Tm" | "pHt"
datafile = paste0(wk.dir, "/1_Optimus/DATA.txt")
out.dir = paste0(wk.dir, "/1_Optimus/out_zeroT_plot")
optimodel.dir = paste0(wk.dir, "/1_Optimus/out_zeroT_opti/", target)
### OTHER SETTINGS #############################################################
neg.base = NULL
out.name = paste0("zeroT_", target)
opti.model = list(Tm=c(tradOpti="SA_iMotif_Hunter_opti_zeroT1_model", 
                       extOpti="SA_iMotif_Hunter_optiExt_zeroT1_model"),
                  pHt=c(tradOpti="SA_iMotif_Hunter_opti_zeroT1_model", 
                        extOpti="SA_iMotif_Hunter_optiExt_zeroT1_model"))
col.v = c(tradG4="#F8766D", extG4="#00BFC4", tradOpti="#C77CFF", extOpti="#508aba")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(S4Vectors)
library(reshape)
library(ggplot2)
library(ggpubr)
source(paste0(wk.dir, "/1_Optimus/lib/HunterScore.R"))
source(paste0(wk.dir, "/lib/lmEqn_string.R"))
source(paste0(wk.dir, "/lib/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
dta <- read.table(file=datafile, header=TRUE, stringsAsFactors=FALSE)
len <- length(dta[,1])
pdata <- list()
#---------------------------------------
# Traditional/Extended G4Hunter scores
K.lst <- list(tradG4=1:4, extG4=1:6)

for( x in c("tradG4", "extG4") ){
  pdata[[x]] <- sapply(X=dta$Sequence, FUN=HunterScore, K=K.lst[[x]], pos.base="C", 
                       neg.base=NULL, USE.NAMES=FALSE, simplify=TRUE)
}
#---------------------------------------
# Traditional/Extended Optimus scores
for( x in c("tradOpti", "extOpti") ){
  load(file=paste0(optimodel.dir, "/", opti.model[[target]][x], "_K.Rdata"))
  pdata[[x]] <- sapply(X=dta$Sequence, FUN=HunterScore, K=K.stored, pos.base="C", 
                       neg.base=NULL, USE.NAMES=FALSE, simplify=TRUE)
  rm(K.stored); gc()
}
#---------------------------------------
pdata <- cbind(do.call("cbind.data.frame", pdata), target=dta[[target]])
rm(dta); gc()

# Plot
p.lst <- list()
ylb <- list(Tm=expression(bold("T"["m"])), pHt=expression(bold("pH"["T"])))
target.min <- min(pdata$target)
scoring.v <- setdiff(names(pdata),"target")
for(scoring in scoring.v){
  eval(parse(text=paste0(
    "p.lst[[scoring]] <- ggplot(data=pdata, aes(x=", scoring, ", y=target))"
  )))
  p.lst[[scoring]] <- p.lst[[scoring]] +
    geom_point(colour=col.v[scoring], shape=1, size=3, stroke=0.6) +
    labs(title=paste0(target, "_", scoring), 
         y=ylb[[target]], 
         x=expression(bold("Score")), 
         colour=NULL
    ) +
    annotate(geom="text", x=max(pdata[[scoring]]), y=target.min,
             size=4, parse=TRUE, hjust=1, vjust=0,
             label=lmEqn_string(x=scoring, y="target", data=pdata[,c(scoring, "target")])
    ) + 
    bgr1 + 
    theme(legend.text=element_text(size=20, face="bold"),
          legend.title=element_text(size=25, face="bold"),
          plot.title=element_text(size=15, face="bold"))
  print(paste0(scoring, " done!"), quote=FALSE)
}
rm(pdata)
p.lst <- ggarrange(plotlist=p.lst, nrow=2, ncol=2, 
                   legend=NULL, common.legend=FALSE)
ggexport(p.lst, width=12, height=9.6,
         filename=paste0(out.dir, "/", out.name, "_plots.pdf"))
# rm(list=ls())


