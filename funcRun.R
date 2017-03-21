library(partykit)
setwd(
  "C:/Users/Afrooz/Google Drive/education/Computational Science/matching/new data/proximity matching/R/Gmatch"
)
source("x2Distance.R")
source("grow.R")
source("growtree.R")
source("partition.R")
source("ordinalize.R")
source("splitrule.R")
source("printpval.R")
source("rfConst.R")
source("smd.R")
#source("Gmatch.R")
source("x2Dist.R")
#BRIEF = read.csv("BRIEF.csv", header = T)
nonOutlierLowADOS = read.csv("lowADOSNonOutlierBefore.csv", header = T)
#nonOutlierLowADOSAftMat = read.csv("After Matching.csv", header = T)
#nonOutlierLowADOS = read.csv("MMC.csv", header = T)
#nonOutlierLowADOS = read.csv("AGE.csv", header = T)
#data = read.csv("exclude low ADOS.csv", header = T)
varList <- c("group",
             "Gender",
             "Handedness" ,
             "RMSD.PRE.censoring" ,
             "Age",
             "WASI.NVIQ")
#Gmatch (BRIEF,varList,30)


#debug(x2Dist)
debug(Gmatch)
#Gmatch (nonOutlierLowADOS, varList, 1000)#to get propensity score
Gmatch (nonOutlierLowADOS, varList, 10,"inverse")
