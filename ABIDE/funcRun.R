library(partykit)
setwd(
  "C:/Users/Afrooz/Google Drive/education/Computational Science/matching/GmatchRepo"
)
source("grow.R")
source("growtree.R")
source("partition.R")
source("ordinalize.R")
source("splitrule.R")
source("rfConst.R")
source("smd.R")
source("x2Dist.R")

#source("x2Distance.R")
#source("0-1Distance.R")
#source("Gmatch.R")
source("optMatch.R")
#source("oneToneMax.R")
#source("reduceVar.R")
source("iterOpt.R")
#BRIEF = read.csv("BRIEF.csv", header = T)
nonOutlierLowADOS = read.csv("lowADOSNonOutlierBefore.csv", header = T)
#data = read.csv("data4Algorith.csv",header = T)
#nonOutlierLowADOSAftMat = read.csv("After Matching.csv", header = T)
#nonOutlierLowADOS = read.csv("MMC.csv", header = T)
#nonOutlierLowADOS = read.csv("AGE.csv", header = T)
#data = read.csv("exclude low ADOS.csv", header = T)
form <-group ~ RMSD.PRE.censoring + Age + WASI.NVIQ +Gender + Handedness 
# varList <- c("Gender",
#              "Handedness" ,
#              "RMSD.PRE.censoring" ,
#              "Age",
#              "WASI.NVIQ")
#Gmatch (BRIEF,varList,30)


#debug(x2Dist)
#debug(rfConst)
debug(Gmatch)
#debug(x2Dist)
#Gmatch (nonOutlierLowADOS, varList, 1000)#to get propensity score
Gmatch (nonOutlierLowADOS, form, 10,"opt-coars-exact-rev","chi")
#Gmatch (data, form, 1000,"1To3GH")

