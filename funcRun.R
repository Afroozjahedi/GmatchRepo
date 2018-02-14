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
source("x2Distance.R")
#source("0-1Distance.R")
#source("Gmatch.R")
source("optMatch.R")
source("summaryGmatch.R")
#source("reduceVar.R")
source("iterOpt.R")

#### Data ####
#BRIEF = read.csv("BRIEF.csv", header = T)
nonOutlierLowADOS = read.csv("lowADOSNonOutlierBefore.csv", header = T)
nonOutlierLowADOSNAVar= as.data.frame(lapply(nonOutlierLowADOS[c(-2,-1)],function(x) x[sample(c(TRUE,NA),prob=c(0.85,0.10), size=length(x), replace=TRUE)]))  
nonOutlierLowADOS=data.frame(nonOutlierLowADOS[,c(1,2)],nonOutlierLowADOSNAVar)
#data = read.csv("data4Algorith.csv",header = T)
#nonOutlierLowADOSAftMat = read.csv("After Matching.csv", header = T)
#nonOutlierLowADOS = read.csv("MMC.csv", header = T)
#nonOutlierLowADOS = read.csv("AGE.csv", header = T)
#data = read.csv("exclude low ADOS.csv", header = T)

### Formula ####
form <-group ~ RMSD.PRE.censoring + Age + WASI.NVIQ +Gender + Handedness 
# varList <- c("Gender",
#              "Handedness" ,
#              "RMSD.PRE.censoring" ,
#              "Age",
#              "WASI.NVIQ")
#Gmatch (BRIEF,varList,30)

### Debug function ####
#debug(x2Dist)
#debug(rfConst)
debug(Gmatch)

#Gmatch (nonOutlierLowADOS, form, 1000,"1To3Dist")#1:3
#Gmatch (nonOutlierLowADOS, form, 1000,"propensity" )
#Gmatch (Gmatch (nonOutlierLowADOSNA, form, 1000,"propensity" ), form, 1000,"propensity" )
#Gmatch (nonOutlierLowADOS, form, 1000, "propensity")#to get propensity score
Gmatch (nonOutlierLowADOS, form, 1000,"opt-one-to-one","chi")
##Gmatch (nonOutlierLowADOS, form, 1000,"opt-coars-exact-rev","chi")
#Gmatch (nonOutlierLowADOS, form, 10,"opt-coars-exact-rev")
### Statistics of the dataset ####
mean(subset(optData, group==0)$Age)
sd(subset(optData, group==0)$Age)
range(subset(optData, group==0)$Age)

mean(subset(optData, group==0)$RMSD.PRE.censoring)
sd(subset(optData, group==0)$RMSD.PRE.censoring)
range(subset(optData, group==0)$RMSD.PRE.censoring)

mean(subset(optData, group==0)$WASI.NVIQ)
sd(subset(optData, group==0)$WASI.NVIQ)
range(subset(optData, group==0)$WASI.NVIQ)

table(subset(optData, group==0)$Handedness)
table(subset(dataNewPro, group==1)$Gender)
