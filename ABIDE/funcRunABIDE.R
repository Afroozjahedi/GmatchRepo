library(partykit)
setwd(
  "C:/Users/Afrooz/Google Drive/education/Computational Science/matching/GmatchRepo/ABIDE/"
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

abideData = read.csv("ABIDE_I_II.csv", header = T)
SelAbideData =abideData[,c("SUB_ID","SITE_ID","ABIDE","GROUP","RMSD","FIQ","AGE" )]
SelAbideData =subset(SelAbideData,RMSD<=0.16 & FIQ>0)
table(complete.cases(SelAbideData))
set.seed(13)
training <- ASD[sample(1:nrow(ASD),130),]
validation <- SelAbideData[training,]

form <- GROUP ~ RMSD + AGE + FIQ 

debug(Gmatch)
Gmatch (ABIDE, form, 10,"opt-coars-exact-rev","chi")
ASD <- subset(SelAbideData,GROUP == 1)
