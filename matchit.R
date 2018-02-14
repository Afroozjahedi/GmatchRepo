 install.packages("MatchIt")
 library(MatchIt)
 update.packages()
 data(lalonde)
 m.out1 <- matchit(treat ~ re74 + re75 + age + educ, data = lalonde,
                    method = "nearest", distance = "logit")
 summary(m.out1)
 m.data1 <- match.data(m.out1)
 chisq.test(m.data1$treat, m.data1$black)
 chisq.test(m.data1$treat, m.data1$hispanic)
 chisq.test(m.data1$treat, m.data1$married)
 pvalre75 <-(t.test(re75 ~ treat, data = m.data1))$p.value
 pvalage <- (t.test(age ~ treat, data = m.data1))$p.value
 pvalre78 <- (t.test(re78 ~ treat, data = m.data1))$p.value
 m.out1 <- matchit(treat ~ educ + black + hispan, data = lalonde,
                  method = "exact")
 m.out <- matchit(treat ~ re74 + re75 + age + educ, data = lalonde, 
                  method = "optimal", ratio = 2)
 #====
 library(optmatch)
 library(RItools)
 data(nuclearplants)
 nuke.nopt <- subset(nuclearplants, pt == 0)
 summary(fullmatch(pr ~ date , data = nuke.nopt))
 summary(fullmatch(pr ~ date, data = nuke.nopt, min = 1))
 summary(fullmatch(pr ~ date , data = nuke.nopt, min = 2, max = 5))
 
 #===== My data ====
 source("smd.R")
# setwd(
 #  "C:/Users/Afrooz/Google Drive/education/Computational Science/matching/new data/proximity matching/R/Gmatch"
 #)
 setwd(
         "C:/Users/Afrooz/Google Drive/education/Computational Science/matching/GmatchRepo"
 )
 source("smd.R")
 library(MatchIt)
 nonOutlierLowADOS = read.csv("lowADOSNonOutlierBefore.csv", header = T)
 nonOutlierLowADOS = (nonOutlierLowADOS)[-11]
 m.out2 = matchit(group ~ Gender + Handedness + RMSD.PRE.censoring +
     Age +  WASI.NVIQ, data = nonOutlierLowADOS,   method = "nearest" )
 m.data1 <- match.data(m.out2)
 pvalGender <- chisq.test(m.data1$group, m.data1$Gender)$p.value
 pvalHandedness <- chisq.test(m.data1$group, m.data1$Handedness)$p.value
 pvalMotion <- (t.test(RMSD.PRE.censoring ~  group, data = m.data1))$p.value
 pvalAge <- (t.test(Age ~  group, data = m.data1))$p.value
 pvalNVIQ <- (t.test(WASI.NVIQ ~  group, data = m.data1))$p.value
 pvals <-
         matrix(c(
                 pvalMotion,
                 pvalAge,
                 pvalNVIQ,
                 pvalGender,
                 pvalHandedness
         ),
         ncol = 5)
 colnames(pvals) <-
         c("RMSD.PRE.censoring",
           "Age",
           "WASI.NVIQ",
           "Gender",
           "Handedness"
           )
 rownames(pvals) <- "Pvals"
 pvals <- as.table(pvals)
 pvals
 debug(smd)
 smdMotion <- smd(m.data1, RMSD.PRE.censoring)
 smdAge <- smd(m.data1, "Age")
 smdNVIQ <- smd(m.data1, "WASI.NVIQ")
 smdGender <- smd(m.data1, "Gender")
 smdHandedness <- smd(m.data1, "Handedness") 
 # Tabling smd for output
 smd <-
         matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
                ncol = 5)
 colnames(smd) <-
         c("RMSD.PRE.censoring",
           "Age",
           "WASI.NVIQ",
           "Gender",
           "Handedness"
           )
 rownames(smd) <- "SMD"
 smd <- as.table(smd)
 smd
 table(m.data1$group)
 plot(m.out2)
 debug(matchit)
 m.out3 <- matchit(group ~ Gender + Handedness + RMSD.PRE.censoring +
                           Age +  WASI.NVIQ, data = nonOutlierLowADOS,  
                   method = "optimal", ratio = 1 )
 ################