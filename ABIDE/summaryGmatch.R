summaryGmatch <- function(DATA){
#=== Tabling pvalues and SMD for output ====
pvalMotion <-
        (t.test(RMSD ~ GROUP, data = DATA))$p.value
meanRMSD <- t.test(RMSD ~ GROUP, data = DATA)$estimate
pvalAge <- (t.test(AGE ~ GROUP, data = DATA))$p.value
pvalFIQ <- (t.test(FIQ ~ GROUP, data = DATA))$p.value
#pvalGender <- chisq.test(DATA$GROUP, DATA$Gender)$p.value
#pvalHandedness <-
#        chisq.test(DATA$GROUP, DATA$Handedness)$p.value

pvalsOpt <-
        matrix(c(
                pvalMotion,
                pvalAge,
                pvalFIQ
        ),
        ncol = 3)
colnames(pvalsOpt) <-
        c("RMSD",
          "Age",
          "FIQ"
        )
rownames(pvalsOpt) <- "Pvals"
pvalsOpt <- as.table(pvalsOpt)
pvalsOpt

meanRMSD<- t.test(RMSD ~ GROUP, data = DATA)$estimate
meanAGE <- t.test(AGE ~ GROUP, data = DATA)$estimate
meanFIQ <- t.test(FIQ ~ GROUP, data = DATA)$estimate
mean <-   matrix(c(
        meanRMSD,
        meanAGE,
        meanFIQ
),
ncol = 3)
colnames(mean) <-
        c("RMSD",
          "Age",
          "FIQ"
        )
rownames(mean) <- c("meanASD","meanTD")
mean <- as.table(mean)
# Calculate standardized mean difference
smdMotion <- smd(DATA, RMSD)
smdAge <- smd(DATA, "AGE")
smdFIQ <- smd(DATA, "FIQ")


#===Tabling smd for output ===
smdOpt <-
        matrix(c(smdMotion, smdAge, smdFIQ),ncol = 3)
colnames(smdOpt) <-
        c("RMSD",
          "AGE",
          "FIQ")
rownames(smdOpt) <- "SMD"
smdOpt <- as.table(smdOpt)
#return(smdMotion )#= smdMotion)#,smdAge = smdAge, smdOpt = smdOpt))
return(list(smdMotion = smdMotion,smdAge = smdAge,smdFIQ = smdFIQ,
                       print(pvalsOpt),print(smdOpt),print(table(DATA$GROUP)),print(mean)))
#cat("pvalsOpt", "smdOpt","table(DATA$GROUP)","\n")
}