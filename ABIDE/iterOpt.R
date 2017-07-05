iterOpt <- function(DM,data){
library("optmatch")
optMatch <-
        fullmatch(
                DM,
                min.controls = 1,
                max.controls = 1,
                omit.fraction = NULL,
                tol = 0.001,
                subclass.indices = NULL
        )
summary(optMatch)
#print(optMatch, responseed = T)

#=== write output to a file to be used again as an obj in R ====
capture.output(print(optMatch, grouped = T), file = "subjMatch.txt")
subjMatch <- read.table("subjMatch.txt", header = T)
splitSubj <-
        do.call("rbind", strsplit(c(
                as.character(subjMatch$Group),
                as.character(subjMatch$Members)
        ), ","))
optData <- data[splitSubj, ]

#=== Tabling pvalues and SMD for output ====
pvalMotion <-
        (t.test(RMSD.PRE.censoring ~ group, data = optData))$p.value
pvalAge <- (t.test(Age ~ group, data = optData))$p.value
pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = optData))$p.value
pvalGender <- chisq.test(optData$group, optData$Gender)$p.value
pvalHandedness <-
        chisq.test(optData$group, optData$Handedness)$p.value

pvalsOpt <-
        matrix(c(
                pvalMotion,
                pvalAge,
                pvalNVIQ,
                pvalGender,
                pvalHandedness
        ),
        ncol = 5)
colnames(pvalsOpt) <-
        c("RMSD.PRE.censoring",
          "Age",
          "WASI.NVIQ",
          "Gender",
          "Handedness"
        )
rownames(pvalsOpt) <- "Pvals"
pvalsOpt <- as.table(pvalsOpt)
pvalsOpt

# Calculate standardized mean difference
smdMotion <- smd(optData, RMSD.PRE.censoring)
smdAge <- smd(optData, "Age")
smdNVIQ <- smd(optData, "WASI.NVIQ")
smdGender <- smd(optData, "Gender")
smdHandedness <- smd(optData, "Handedness")

#===Tabling smd for output ===
smdOpt <-
        matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
               ncol = 5)
colnames(smdOpt) <-
        c("RMSD.PRE.censoring",
          "Age",
          "WASI.NVIQ",
          "Gender",
          "Handedness"
        )
rownames(smdOpt) <- "SMD"
smdOpt <- as.table(smdOpt)
smdOpt
table(optData$group)

#=== itteration ===== 
# Check SMD for all variables. if it is not below 10% keep improving
while (abs(smdMotion) > 0.1 | abs(smdAge) > 0.1| abs(smdNVIQ) > 0.1 | abs(smdGender) > 0.1  |abs(smdHandedness) > 0.1 )
{
        #= Can we improve SMD? First =
        # who has the largest distance? 
        distance <- NULL
        for (i in 1:(table(optData$group)[2])){
                distance <- c(distance,(DM[as.character(splitSubj[i,1]),as.character(splitSubj[(i+(table(optData$group)[2])),1])]) )    
        }
        splitSubj <- data.frame(splitSubj,distance)
        # Find the name of ASD subject has the largest distance
        #excSubj <- (splitSubj[apply(splitSubj,2,which.max)$distance, ])[1]
        if (all(splitSubj[,"distance"]==10)){
                maxSubjs <- (splitSubj[ splitSubj[,"distance"]==10 ,]) 
                excSubj <- maxSubjs[which.max(maxSubjs[,2]),1]
        }else {
                excSubj <- splitSubj[which.max(splitSubj[,2]),1]
        }
        print(excSubj)
        #Get the list of ASD subject that partcipate in the matching
        #allMatchSubj =rownames(DM[splitSubj[1:(nrow(splitSubj)), 1], ])
        allMatchSubj =splitSubj[1:(table(optData$group)[2]),1]
        
        
        #remaining ASD subjects
        remainedASD <- as.character(allMatchSubj[allMatchSubj!= excSubj])
        print(remainedASD)
        
        #trimed distance matrix
        trimedDM <- DM[remainedASD,]
        
        #Redo the optimal matching
        optMatch <-
                fullmatch(
                        trimedDM,
                        min.controls = 1,
                        max.controls = 1,
                        omit.fraction = NULL,
                        tol = 0.001,
                        subclass.indices = NULL
                )
        summary(optMatch)
        print(optMatch, responseed = T)
        
        #==== write the itterated matched output to a text file ====
        capture.output(print(optMatch, grouped = T), file = "subjMatch.txt")
        
        subjMatch <- read.table("subjMatch.txt",header = T)
        splitSubj <- do.call("rbind",strsplit(c(as.character(subjMatch$Group),as.character(subjMatch$Members)), ","))
        optData <- data[splitSubj, ]
        #==== Tabling pvalues and SMD after matching ====
        
        pvalMotion <-
                (t.test(RMSD.PRE.censoring ~ group, data = optData))$p.value
        pvalAge <- (t.test(Age ~ group, data = optData))$p.value
        pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = optData))$p.value
        pvalGender <- chisq.test(optData$group, optData$Gender)$p.value
        pvalHandedness <-
                chisq.test(optData$group, optData$Handedness)$p.value
        
        
        pvalsOpt <-
                matrix(c(
                        pvalMotion,
                        pvalAge,
                        pvalNVIQ,
                        pvalGender,
                        pvalHandedness
                ),
                ncol = 5)
        colnames(pvalsOpt) <-
                c("RMSD.PRE.censoring",
                  "Age",
                  "WASI.NVIQ",
                  "Gender",
                  "Handedness"
                )
        rownames(pvalsOpt) <- "Pvals"
        pvalsOpt <- as.table(pvalsOpt)
        print(pvalsOpt)
        # Calculate standardized mean difference
        smdMotion <- smd(optData, RMSD.PRE.censoring)
        smdAge <- smd(optData, "Age")
        smdNVIQ <- smd(optData, "WASI.NVIQ")
        smdGender <- smd(optData, "Gender")
        smdHandedness <- smd(optData, "Handedness")
        
        # Tabling smd for output
        smdOpt <-
                matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
                       ncol = 5)
        colnames(smdOpt) <-
                c("RMSD.PRE.censoring",
                  "Age",
                  "WASI.NVIQ",
                  "Gender",
                  "Handedness"
                )
        rownames(smdOpt) <- "SMD"
        smdOpt <- as.table(smdOpt)
        print(smdOpt)
        print(table(optData$group))
}
}