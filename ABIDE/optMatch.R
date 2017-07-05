library(optmatch)
#==============================================================================
# Copyright statement comment
# Author: Afrooz Jahedi
# Goal: Create groupwise (1-1) matching.
# Inputs: data, variables to match on, number of trees, method of matching
# Outputs:pvalues and standardized mean difference for each variable.number
#  of subjects in each group
# Description: build a distance matrix based on random forest with 0-1 defnition or 
# chi-square defnition.
#==============================================================================
ptm <- proc.time()
sink("ala.ABIDE.1000.txt")

Gmatch <- function (data, formula, nTree, method, distance) {
        #==== setting parameters ====
        response <- all.vars(formula)[1]
        # response <- deparse(substitute(response))
        data[[response]]
        #==== Cleaning data ====
        # It is important to have ASD first and then TD for distance matrix
        #?????????????????????????? how to not hard code variables.
#sink()        
        
        data[[response]] <- as.factor (data[[response]])
        
        # Number of subjects in each response for original data
        table(data[[response]])
       
  
        complete.cases(data)
        # screening for outlier in motion variable
        outlierVal <-
                mean(data$RMSD) + 3 * sd(data$RMSD)
        
        # Create the dataframe for matching without outlier and label rows
        Gdata <-
                data[which(data$RMSD <= outlierVal), all.vars(formula)]
        
        rownames(Gdata) <- data[, "SUB_ID"]
        
        # Number of subjects in each group for filtered data
        table (Gdata[[response]])
        
        # Outlier Subjects
        outliers <-
                as.character(data[which(data$RMSD >= outlierVal), ]$subj)
        cat("subjects with extreme motion are", outliers, "\n")
        
        
        #==== Tree building ====
        # Building the tree using covariates using greedy search for binary response
        #  which uses gini index for split. stop criteria for splitting is thirthy
        # subjects. At each split we use randomly three variables.
        #????????????????Hard coding variables
        # form <-
        #   group ~ Gender + Handedness + RMSD + Age + WASI.NVIQ
        # mydata <- grow(
        #         formula,
        #         data = Gdata ,
        #         search = "greedy" ,
        #         method = "class" ,
        #         split = "gini" ,
        #         minsplit = 20,
        #         mtry = 3,
        #         nsplit = NULL
        # )
        # # Plot the tree
        # plot(mydata)
        
        #==== variable pvalue before matching ====
        summaryGmatch(Gdata)
       
        # #==== Random Forest growing ====
        # set.seed(184684)
         ntrees <- nTree
        # 
        # Create random forest using the same matching covariates using all data
        rfRF100 <-
                rfConst(
                        ntrees = ntrees,
                        formula = form,
                        training =  Gdata,
                        search = "greedy",
                        method = "class",
                        split = "gini",
                        mtry = 3,
                        nsplit = NULL,
                        minsplit = 20,
                        maxdepth = 10,
                        minbucket = 5,
                        bootstrap = FALSE
                )
        
        #==== Reference group for matching to start====
        
        # Which group has less subject. That group will be our reference for matching.
        ASD <- Gdata$GROUP == 1
        if (NROW (Gdata[ASD,]) < NROW (Gdata[!ASD, ])) {
                (nMin = NROW (Gdata[ASD, ]))
        } else {
                nMin <- NROW (Gdata[!ASD, ])
        }
        #==== distance matrix =====
        if (distance=="chi"){
                # Extracting terminal node for subjects across all trees.
                nodeResponse <- sapply(rfRF100, function(x) {
                predict(x$tree, newdata = Gdata , type = "node")
        })
        # Preparing distance matrix
        sumXDistMat = matrix(0,
                             nrow = NROW(nodeResponse),
                             ncol = NROW(nodeResponse))
        colnames(sumXDistMat) <- rownames(nodeResponse)
        rownames(sumXDistMat) <- rownames(nodeResponse)
        
        ntrees <- nTree
        
        # Define a function x2Dist
        #x2Dist(rfRF100[[treeNum]]$tree,nodeResponse,treeNum)
        
        # Add up tree p-values
        for (treeNum in 1:ntrees) {
                sumXDistMat <-
                        sumXDistMat + x2Dist(rfRF100[[treeNum]]$tree, nodeResponse, treeNum)
        }
        
        #average the distance matrix and label rows and columns.
        xDistFor <- sumXDistMat / ntrees
        rownames(xDistFor) <- rownames(Gdata)
        colnames(xDistFor) <- rownames(Gdata)
   
        # ASD and TD distance only and label them.
        x2distForSel <- xDistFor[1:nMin, (nMin + 1):ncol(xDistFor)]
        rownames(x2distForSel) <- rownames(Gdata[1:nMin, ])
        colnames(x2distForSel) <-
                rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
        #replace 0 with 0.1 and 1 with 0.9.
        
        DM <- x2distForSel
        
        }else if(distance=="0-1"){
                # predict function can provide node info for each obs. Using sapply use all
                # info from all trees.Node response contains a mtrix of number subjectsx #trees
                # it shows each subject ended up in which terminal node. #we get thiese terminal
                # node information for all trees.
                nodeResponse <-
                        sapply(rfRF100, function(x) {
                                predict(x$tree, newdata = Gdata , type = "node")
                        })
                
                # For creating distance matrix we should set a default matrix first
                distMat <-
                        matrix(ntrees,
                               nrow = NROW(nodeResponse),
                               ncol = NROW(nodeResponse))
                # Distance defines as 1 if subjects are in same terminal node,0 otherwise.
                for (tree in 1:ntrees) {
                        for (row in 1:NROW(nodeResponse)) {
                                distMat[, row] <-
                                        distMat[, row] - (nodeResponse[row, tree] == nodeResponse[, tree])
                        }
                }
                
                # Average of distance across all trees.
                distForest <- distMat / ntrees
                # average distance in the forest
                
                rownames(distForest) <- rownames(Gdata)
                colnames(distForest) <- rownames(Gdata)
                distForSel <- distForest[1:nMin, (nMin + 1):ncol(distForest)]
                rownames(distForSel) <- rownames(Gdata[1:nMin, ])
                colnames(distForSel) <-
                        rownames(Gdata[(nMin + 1):ncol(distForest), ])
                DM <- distForSel
        }
        #==== Methods ====
        if (method == "opt-coars-exact") {
                #=== Distance prepration ====
                # ASD and TD distance only and label them.
                x2distForSel <-
                        xDistFor[1:nMin, (nMin + 1):ncol(xDistFor)]
                rownames(x2distForSel) <- rownames(Gdata[1:nMin, ])
                colnames(x2distForSel) <-
                        rownames(Gdata[(nMin + 1):ncol(xDistFor),])
                
                # All continious variables of data,dimension and name
                numData <-
                        Gdata[, split(names(Gdata), sapply(Gdata, function(x)
                                paste(class(x), collapse = " ")))$numeric]
                nv <- dim(numData)[2]
                vnames <- colnames(numData)
                
                #=== Discritize the selected var ===
                for (i in (1:nv)) {
                        tmp <- reduceVar(numData[[i]], "quantile")
                        numData[[i]] <- tmp$x
                        numData[[i]] <- as.factor(numData[[i]])
                }
                #=== Exact matching ====
                library(foreach)
                # Exact matching of ASD subjects based on categorical variables. If exact match, it's distance, otherwise 1.
                filtDataCat <- foreach (i = 1:nMin) %do%
                        ifelse(Gdata[rownames(x2distForSel)[i], "Gender"]
                               == Gdata[colnames(x2distForSel), "Gender"] &
                                       Gdata[rownames(x2distForSel)[i], "Handedness"]
                               == Gdata[colnames(x2distForSel), "Handedness"]  , x2distForSel, 1)
                
                filtDataCat <-
                        matrix(unlist(filtDataCat),
                               ncol = ncol(x2distForSel),
                               byrow = T)
                colnames(filtDataCat) <-
                        rownames(Gdata[(nMin + 1):ncol(xDistFor),])
                rownames(filtDataCat) <- rownames(Gdata[1:nMin,])
                
                # Exact matching of ASD subjects based on continious variables. If exact match, it's distance,otherwise 1.
                 thr <- 100
                 filtDataCont <- foreach (i = 1:nMin) %do%
                        ifelse(
                                numData[rownames(x2distForSel)[i], "RMSD"]
                                == numData[colnames(x2distForSel), "RMSD"] &
                                        numData[rownames(x2distForSel)[i], "Age"]
                                == numData[colnames(x2distForSel), "Age"] &
                                        numData[rownames(x2distForSel)[i], "WASI.NVIQ"]
                                == numData[colnames(x2distForSel), "WASI.NVIQ"],
                                x2distForSel,
                                thr
                        )
                
                # Build a distance matrix based on exact matched.label cols and rows.
                filtDataCont <- matrix(unlist(filtDataCont),
                                       ncol = ncol(x2distForSel),
                                       byrow = T)
                colnames(filtDataCont) <-
                        rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
                rownames(filtDataCont) <- rownames(Gdata[1:nMin, ])
                write.csv(filtDataCont,file = "filtDataCont.csv")
                #=== one-to-one optimal matching ====
                library("optmatch")
                DM <- filtDataCont
                optMatch <-
                        fullmatch(
                                DM,
                                min.controls = 1,
                                max.controls = 1,
                                omit.fraction = NULL,
                                tol = 0.001,
                                subclass.indices = NULL
                        )
                matchSubj <- print(optMatch, grouped = T)
                #print(optMatch)
                print(summary(optMatch))
                
                #=== write output into a object and read it from a file ====
                capture.output(print(optMatch, grouped = T), file = "subjMatch.txt")
                subjMatch <- read.table("subjMatch.txt")
                splitSubj <-
                        do.call("rbind", strsplit(c(
                                as.character(subjMatch$Group),
                                as.character(subjMatch$Members)
                        ), ","))
                optData <-
                        Gdata[c(as.character(gsub(',', '', (subjMatch[[1]]))),
                                as.character(subjMatch[[2]])), ]
                # Set the matched data
                matchedData <- optData
                
                #=== Tabling pvalues and SMD after matching ====
                
                pvalMotion <-
                        (t.test(RMSD ~ group, data = matchedData))$p.value
                pvalAge <-
                        (t.test(Age ~ group, data = matchedData))$p.value
                pvalNVIQ <-
                        (t.test(WASI.NVIQ ~ group, data = matchedData))$p.value
                pvalGender <-
                        chisq.test(matchedData$group, matchedData$Gender)$p.value
                pvalHandedness <-
                        chisq.test(matchedData$group, matchedData$Handedness)$p.value
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
                        c("RMSD",
                          "Age",
                          "WASI.NVIQ",
                          "Gender",
                          "Handedness")
                rownames(pvalsOpt) <- "Pvals"
                pvalsOpt <- as.table(pvalsOpt)
                print(pvalsOpt)
                
                # Calculate standardized mean difference
                smdMotion <- smd(matchedData, RMSD)
                smdAge <- smd(matchedData, "Age")
                smdNVIQ <- smd(matchedData, "WASI.NVIQ")
                smdGender <- smd(matchedData, "Gender")
                smdHandedness <- smd(matchedData, "Handedness")
                
                # Tabling smd for output
                smdOpt <-
                        matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
                               ncol = 5)
                colnames(smdOpt) <-
                        c("RMSD",
                          "Age",
                          "WASI.NVIQ",
                          "Gender",
                          "Handedness")
                rownames(smdOpt) <- "SMD"
                smdOpt <- as.table(smdOpt)
                print(smdOpt)
                print(table(matchedData$group))
                print(rownames(matchedData))
        
                #=== itteretive optimal matching ====
                # Check SMD for all variables. if it is not below 10% keep improving 
                while (abs(smdMotion) > 0.1 |
                       abs(smdAge) > 0.1 |
                       abs(smdNVIQ) > 0.1 |
                       abs(smdGender) > 0.1  |
                       abs(smdHandedness) > 0.1)
                {
                        #=== Obtain matched ASD subject distance ====
                        distance <- NULL
                        for (i in 1:(table(optData$group)[2])){
                                distance <- c(distance,(DM[splitSubj[i,1],splitSubj[(i+(table(optData$group)[2])),1]]) )    
                        }
                        splitSubj <- data.frame(splitSubj,distance)
                        
                        # Find the name of ASD subject has the largest distance
                        # subset(splitSubj,distance!=100,select = which.max(distance))
                        # (subset(splitSubj,distance!=100,select =(distance)))
                        maxSubjs <- (splitSubj[ splitSubj[,"distance"]!=100 ,]) 
                        excSubj <- maxSubjs[which.max(maxSubjs[,2]),1]
                        print(excSubj)
                        #Get the list of ASD subject that partcipate in the matching
                        allMatchSubj =splitSubj[1:(table(optData$group)[2]),1]
                        # 
                        # for (i in 1:(table(matchedData$group)[2])) {
                        #         ASDdistance <-
                        #                 c(ASDdistance, (DM[rownames((matchedData)[i,]),rownames(matchedData[(i +table(matchedData$group)[2]), ])]))
                        # }
                        # 
                        # # Distance of matched ASD subjects
                        # ASDdistance <- data.frame(subjId=(rownames(matchedData))[1:(table(matchedData$group))[2]], dist=ASDdistance)
                        # 
                        # #=== Find the largest distance ASD name that is not exclude due to exact matching.====
                        # excSubj <- (ASDdistance[which.max(ASDdistance$dist!=thr),])$subjId
                        # print(excSubj)
                        # print(ASDdistance[excSubj,])
                        # #Get the list of ASD subject that partcipate in the matching
                        # allMatchSubj = rownames(matchedData[1:(table(matchedData$group)[2]), ])
                        # 
                        #remaining ASD subjects
                        remainedASD <- as.character(allMatchSubj[allMatchSubj!= excSubj])
                        print(remainedASD)
                        
                        #trimed distance matrix
                        trimedDM <- DM[remainedASD,]
                        

                        #====Redo the optimal matching ======
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
                        #print(optMatch, responseed = T)
                        
                        # #==== write the itteretive matched output to a text file ====
                        capture.output(print(optMatch, grouped = T), file = "iterSubjMatch.txt")
                        
                        subjMatch <- read.table("iterSubjMatch.txt",header = T)
                        splitSubj <- do.call("rbind",strsplit(c(as.character(subjMatch$Group),as.character(subjMatch$Members)), ","))
                        optData <- Gdata[splitSubj, ]
                        # capture.output(print(optMatch, grouped = T), file = "iterSubjMatch.txt")
                        # subjMatch <- read.table("iterSubjMatch.txt")
                        # optData <-
                        #         Gdata[c(as.character(gsub(',', '', (subjMatch[[1]]))),
                        #                 as.character(subjMatch[[2]])), ]
                         matchedData <- optData
                        #====  # Tabling pvalues and SMD after matching ====
                        summaryGmatch(optData)
                        
                                      }

        }else if (method == "opt-coars-exact-rev"){
                #==== exact revised matching ====
                library(fields)
                motionDist <- rdist(Gdata$RMSD)
                motionDist[lower.tri(motionDist)] <- NA
                diag(motionDist) <- NA
                colnames(motionDist) <- rownames(Gdata)
                rownames(motionDist) <- rownames(Gdata)
                motionDist <- motionDist[1:nMin,(nMin+1):nrow(Gdata)]
                #motionDist <- ifelse(motionDist >.05, 100,motionDist)
                filtMotionDist <- motionDist[,]<0.151
                
                # distance matrix for age
                ageDist <- rdist(Gdata$Age)
                ageDist[lower.tri(ageDist)] <- NA
                diag(ageDist) <- NA
                colnames(ageDist) <- rownames(Gdata)
                rownames(ageDist) <- rownames(Gdata)
                ageDist <- ageDist[1:nMin,(nMin+1):nrow(Gdata)]
                #ageDist <- ifelse(ageDist>2,100,ageDist)
                filtAgeDist <- ageDist[,]<9.9
                
                # NVIQDist <- rdist(Gdata$WASI.NVIQ)
                # NVIQDist[lower.tri(NVIQDist)] <- NA
                # diag(NVIQDist) <- NA
                # colnames(NVIQDist) <- rownames(Gdata)
                # rownames(NVIQDist) <- rownames(Gdata)
                # NVIQDist <- NVIQDist[1:nMin,(nMin+1):nrow(Gdata)]
                # #ageDist <- ifelse(ageDist>2,100,ageDist)
                # filtNVIQDist <- NVIQDist[,]<10
                
                #combine exact matching
                combineExactDist <- filtAgeDist&filtMotionDist#&filtNVIQDist
                combineExactDist <- ifelse(combineExactDist,DM, 10) 
                a <- iterOpt(combineExactDist,Gdata) 
        } else if (method == "opt-one-to-one"){

                #=== Optimal one-to-one matching ====
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
                optData <- Gdata[splitSubj, ]
               
                #=== Tabling pvalues and SMD for output ====
                S <- summaryGmatch(optData)
               
                #=== itteration ===== 
                # Check SMD for all variables. if it is not below 10% keep improving
                while (abs(S$smdMotion) > 0.1 | abs(S$smdAge) > 0.1| abs(S$smdFIQ) > 0.1 )
                        {
                        #= Can we improve SMD? First =
                        # who has the largest distance? 
                        distance <- NULL
                        for (i in 1:(table(optData$GROUP)[2])){
                                distance <- c(distance,(DM[as.character(splitSubj[i,1]),as.character(splitSubj[(i+(table(optData$GROUP)[2])),1])]) )    
                        }
                        splitSubj <- data.frame(splitSubj,distance)
                        # Find the name of ASD subject has the largest distance
                        #excSubj <- (splitSubj[apply(splitSubj,2,which.max)$distance, ])[1]
                        
                        maxSubjs <- (splitSubj[ splitSubj[,"distance"]!=100 ,]) 
                        excSubj <- maxSubjs[which.max(maxSubjs[,2]),1]
                        print(excSubj)
                        #Get the list of ASD subject that partcipate in the matching
                        #allMatchSubj =rownames(DM[splitSubj[1:(nrow(splitSubj)), 1], ])
                        allMatchSubj =splitSubj[1:(table(optData$GROUP)[2]),1]
                        
                        
                        #remaining ASD subjects
                        remainedASD <- as.character(allMatchSubj[allMatchSubj!= excSubj])
                        #print(remainedASD)
                        
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
                        #print(optMatch, responseed = T)
                        
                        #==== write the itterated matched output to a text file ====
                        capture.output(print(optMatch, grouped = T), file = "subjMatch.txt")
                        
                        subjMatch <- read.table("subjMatch.txt",header = T)
                        splitSubj <- do.call("rbind",strsplit(c(as.character(subjMatch$Group),as.character(subjMatch$Members)), ","))
                        optData <- Gdata[splitSubj, ]
                        #==== Tabling pvalues and SMD after matching ====
                        S <- summaryGmatch(optData)
                        
                }
        }                    
      
        
        proc.time() - ptm
        sink()
          }
