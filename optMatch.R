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

Gmatch <- function (data, formula, nTree, method, distance) {
        #==== setting parameters ====
        response <- all.vars(formula)[1]
        # response <- deparse(substitute(response))
        data[[response]]
        #==== Cleaning data ====
        # It is important to have ASD first and then TD for distance matrix
        #?????????????????????????? how to not hard code variables.
        
        
        data[[response]] <- as.factor (data[[response]])
        
        # Number of subjects in each response for original data
        table(data[[response]])
        
        # screening for outlier in motion variable
        outlierVal <-
                mean(data$RMSD.PRE.censoring) + 3 * sd(data$RMSD.PRE.censoring)
        
        # Create the dataframe for matching without outlier and label rows
        Gdata <-
                data[which(data$RMSD.PRE.censoring <= outlierVal), all.vars(formula)]
        
        rownames(Gdata) <- data[, "subj"]
        
        # Number of subjects in each group for filtered data
        table (Gdata[[response]])
        
        # Outlier Subjects
        outliers <-
                as.character(data[which(data$RMSD.PRE.censoring >= outlierVal), ]$subj)
        cat("subjects with extreme motion are", outliers, "\n")
        
        
        #==== Tree building ====
        # Building the tree using covariates using greedy search for binary response
        #  which uses gini index for split. stop criteria for splitting is thirthy
        # subjects. At each split we use randomly three variables.
        #????????????????Hard coding variables
        # form <-
        #   group ~ Gender + Handedness + RMSD.PRE.censoring + Age + WASI.NVIQ
        mydata <- grow(
                formula,
                data = Gdata ,
                search = "greedy" ,
                method = "class" ,
                split = "gini" ,
                minsplit = 10,
                mtry = 3,
                nsplit = NULL
        )
        # Plot the tree
        plot(mydata)
        
        #==== variable pvalue before matching ====
        
        pvalNumeric <-
                (lapply(Gdata[, which(sapply(Gdata, is.numeric))],
                        function(x)
                                (t.test(x ~ Gdata[[all.vars(formula)[1]]]))$p.value))
        pvalFact <-
                (lapply(Gdata[, which(sapply(Gdata, is.factor))],
                        function(x)
                                (chisq.test(Gdata[[all.vars(formula)[1]]], x))$p.value))
        pval <- append(pvalNumeric, pvalFact)
        pval[[4]] <- NULL
        print(unlist(pval))
        
        # Calculate standardized mean difference
        # smd(nonOutlierLowADOS,all.vars(formula)[[2]])
        #
        # smdNumeric <- (lapply(Gdata[, which(sapply(Gdata, is.numeric))],
        #                        function(x) smd(Gdata, RMSD.PRE.censoring)))
        #
        smdMotion <- smd(Gdata, RMSD.PRE.censoring)
        smdAge <- smd(Gdata, "Age")
        smdNVIQ <- smd(Gdata, "WASI.NVIQ")
        smdGender <- smd(Gdata, "Gender")
        smdHandedness <- smd(Gdata, "Handedness")
        
        # Tabling smd for output
        smdBef <-
                matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
                       ncol = 5)
        colnames(smdBef) <-
                c("RMSD.PRE.censoring",
                  "Age",
                  "WASI.NVIQ",
                  "Handedness",
                  "Gender")
        rownames(smdBef) <- "BEFORE"
        smdBef <- as.table(smdBef)
        smdBef
        #==== Random Forest growing ====
        set.seed(184684)
        ntrees <- nTree
        
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
                        minbucket = 10,
                        bootstrap = FALSE
                )
        
        #==== Reference group for matching to start====
        
        # Which group has less subject. That group will be our reference for matching.
        ASD <- Gdata$group == 1
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
        # 1- 1To3GH: 1-3 filtering subject on gender and handedness.
        #           Based on distance percentile, certain threshold will be picked.
        if (method == "1To3GH") {
                #Tight matching on gender and handedness using certain distance thershold
                
                # ASD and TD distance only and label them.
                x2distForSel <-
                        xDistFor[1:nMin, (nMin + 1):ncol(xDistFor)]
                rownames(x2distForSel) <- rownames(Gdata[1:nMin, ])
                colnames(x2distForSel) <-
                        rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
                
                #filter TD gender that are equal to ASD gender label to 1 otherwise 0.
                library(foreach)
                filtGen <- foreach (i = 1:nMin) %do%
                        ifelse(Gdata[rownames(x2distForSel)[i],
                                     "Gender"] == Gdata[colnames(x2distForSel), "Gender"] , 1, 0)
                
                filtGen <-
                        matrix(unlist(filtGen),
                               ncol = ncol(x2distForSel),
                               byrow = T)
                colnames(filtGen) <-
                        rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
                rownames(filtGen) <- rownames(Gdata[1:nMin, ])
                #filter TD handedness that are equal to ASD handedness label to 1 otherwise 0
                filtHand <- foreach (i = 1:nMin) %do%
                        ifelse(Gdata[rownames(x2distForSel)[i],
                                     "Handedness"] == Gdata[colnames(x2distForSel), "Handedness"] , 1, 0)
                
                filtHand <-
                        matrix(
                                unlist(filtHand),
                                ncol = ncol(x2distForSel),
                                byrow = T
                        )
                colnames(filtHand) <-
                        rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
                rownames(filtHand) <- rownames(Gdata[1:nMin, ])
                
                # Add these two conditions
                filtGenHand <- filtGen + filtHand
                rownames(filtGenHand) <- rownames(Gdata[1:nMin,])
                colnames(filtGenHand) <-
                        rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
                
                # Distance thereshold
                distFiltVal <- 0.7
                
                # Select TD subjects tha t qualify these 3 conditions.
                filtGenHandDist <- matrix(
                        ifelse(
                                x2distForSel[,] <= distFiltVal & filtGenHand[,] == 2,
                                x2distForSel,
                                1
                        ),
                        nrow = nMin ,
                        byrow = F
                )
                rownames(filtGenHandDist) <-
                        rownames(Gdata[1:nMin, ])
                colnames(filtGenHandDist) <-
                        rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
                
                # Create distance matrix csv file.
                outFileName = paste("output//x2DistMat",
                                    distFiltVal,
                                    ".csv",
                                    sep = "")
                write.table(
                        filtGenHandDist,
                        outFileName,
                        sep = ",",
                        col.names = TRUE
                )
                library("optmatch")
                DM <- filtGenHandDist
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
                print(optMatch)
                print(summary(optMatch))
                capture.output(print(optMatch, grouped = T), file = "subjMatch.txt")
                subjMatch <- read.table("subjMatch.txt")
                optData <-
                        Gdata[c(as.character(subjMatch[[1]]),
                                as.character(subjMatch[[2]])),]
                #==== 1-1 optimal matching =====
        } else if (method == "1-1") {
                #==== Optimal matching ====
                library("optmatch")
                # Optimal one-to-one matching
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
                print(optMatch, responseed = T)
                
                # write output to a file to be used again as an obj in R
                capture.output(print(optMatch, grouped = T), file = "subjMatch.txt")
                subjMatch <- read.table("subjMatch.txt",header = T)
                splitSubj <- do.call("rbind",strsplit(c(as.character(subjMatch$Group),as.character(subjMatch$Members)), ","))
                optData <- Gdata[splitSubj, ]
               
                #==== Tabling pvalues for output ====
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
                
                #==== Tabling smd for output ====
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
                
                # Check SMD for all variables. if it is not below 10% keep improving
                while (abs(smdMotion) > 0.1 | abs(smdAge) > 0.1| abs(smdNVIQ) > 0.1 | abs(smdGender) > 0.1  |abs(smdHandedness) > 0.1 ){
                        #= Can we improve SMD? First =
                        # who has the largest distance? 
                        distance <- NULL
                        for (i in 1:(table(optData$group)[2])){
                                distance <- c(distance,(DM[as.character(splitSubj[i,1]),as.character(splitSubj[(i+(table(optData$group)[2])),1])]) )    
                        }
                        splitSubj <- cbind(matrix(splitSubj,nrow = table(optData$group)[2]),distance)
                        # Find the name of ASD subject has the largest distance
                        excSubj <- (splitSubj[apply(splitSubj,2,which.max)$distance, ])[1]
                        
                        #Get the list of ASD subject that partcipate in the matching
                        allMatchSubj =rownames(DM[splitSubj[1:(nrow(splitSubj)), 1], ])
                        
                        #remaining ASD subjects
                        remainedASD <- allMatchSubj[allMatchSubj!=excSubj ]
                        
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
                        
                        # write the output to a text file
                        capture.output(print(optMatch, grouped = T), file = "subjMatch.txt")
                        
                        subjMatch <- read.table("subjMatch.txt",header = T)
                        splitSubj <- do.call("rbind",strsplit(c(as.character(subjMatch$Group),as.character(subjMatch$Members)), ","))
                        optData <- Gdata[splitSubj, ]
                        #====  # Tabling pvalues and SMD after matching ====
                        
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
                print(rownames(optData))
                                       
        #====
        #         matchASD <-c( "022A","096A","111A","115A","123A","131A","136A","138A","141A",
        #               "144A","147A","038A","148A","150A","153A","155A","157A","159A_V2",
        #               "160A","162A","167A","174A","045A","175A","181A","191A","192A",
        #               "194A","198A_V2","205A","209A","224A","225A","056A","226A","404A",
        #               "405A","406A_V2","409A","410A","413A","420A","071A","082A","091A",
        #               "092A","095A")
        #
        #         matchTD <- c("415C","064C","104C","090C","122C","172C","088C","100C","154C","151C",
        #              "197C","121C","407C","106C","146C","119C","109C","097C","166C","185C",
        #              "8C_V2","183C","078C","067C","179C","402C","049C","086C","133C","065C",
        #              "152C","128C","182C","098C","101C","108C","178C","195C","186C","165C",
        #              "048C","168C","110C","125C","211C","176C","173C"
        # )
        
        #        optData<- Gdata[c(matchASD,matchTD), ]
        
        }
}
