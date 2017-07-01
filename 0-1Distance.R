#==============================================================================
# Copyright statement comment
# Author: Afrooz Jahedi
# Goal: Create groupwise (1-3) matching.
# Inputs: data, variables to match on, number of trees, method of matching
# Outputs:pvalues and standardized mean difference for each variable.number 
#  of subjects in each group
# Description: build a distance matrix based on random forest where subject distance  
# in the same terminal node is zero and subjects in different terminal nodes have
#  distance equal to 1.
#==============================================================================


Gmatch <- function (data, formula, nTree, method) {
        response <- all.vars(formula)[1]
        # response <- deparse(substitute(response))
        data[[response]]
        #==== Cleaning data ####
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
                as.character(data[which(data$RMSD.PRE.censoring >= outlierVal),]$subj)
        cat("subjects with extreme motion are", outliers, "\n")
        
        
        ##### Tree building ####
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
                minsplit = 20,
                mtry = 3,
                nsplit = NULL
        )
        # Plot the tree
        plot(mydata)
        
        ##### variable pvalue before matching ####
        summaryGmatch(Gdata)
        # pvalNumeric <- (lapply(Gdata[, which(sapply(Gdata, is.numeric))],
        #                        function(x)
        #                                (t.test(x ~ Gdata[[all.vars(formula)[1]]]))$p.value))
        # pvalFact <- (lapply(Gdata[, which(sapply(Gdata, is.factor))],
        #                     function(x)
        #                             (chisq.test(Gdata[[all.vars(formula)[1]]], x))$p.value))
        # pval <- append(pvalNumeric, pvalFact)
        # pval[[4]] <- NULL
        # print(unlist(pval))
        # 
        # # Calculate standardized mean difference
        # # debug(smd)
        # # smd(nonOutlierLowADOS,all.vars(formula)[[2]])
        # #
        # # smdNumeric <- (lapply(Gdata[, which(sapply(Gdata, is.numeric))],
        # #                        function(x) smd(Gdata, RMSD.PRE.censoring)))
        # #
        # smdMotion <- smd(Gdata, RMSD.PRE.censoring)
        # smdAge <- smd(Gdata, "Age")
        # smdNVIQ <- smd(Gdata, "WASI.NVIQ")
        # smdGender <- smd(Gdata, "Gender")
        # smdHandedness <- smd(Gdata, "Handedness")
        # 
        # # Tabling smd for output
        # smdBef <-
        #         matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
        #                ncol = 5)
        # colnames(smdBef) <-
        #         c("RMSD.PRE.censoring",
        #           "Age",
        #           "WASI.NVIQ",
        #           "Handedness",
        #           "Gender")
        # rownames(smdBef) <- "BEFORE"
        # smdBef <- as.table(smdBef)
        # smdBef
        ##### Random Forest growing ####
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
                        minbucket = 5,
                        bootstrap = FALSE
                )
        
        ##### Reference group for matching to start####
        
        # Which group has less subject. That group will be our reference for matching.
        ASD <- Gdata$group == 1
        if (NROW (Gdata[ASD, ]) < NROW (Gdata[!ASD,])) {
                (nMin = NROW (Gdata[ASD,]))
        } else {
                nMin <- NROW (Gdata[!ASD,])
        }
        ##### Disstance matrix ####
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
        #write.table(xDistFor,  sep=",", col.names=TRUE)
        
        ##### Group-wise ,matching methods ####
        # 1- 1To3GH: exact 1-3 matching e.g., filtering subject X2-distance for Gender and handedness.
        #           Based on distance percentile, certain threshold will be picked.
        
        if(method =="1To3Dist") {
                selCandid <- 3
        candidVal <-
                apply(distForest[1:nMin , (nMin + 1):ncol(distForest)], 1, function(x)
                        return((sort(x))[1:selCandid]))
        
        candid <-
                apply(distForest[1:nMin, (nMin + 1):ncol(distForest)], 1, function(x)
                        return((order(x))[1:selCandid] + nMin))
        
        distCandid <- (table(candid))
        names(distCandid) <-
                rownames(Gdata[as.numeric(names(table(candid))), ])
        hist(distCandid, main = "Number of times a subject is selected", col="red")
        
        # Matched data using definde ditance
        dataNew <-
                rbind (Gdata[1:nMin, ], (Gdata[rownames(distCandid), ]))
        #### distance Matrix p-values after matching for  ####
        summaryGmatch(dataNew)
        # pvalHandedness <-
        #         chisq.test(dataNew$group, dataNew$Handedness)$p.value
        # pvalGender <- chisq.test(dataNew$group, dataNew$Gender)$p.value
        # pvalMotion <-
        #         (t.test(RMSD.PRE.censoring ~ group, data = dataNew))$p.value
        # pvalAge <- (t.test(Age ~ group, data = dataNew))$p.value
        # pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = dataNew))$p.value
        # 
        # pvals1_3 <-
        #         matrix(c(pvalMotion, pvalAge, pvalNVIQ, pvalGender, pvalHandedness),
        #                ncol = 5)
        # colnames(pvals1_3) <-
        #         c("RMSD.PRE.censoring",
        #           "Age",
        #           "WASI.NVIQ",
        #           "Handedness",
        #           "Gender")
        # rownames(pvals1_3) <- "1-3"
        # pvals1_3 <- as.table(pvals1_3)
        # print(pvals1_3)
        # 
        # 
        # smdMotion <- smd(dataNew, RMSD.PRE.censoring)
        # smdAge <- smd(dataNew, Age)
        # smdNVIQ <- smd(dataNew, WASI.NVIQ)
        # smdGender <- smd(dataNew, Gender)
        # smdhandedness <- smd(dataNew, Handedness)
        # 
        # smd1_3 <-
        #         matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
        #                ncol = 5)
        # colnames(smd1_3) <-
        #         c("RMSD.PRE.censoring",
        #           "Age",
        #           "WASI.NVIQ",
        #           "Handedness",
        #           "Gender")
        # rownames(smd1_3) <- "1_3"
        # smd1_3 <- as.table(smd1_3)
        # print(smd1_3)
        print(table(dataNew$group))
        }else if (method == "1To3GH") {
        #### Exact matching?####
        distForSel <- distForest[1:nMin, (nMin + 1):ncol(distForest)]
        rownames(distForSel) <- rownames(Gdata[1:nMin,])
        colnames(distForSel) <-
                rownames(Gdata[(nMin + 1):ncol(distForest),])
        outFileName = paste("output//distance matrix", ".csv", sep = "")
        write.table(distForSel,
                    outFileName,
                    sep = ",",
                    col.names = TRUE)
        
        #filter TD gender that are equal to ASD gender to 1 otherwise 0.
        library(foreach)
        filtGen <- foreach (i = 1:nMin) %do%
                ifelse(Gdata[rownames(distForSel)[i],
                             "Gender"] == Gdata[colnames(distForSel), "Gender"] , 1, 0)
        filtGen <-
                matrix(unlist(filtGen),
                       ncol = ncol(distForSel),
                       byrow = T)
        colnames(filtGen) <- rownames(Gdata[(nMin + 1):ncol(distForest),])
        rownames(filtGen) <- rownames(Gdata[1:nMin,])
        
        filtHand <- foreach (i = 1:nMin) %do%
                ifelse(Gdata[rownames(distForSel)[i],
                             "Handedness"] == Gdata[colnames(distForSel), "Handedness"] , 1, 0)
        filtHand <-
                matrix(unlist(filtHand),
                       ncol = ncol(distForSel),
                       byrow = T)
        colnames(filtHand) <- rownames(Gdata[(nMin + 1):ncol(distForest),])
        rownames(filtHand) <- rownames(Gdata[1:nMin,])
        
        filtGenHand <- filtGen + filtHand
        rownames(filtGenHand) <- rownames(Gdata[1:nMin,])
        colnames(filtGenHand) <-
                rownames(Gdata[(nMin + 1):ncol(distForest), ])
        
        distFiltVal <- 0.4
        filtGenHandDist <-
                matrix(
                        ifelse(distForSel[, ] <= distFiltVal & filtGenHand[, ] == 2,
                               distForSel, 100),
                        nrow = nMin ,
                        byrow = F
                )
        rownames(filtGenHandDist) <- rownames(Gdata[1:nMin, ])
        colnames(filtGenHandDist) <-
                rownames(Gdata[(nMin + 1):ncol(distForest), ])
        
        #outFileName = paste("output//Mat", distFiltVal, ".csv", sep="")
        #write.table(filtGenHandDist, outFileName, sep=",", col.names=TRUE)
        
        #ASD trouble makers
        #table((rowSums(filtGenHandDist) >= 5200))
        # Useless TD's
        #table(colSums(filtGenHandDist) >= 4700)
        
        remFilt = filtGenHandDist[!apply(filtGenHandDist, 1, function(x) {
                all(x == 100)
        }), !apply(filtGenHandDist, 2, function(x) {
                all(x == 100)
        })]
        remRowName <- rownames(remFilt)
        remColName <- colnames(remFilt)
        remained <- distForSel[remRowName, remColName]
        outFileName = paste("output//distance matrix 0.5GH", ".csv", sep = "")
        write.table(remained, outFileName, sep = ",", col.names = TRUE)
        
        selCandid <- 3
        candidValFilt <-
                apply(remained, 1, function(x)
                        return((sort(x))[1:selCandid]))
        # who are the TD selected candidates. not TD labels just their index in distForest mat.
        candidFilt <-
                apply(remained, 1, function(x)
                        return(((order(
                                x
                        ))[1:selCandid])))
        
        # Vector of TD candidates
        distFilt <- (table(candidFilt))
        namesTD <-
                colnames(remained[, as.numeric(names(table(candidFilt)))])
        hist(distFilt, main = "Number of times a subject is selected")
        
        ###  Matched data using defined ditance deleting trouble ASD
        
        # Matched data using definde ditance
        dataNewFilt <- rbind (Gdata[remRowName, ], (Gdata[namesTD, ]))
        
        ### p-values after matching 
        
        pvalHandedness <-
                chisq.test(dataNewFilt$group, dataNewFilt$Handedness)$p.value
        pvalGender <-
                chisq.test(dataNewFilt$group, dataNewFilt$Gender)$p.value
        pvalMotion <-
                (t.test(RMSD.PRE.censoring ~ group, data = dataNewFilt))$p.value
        pvalAge <- (t.test(Age ~ group, data = dataNewFilt))$p.value
        pvalNVIQ <-
                (t.test(WASI.NVIQ ~ group, data = dataNewFilt))$p.value
        
        # Tabling pvalues for output
        pvals1_3GH <-
                matrix(c(pvalMotion, pvalAge, pvalNVIQ, pvalGender, pvalHandedness),
                       ncol = 5)
        colnames(pvals1_3GH) <-
                c("RMSD.PRE.censoring",
                  "Age",
                  "WASI.NVIQ",
                  "Handedness",
                  "Gender")
        rownames(pvals1_3GH) <- "1-3GH"
        pvals1_3GH <- as.table(pvals1_3GH)
        print(pvals1_3GH)
        
        smdMotion <- smd(dataNewFilt, RMSD.PRE.censoring)
        smdAge <- smd(dataNewFilt, Age)
        smdNVIQ <- smd(dataNewFilt, WASI.NVIQ)
        smdGender <- smd(dataNewFilt, Gender)
        smdhandedness <- smd(dataNewFilt, Handedness)
        
        # Tabling smd for output
        smd1_3GH <-
                matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
                       ncol = 5)
        colnames(smd1_3GH) <-
                c("RMSD.PRE.censoring",
                  "Age",
                  "WASI.NVIQ",
                  "Handedness",
                  "Gender")
        rownames(smd1_3GH) <- "1_3GH"
        smd1_3GH <- as.table(smd1_3GH)
        print(smd1_3GH)
        print(table(dataNewFilt$group))
        }
        
}
