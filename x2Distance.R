#==============================================================================
# Copyright statement comment
# Author: Afrooz Jahedi
# Goal: Create groupwise (1-3) matching.
# Inputs: data, variables to match on, number of trees, method of matching
# Outputs:pvalues and standardized mean difference for each variable.number 
#  of subjects in each group
# Description: build a distance matrix based on random forest where subject distance  
# in the same terminal node is zero and subjects in different terminal nodes have
#  distance equal to 1-chi-square pvalue.
#==============================================================================

Gmatch <- function (data, matchList, nTree, method) {
  #==== Cleaning data ====
  # It is important to have ASD first and then TD for distance matrix
  #?????????????????????????? how to not hard code variables.
  data$group <- as.factor (data$group)
  
  # Number of subjects in each group for original data
  table(data$group)
  
  # screening for outlier in motion variable
  outlierVal <-
    mean(data$RMSD.PRE.censoring) + 3 * sd(data$RMSD.PRE.censoring)
  
  # Create the dataframe for matching without outlier and label rows
  Gdata <-
    data[which(data$RMSD.PRE.censoring <= outlierVal), matchList]
  rownames(Gdata) <- data$subj
  rownames(Gdata) <- data[, "subj"]
  
  # Number of subjects in each group for filtered data
  table (Gdata$group)
  
  # Outlier Subjects
  outliers <-
    as.character(data[which(data$RMSD.PRE.censoring >= outlierVal), ]$subj)
  cat("subjects with extreme motion are", outliers, "\n")
  
  
  #==== Tree building ====
  # Building the tree using covariates using greedy search for binary response
  #  which uses gini index for split. stop criteria for splitting is thirthy
  # subjects. At each split we use randomly three variables.
  #????????????????Hard coding variables
  form <-
    group ~ Gender + Handedness + RMSD.PRE.censoring + Age + WASI.NVIQ
  mydata <- grow(
    form,
    data = Gdata ,
    search = "greedy" ,
    method = "class" ,
    split = "gini" ,
    minsplit = 30,
    mtry = 3,
    nsplit = NULL
  )
  # Plot the tree
  plot(mydata)
  
  #==== variable pvalue before matching ====
  #????????????????hard coding variables
  pvalMotion <-
    (t.test(RMSD.PRE.censoring ~ group, data = Gdata))$p.value
  pvalAge <- (t.test(Age ~ group, data = Gdata))$p.value
  pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = Gdata))$p.value
  pvalGender <- chisq.test(Gdata$group, Gdata$Gender)$p.value
  pvalHandedness <-
    chisq.test(Gdata$group, Gdata$Handedness)$p.value
  
  # Tabling pvalues for output
  pvalsBef <-
    matrix(c(pvalMotion, pvalAge, pvalNVIQ, pvalGender, pvalHandedness),
           ncol = 5)
  colnames(pvalsBef) <-
    c("RMSD.PRE.censoring",
      "Age",
      "WASI.NVIQ",
      "Handedness",
      "Gender")
  rownames(pvalsBef) <- "BEFORE"
  pvalsBef <- as.table(pvalsBef)
  
  # Calculate standardized mean difference
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
      minsplit = 30,
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
  #=================== Creating chi-sqr distance matrix ================================
  # Extracting terminal node for subjects across all trees.
  nodeResponse <- sapply(rfRF100, function(x) {
    predict(x$tree, newdata = Gdata , type = "node")
  })
  # Preparing distance matrix 
  sumXDistMat = matrix(0, nrow = NROW(nodeResponse),ncol = NROW(nodeResponse))
  colnames(sumXDistMat) <- rownames(nodeResponse)
  rownames(sumXDistMat) <- rownames(nodeResponse)
  
  ntrees <- nTree
  
  # Define a function x2Dist
  #x2Dist(rfRF100[[treeNum]]$tree,nodeResponse,treeNum)
  
  # Add up tree p-values 
  for (treeNum in 1:ntrees){
    sumXDistMat <- sumXDistMat + x2Dist(rfRF100[[treeNum]]$tree,nodeResponse,treeNum)
  }
  
  #average the distance matrix and label rows and columns.
  xDistFor <- sumXDistMat / ntrees
  rownames(xDistFor) <- rownames(Gdata)
  colnames(xDistFor) <- rownames(Gdata)
  
  #==== percentile info from distance matrix =====
  quantile(xDistFor, c(0.25, 0.35, 0.45, 0.55, 0.6, 0.66, 0.7, 0.8))
  #      25%       35%       45%       55%       60%       66%       70%       80%
  # 0.6075259 0.6921892 0.7658434 0.8012043 0.8428033 0.8804234 0.8930838 0.9271538
  
  # percentile for just ASD vs TD distance matrix
  selXDistFor <- xDistFor[1:nMin, (nMin + 1):NROW(Gdata)]
  quantile(selXDistFor, c(0.25, 0.33, 0.35, 0.45, 0.55, 0.6, 0.66, 0.7, 0.8))
  #   25%       33%       35%       45%       55%       60%       66%       70%       80%
  #0.6351190 0.7066034 0.7289825 0.7766001 0.8459433 0.8761321 0.8987989 0.9112902 0.9480292
  
  
  #==== 1To3GH ====
  # 1- 1To3GH: 1-3 filtering subject x2-distance for Gender and handedness.
  #           Based on distance percentile, certain threshold will be picked.
  if(method =="1To3GH") {
    #Tight matching on gender and handedness using certain distance thershold
    
    # ASD and TD distance only and label them.
    x2distForSel <- xDistFor[1:nMin, (nMin + 1):ncol(xDistFor)]
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
    colnames(filtGen) <- rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
    rownames(filtGen) <- rownames(Gdata[1:nMin, ])
    #filter TD handedness that are equal to ASD handedness label to 1 otherwise 0
    filtHand <- foreach (i = 1:nMin) %do%
      ifelse(Gdata[rownames(x2distForSel)[i],
                   "Handedness"] == Gdata[colnames(x2distForSel), "Handedness"] , 1, 0)
    
    filtHand <-
      matrix(unlist(filtHand),
             ncol = ncol(x2distForSel),
             byrow = T)
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
        100
      ),
      nrow = nMin ,
      byrow = F
    )
    rownames(filtGenHandDist) <- rownames(Gdata[1:nMin, ])
    colnames(filtGenHandDist) <-
      rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
    
    # Create distance matrix csv file.
    outFileName = paste("output//x2DistMat", distFiltVal, ".csv", sep = "")
    write.table(filtGenHandDist,
                outFileName,
                sep = ",",
                col.names = TRUE)
    
    #=== what????? ====
    remFilt = filtGenHandDist[!apply(filtGenHandDist, 1, function(x) {
      all(x == 100)
    }),!apply(filtGenHandDist, 2, function(x) {
      all(x == 100)
    })]
    remRowName <- rownames(remFilt)
    remColName <- colnames(remFilt)
    
    remained <- x2distForSel[remRowName, remColName]
    
    #outFileName = paste("output//x2distance matrix 0.4GH", ".csv", sep = "")
    #write.table(remained, outFileName, sep = ",", col.names = TRUE)
    
    # Select 3 candidates for group-wise matching.
    selCandid <- 3
    
    # candidate values
    candidValFilt <- apply(remained, 1, function(x)
      return((sort(x))[1:selCandid]))
    
    # who are the TD selected candidates. not TD labels just their index in xDistFor mat.
    candidFilt <- apply(remained, 1, function(x)
      return(((order(
        x
      ))[1:selCandid])))
    
    # Vector of TD candidates and their index
    distFilt <- (table(candidFilt))
    namesTD <-
      colnames(remained[, as.numeric(names(table(candidFilt)))])
    
    # Shows the distribution of number of times each TD is selected
    hist(distFilt, main = "Number of times a subject is selected")
    dataNewExaX2 <- rbind (Gdata[remRowName, ], (Gdata[namesTD, ]))
  
  pvalHandedness <- chisq.test(dataNewExaX2$group, dataNewExaX2$Handedness)$p.value
  pvalGender <- chisq.test(dataNewExaX2$group, dataNewExaX2$Gender)$p.value
  pvalMotion <- (t.test(RMSD.PRE.censoring ~ group, data = dataNewExaX2))$p.value
  pvalAge <- (t.test(Age ~ group, data = dataNewExaX2))$p.value
  pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = dataNewExaX2))$p.value
  
  # Tabling pvalues for output
  pvals1_GH <-
    matrix(c(pvalMotion, pvalAge, pvalNVIQ, pvalGender, pvalHandedness),
           ncol = 5)
  colnames(pvals1-3GH) <-
    c("RMSD.PRE.censoring",
      "Age",
      "WASI.NVIQ",
      "Handedness",
      "Gender")
  rownames(pvals1_3GH) <- "1-3GH"
  pvals1_3GH <- as.table(pvals1_3GH)
  
  smdMotion <- smd(dataNewExaX2, RMSD.PRE.censoring)
  smdAge <- smd(dataNewExaX2, Age)
  smdNVIQ <- smd(dataNewExaX2, WASI.NVIQ)
  smdGender <- smd(dataNewExaX2, Gender)
  smdhandedness <- smd(dataNewExaX2, Handedness)
  
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
  print(table(dataNewExaX2$group))
  
  } else if (method == "1To3Dist") {
    ### x2-defined distance using 1-3 match
    # 2- 1To3Dist: 1 to 3 based on x2-distance matrix
    nodeResponse <- sapply(rfRF100, function(x) {
      predict(x$tree, newdata = Gdata , type = "node")
    })
    
    #xDistFor <- sumXDistMat/ntrees
    #rownames(xDistFor) <- rownames(Gdata)
    #colnames(xDistFor) <- rownames(Gdata)
    #mean(xDistFor)
    
    candidVal <-
      apply(xDistFor[1:nMin , (nMin + 1):ncol(xDistFor)], 1, function(x)
        return((sort(x))[1:3]))
    
    candid <-
      apply(xDistFor[1:nMin, (nMin + 1):ncol(xDistFor)], 1, function(x)
        return((order(x))[1:3] + nMin))
  
  distCandid <- (table(candid))
  names(distCandid) <- rownames(Gdata[as.numeric(names(table(candid))),])
  hist(distCandid, main = "Number of times a subject is selected")
  
  # Matched data using definde ditance
  dataNew1to3X2 <-
    rbind (Gdata[1:nMin, ], (Gdata[rownames(distCandid), ]))
  
  # 1-3 matching #
  pvalHandedness <-
    chisq.test(dataNew1to3X2$group, dataNew1to3X2$Handedness)$p.value
  pvalGender <-
    chisq.test(dataNew1to3X2$group, dataNew1to3X2$Gender)$p.value
  pvalMotion <-
    (t.test(RMSD.PRE.censoring ~ group, data = dataNew1to3X2))$p.value
  pvalAge <- (t.test(Age ~ group, data = dataNew1to3X2))$p.value
  pvalNVIQ <-
    (t.test(WASI.NVIQ ~ group, data = dataNew1to3X2))$p.value
  
  
  smdMotion <- smd(dataNew1to3X2, RMSD.PRE.censoring)
  smdAge <- smd(dataNew1to3X2, Age)
  smdNVIQ <- smd(dataNew1to3X2, WASI.NVIQ)
  smdGender <- smd(dataNew1to3X2, Gender)
  smdhandedness <- smd(dataNew1to3X2, Handedness)
  
  print(table(dataNew1to3X2$group))
    }else if(method == "1To3Caliper"){
      
      # 3- 1To3Caliper: 1 to 3 matching filtering distance based on propensity score and 
      #    within those find smallest distance.
      
    #### Distance within calipers by the propensity score############
    #### Propenity Score ###
    # predict(rfRF100[[1]]$tree[[1]],type = "prob")
    # extracting prob of being ASD for each subject across all trees.
    #
    pro <- rowMeans (sapply(rfRF100, function(x) {
        predict(x$tree, newdata = Gdata , type = "prob")[, 2]}))
    
    #logistic propensity
    logitPro <- log(pro / (1 - pro))
    
    spread <- 1/2 * (sd(logitPro))
    # calculate distance matrix for propensity vector
    library(fields)
    distPro <- rdist(logitPro)
    # We just need the upper tarangular of distance matrix without diagonals.
    distPro[lower.tri(distPro)] <- NA
    diag(distPro) <- NA
    
    distForFilt <- xDistFor
    filtPro = with(data.frame(distPro), subset(data.frame(distPro)<= spread))
    distFodistrFilt <- ifelse(filtPro == FALSE, 100, distForFilt)
    
    # min value of distance forest (only TD distane) which is filtered by propensity across row.
    minValDFP <- apply(distFodistrFilt[1:nMin ,(nMin+1): ncol(distForFilt)], 1, function(x)
        return((sort(x))[1:3]))
    
    minDFP <- apply(distFodistrFilt[ 1:nMin, (nMin+1): ncol(distForFilt)], 1, function(x)
        return((order(x))[1:3]+nMin))
    
    # Select TD subject equal to the number of ASD subjects from the table pool of selected candidates. 
    DFPCand <- sort(table(minDFP), decreasing = TRUE)
    
    names(DFPCand) <- rownames(Gdata[as.numeric(names(table(minDFP))),])
    
    # Numbers of times TD candidates can be selected histogram.
    hist(DFPCand, main = "Number of times a subject is selected")
    
    # number of unique TD candidates
    length((which(!is.na(unique(rownames(DFPCand))))))
    
    # Matched data using within caliper propensity score.
    data1To3Caliper <- rbind (Gdata[1:nMin, ], (Gdata [rownames(DFPCand), ]))

    #### p-values after distance within the calipar matching ###
    chisq.test(data1To3Caliper$group, data1To3Caliper$Handedness)
    chisq.test(data1To3Caliper$group, data1To3Caliper$Gender)
    pvalMotion <-(t.test(RMSD.PRE.censoring ~ group, data = data1To3Caliper))$p.value
    pvalAge <- (t.test(Age ~ group, data = data1To3Caliper))$p.value
    pvalNVIQ <-(t.test(WASI.NVIQ ~ group, data = data1To3Caliper))$p.value
    
    DDFP <- printPval(pvalMotion, "data1To3Caliper$RMSD.PRE.censoring")
    EDFP <- printPval(pvalAge, "data1To3Caliper$Age")
    GDFP <- printPval(pvalNVIQ, "data1To3Caliper$WASI.NVIQ")
    JDFP <- printPval(pvalHandedness, "data1To3Caliper$handedness")
    KDFP <- printPval(pvalGender, "data1To3Caliper$Gender")
    
    smdMotion <- smd(data1To3Caliper, RMSD.PRE.censoring)
    smdAge <- smd(data1To3Caliper, Age)
    smdNVIQ <- smd(data1To3Caliper, WASI.NVIQ)
    smdGender <- smd(data1To3Caliper, Gender)
    smdhandedness <- smd(data1To3Caliper, Handedness)
    print ( table(data1To3Caliper$group))
    
  }
  else if (method =="inverse"){
invDistFor <-xDistFor
# We just need the upper tarangular of inverse distance matrix without diagonals.
invDistFor[lower.tri(invDistFor)] <- NA
diag(invDistFor) <- NA
sumInv <-
  apply(invDistFor[1:nMin , (nMin + 1):ncol(invDistFor)], 2, function(x)
    return((sum(sort(x, decreasing = F)[1:3]))))
minValInv <- sort(sumInv, decreasing = F)[1:(nMin)]
candidInv <- (order(sumInv, decreasing = F) + nMin)[1:(nMin)]
#dataNewInv <- rbind (Gdata[1:nMin,], Gdata[candidInv,])
dataNewInv <- rbind(Gdata[1:(nMin),], Gdata[candidInv,])

#### p-values after inverse distance matching ###
chisq.test(dataNewInv$group, dataNewInv$Handedness)
chisq.test(dataNewInv$group, dataNewInv$Gender)
pvalMotion <-(t.test(RMSD.PRE.censoring ~ group, data = dataNewInv))$p.value
pvalAge <- (t.test(Age ~ group, data = dataNewInv))$p.value
pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = dataNewInv))$p.value

DDFPB <- printPval(pvalMotion, "dataNewInv$RMSD.PRE.censoring")
EDFPB <- printPval(pvalAge, "dataNewInv$Age")
GDFPB <- printPval(pvalNVIQ, "dataNewInv$WASI.NVIQ")
JDFPB <- printPval(pvalHandedness, "dataNewInv$handedness")
KDFPB <- printPval(pvalGender, "dataNewInv$Gender")

smdMotion <- smd(dataNewInv, RMSD.PRE.censoring)
smdAge <- smd(dataNewInv, Age)
smdNVIQ <- smd(dataNewInv, WASI.NVIQ)
smdGender <- smd(dataNewInv, Gender)
smdhandedness <- smd(dataNewInv, Handedness)
table(dataNewInv$group)

pvalsInv <-
  matrix(c(pvalMotion, pvalAge, pvalNVIQ, pvalGender, pvalHandedness),
         ncol = 5)
colnames(pvalsInv) <-
  c("RMSD.PRE.censoring",
    "Age",
    "WASI.NVIQ",
    "Handedness",
    "Gender")
rownames(pvalsInv) <- "Inverse"
pvalsInv <- as.table(pvalsInv)
pvalsInv
rbind(pvals,pvalsInv)
  }
}
