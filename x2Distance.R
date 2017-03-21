###
# Copyright statement comment
# Author comment:
# File description comment, including purpose of program, inputs, and outputs
# new defined distance matrix.distance is zero if two subjects are in the same
# terminal node, 1-pvalue if different terminal node.
# source() and library() statements
# Function definitions
# Executed statements, if applicable (e.g., print, plot)
###

Gmatch <- function (data, matchList, nTree, method) {
  # select best candidate from a group(ASD/TD) to be matched for (TD/ASD).
  #
  # Args:
  #   data: unmatched data with no missing values.
  #   matchList: list of variables that are meant to be matched on.,
  #
  #   nTree: Number of trees to build the forest.
  #
  # Returns:
  #   The matched data and p-value and SMD for each variable
  
  # Cleaning data ####
  # it is important to have ASD first and then TD for this code
  data$group <- as.factor (data$group)
  table(data$group)
  
  # screen for outlier in motion variable
  
  outlierVal <-
    mean(data$RMSD.PRE.censoring) + 3 * sd(data$RMSD.PRE.censoring)
  Gdata <-
    data[which(data$RMSD.PRE.censoring <= outlierVal), matchList]
  rownames(Gdata) <- data$subj
  rownames(Gdata) <- data[, "subj"]
  table (Gdata$group)
  #outliers <- as.character(data[which(data$RMSD.PRE.censoring >= outlierVal),]$subj)
  #cat("subjects with extreme motion are",outliers,"\n")
  
  # Gdata = data.frame(data[,matchList[1]] ,data[,matchList[2]], data[,matchList[3]],
  #                   data[,matchList[4]], data[,matchList[5]],  data[,matchList[6]])
  #colnames(Gdata) = varList
  # data_real = data.frame(data[1],group ,Gender, Handedness, RMSD.PRE.censoring, Age, WASI.NVIQ)
  # Tree building ###
  form <-
    group ~ Gender + Handedness + RMSD.PRE.censoring + Age + WASI.NVIQ
  (
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
  )
  
  plot(mydata)
  
  library("plot3D")
  x <- Gdata$WASI.NVIQ
  y <- Gdata$Age
  z <- Gdata$RMSD.PRE.censoring
  scatter3D(  x,  y,  z,  phi = 0,bty = "g",type = "h",
    ticktype = "detailed",  pch = 19,  cex = 0.5, xlab = "NVIQ",
    ylab = "Age",   zlab = "Motion",main = "BRIEF data before matching",
    clab = c("Motion"))
  
  # pvalue before matching ###
  pvalMotion <-(t.test(RMSD.PRE.censoring ~ group, data = Gdata))$p.value
  pvalAge <- (t.test(Age ~ group, data = Gdata))$p.value
  pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = Gdata))$p.value
  pvalHandedness <-chisq.test(Gdata$group, Gdata$Handedness)$p.value
  pvalGender <- chisq.test(Gdata$group, Gdata$Gender)$p.value
  
  A <- printPval(pvalMotion, "Gdata$RMSD.PRE.censoring")
  B <- printPval(pvalAge, "Gdata$Age")
  C <- printPval(pvalNVIQ, "Gdata$WASI.NVIQ")
  H <- printPval(pvalHandedness, "Gdata$Handedness")
  I <- printPval(pvalGender, "Gdata$Gender")
  
 pvals <- matrix(c(A,B,C,H,I),ncol=5)
   colnames(pvals) <- c("RMSD.PRE.censoring","Age","WASI.NVIQ","Handedness", "Gender")
   rownames(pvals) <- c("ntree")
   pvals <- as.table(pvals)
   pvals
 
  smdMotion <- smd(Gdata, RMSD.PRE.censoring)
  smdAge <- smd(Gdata, "Age")
  smdNVIQ <- smd(Gdata, "WASI.NVIQ")
  smdGender <- smd(Gdata, "Gender")
  smdHandedness <- smd(Gdata, "Handedness")
  
  ### Random Forest growing ###
  set.seed(184684)
  ntrees <- nTree
  rfRF100 <- rfConst(ntrees = ntrees,formula = form, training =  Gdata,
      search = "greedy",  method = "class", split = "gini", mtry = 3,
      nsplit = NULL,  minsplit = 30,maxdepth = 10, minbucket = 10,
      bootstrap = FALSE    )
  # return(rfRF100)
  
  ### Reference group for matching to start with ########
  # which group has less subject. That group will be our reference for matching.
  ASD <- Gdata$group == 1
  if (NROW (Gdata[ASD, ]) < NROW (Gdata[!ASD,])) {
    (nMin = NROW (Gdata[ASD,]))
  } else {
    nMin <- NROW (Gdata[!ASD,])
  }
  #terminal node for subjects across all trees.
  nodeResponse <- sapply(rfRF100, function(x) {
    predict(x$tree, newdata = Gdata , type = "node")
  })
  
  #plot(rfRF100[[1]]$tree)
  #debug(x2Dist)
  #a= x2Dist(rfRF100[[1]]$tree,nodeResponse)
  
  sumXDistMat = matrix(0, nrow = NROW(nodeResponse),ncol = NROW(nodeResponse))
  colnames(sumXDistMat) <- rownames(nodeResponse)
  rownames(sumXDistMat) <- rownames(nodeResponse)
  
  #treeNum <- 2
  ntrees <- nTree
  #x2Dist(rfRF100[[treeNum]]$tree,nodeResponse,treeNum)
  
  for (treeNum in 1:ntrees){
    sumXDistMat <- sumXDistMat + x2Dist(rfRF100[[treeNum]]$tree,nodeResponse,treeNum)
  }
  xDistFor <- sumXDistMat/ntrees
  rownames(xDistFor) <- rownames(Gdata)
  colnames(xDistFor) <- rownames(Gdata)
  
  #Descriptive statistics for distance matrix
  quantile(xDistFor)
  #  0%           25%       50%       75%      100% 
  #0.0000000 0.6075259 0.7792734 0.9103751 0.9951918 
  quantile(xDistFor,c(0.25,0.35,0.45,0.55,0.6,0.66, 0.7,0.8))
  #      25%       35%       45%       55%       60%       66%       70%       80% 
  # 0.6075259 0.6921892 0.7658434 0.8012043 0.8428033 0.8804234 0.8930838 0.9271538 
  
  # Descriptive statistics for just ASD vs TD distance matrix
  selXDistFor <- xDistFor[1:nMin,(nMin+1):NROW(Gdata)]
  quantile(selXDistFor)
  #0%           25%       50%       75%      100% 
  #0.0000000 0.6351190 0.8039192 0.9262075 0.9950895 
  quantile(selXDistFor,c(0.25,0.33, 0.35,0.45,0.55,0.6,0.66, 0.7,0.8))
  #   25%       33%       35%       45%       55%       60%       66%       70%       80% 
  #0.6351190 0.7066034 0.7289825 0.7766001 0.8459433 0.8761321 0.8987989 0.9112902 0.9480292 

  #### New deifined matrix-x2distance matrix 
  #### x2-Disstance matrix deleting trouble makers and exact matching ###
  
  #there are different methods
  # 1- 1To3GH: 1 to 3 matching filtering subject x2-distance for Gender and handedness
  # 2- 1To3Dist: 1 to 3 based on x2-distance matrix
  # 3- 1To3Caliper: 1 to 3 matching filtering distance based on propensity score and 
  #    within those find smallest distance.
  
  if(method =="1To3GH"){
    
  x2distForSel <- xDistFor[1:nMin, (nMin + 1):ncol(xDistFor)]
  rownames(x2distForSel) <- rownames(Gdata[1:nMin, ])
  colnames(x2distForSel) <- rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
  
  #filter TD gender that are equal to ASD gender to 1 otherwise 0.
  library(foreach)
  filtGen <- foreach (i = 1:nMin) %do%
    ifelse(Gdata[rownames(x2distForSel)[i],
                 "Gender"] == Gdata[colnames(x2distForSel), "Gender"] , 1, 0)
  
  filtGen <-matrix(unlist(filtGen), ncol = ncol(x2distForSel), byrow = T)
  colnames(filtGen) <- rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
  rownames(filtGen) <- rownames(Gdata[1:nMin, ])
  
  filtHand <- foreach (i = 1:nMin) %do%
    ifelse(Gdata[rownames(x2distForSel)[i],
                 "Handedness"] == Gdata[colnames(x2distForSel), "Handedness"] , 1, 0)
  
  filtHand <- matrix(unlist(filtHand), ncol = ncol(x2distForSel), byrow = T)
  colnames(filtHand) <- rownames(Gdata[(nMin + 1):ncol(xDistFor), ])
  rownames(filtHand) <- rownames(Gdata[1:nMin, ])
  
  filtGenHand <- filtGen + filtHand
  rownames(filtGenHand) <- rownames(Gdata[1:nMin, ])
  colnames(filtGenHand) <- rownames(Gdata[(nMin + 1):ncol(xDistFor),])
  
  distFiltVal <- 0.7
  filtGenHandDist <- matrix(
      ifelse(x2distForSel[,] <= distFiltVal & filtGenHand[,] == 2,
             x2distForSel, 100),nrow = nMin , byrow = F)
  rownames(filtGenHandDist) <- rownames(Gdata[1:nMin,])
  colnames(filtGenHandDist) <- rownames(Gdata[(nMin + 1):ncol(xDistFor),])
  
  outFileName = paste("output//x2DistMat", distFiltVal, ".csv", sep="")
  write.table(filtGenHandDist, outFileName, sep=",", col.names=TRUE)
  
  #ASD trouble makers
  #table((rowSums(filtGenHandDist) >= 5200))
  # Useless TD's
  #table(colSums(filtGenHandDist) >= 4700)
  
  remFilt = filtGenHandDist[!apply(filtGenHandDist, 1, function(x) {
    all(x == 100) }),!apply(filtGenHandDist, 2, function(x) {all(x == 100)})]
  remRowName <- rownames(remFilt)
  remColName <- colnames(remFilt)
  
  remained <- x2distForSel[remRowName, remColName]
  
  #outFileName = paste("output//x2distance matrix 0.4GH", ".csv", sep = "")
  #write.table(remained, outFileName, sep = ",", col.names = TRUE)
  
  selCandid <- 3
  candidValFilt <- apply(remained, 1, function(x)
      return((sort(x))[1:selCandid]))
  
  # who are the TD selected candidates. not TD labels just their index in xDistFor mat.
  candidFilt <- apply(remained, 1, function(x)
      return(((order(x))[1:selCandid])))
  
  # Vector of TD candidates
  distFilt <- (table(candidFilt))
  namesTD <- colnames(remained[, as.numeric(names(table(candidFilt)))])
  hist(distFilt, main = "Number of times a subject is selected")
  
  ###  Matched data using defined ditance deleting trouble ASD
  dataNewExaX2 <- rbind (Gdata[remRowName,], (Gdata[namesTD,]))
  
  ### distance Matrix p-values after matching Matched data deleting trouble ASD
  pvalHandedness <- chisq.test(dataNewExaX2$group, dataNewExaX2$Handedness)$p.value
  pvalGender <- chisq.test(dataNewExaX2$group, dataNewExaX2$Gender)$p.value
  pvalMotion <- (t.test(RMSD.PRE.censoring ~ group, data = dataNewExaX2))$p.value
  pvalAge <- (t.test(Age ~ group, data = dataNewExaX2))$p.value
  pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = dataNewExaX2))$p.value
  
  DDF <- printPval(pvalMotion, "dataNewExaX2$RMSD.PRE.censoring")
  EDF <- printPval(pvalAge, "dataNewExaX2$Age")
  GDF <- printPval(pvalNVIQ, "dataNewExaX2$WASI.NVIQ")
  JDF <- printPval(pvalHandedness, "dataNewExaX2$handedness")
  KDF <- printPval(pvalGender, "dataNewExaX2$Gender")
  
  smdMotion <- smd(dataNewExaX2, RMSD.PRE.censoring)
  smdAge <- smd(dataNewExaX2, Age)
  smdNVIQ <- smd(dataNewExaX2, WASI.NVIQ)
  smdGender <- smd(dataNewExaX2, Gender)
  smdhandedness <- smd(dataNewExaX2, Handedness)
  print(table(dataNewExaX2$group))
  }else if(method=="1To3Dist"){
  
  ### x2-defined distance using 1-3 match ###
  nodeResponse <- sapply(rfRF100, function(x) {
      predict(x$tree, newdata = Gdata , type = "node")
    })
  
  xDistFor <- sumXDistMat/ntrees
  rownames(xDistFor) <- rownames(Gdata)
  colnames(xDistFor) <- rownames(Gdata)
  #mean(xDistFor)  
  
  candidVal <- apply(xDistFor[1:nMin , (nMin + 1):ncol(xDistFor)], 1, function(x)
      return((sort(x))[1:3]))
  
  candid <- apply(xDistFor[1:nMin, (nMin + 1):ncol(xDistFor)], 1, function(x)
      return((order(x))[1:3] + nMin))
  
  distCandid <- (table(candid))
  names(distCandid) <- rownames(Gdata[as.numeric(names(table(candid))),])
  hist(distCandid, main = "Number of times a subject is selected")
  
  # Matched data using definde ditance
  dataNew1to3X2 <- rbind (Gdata[1:nMin,], (Gdata[rownames(distCandid),]))
  
  #### distance Matrix p-values after matching for  ####
  pvalHandedness <- chisq.test(dataNew1to3X2$group, dataNew1to3X2$Handedness)$p.value
  pvalGender <- chisq.test(dataNew1to3X2$group, dataNew1to3X2$Gender)$p.value
  pvalMotion <- (t.test(RMSD.PRE.censoring ~ group, data = dataNew1to3X2))$p.value
  pvalAge <- (t.test(Age ~ group, data = dataNew1to3X2))$p.value
  pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = dataNew1to3X2))$p.value
  
  D <- printPval(pvalMotion, "dataNew1to3X2$RMSD.PRE.censoring")
  E <- printPval(pvalAge, "dataNew1to3X2$Age")
  G <- printPval(pvalNVIQ, "dataNew1to3X2$WASI.NVIQ")
  J <- printPval(pvalHandedness, "dataNew1to3X2$handedness")
  K <- printPval(pvalGender, "dataNew1to3X2$Gender")
  
  smdMotion <- smd(dataNew1to3X2, RMSD.PRE.censoring)
  smdAge <- smd(dataNew1to3X2, Age)
  smdNVIQ <- smd(dataNew1to3X2, WASI.NVIQ)
  smdGender <- smd(dataNew1to3X2, Gender)
  smdhandedness <- smd(dataNew1to3X2, Handedness)
  print(table(dataNew1to3X2$group))
    }else if(method == "1To3Caliper"){
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
  }
}
