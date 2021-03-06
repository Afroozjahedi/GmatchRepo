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

Gmatch <- function (data, matchList, nTree) {
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
  scatter3D(
    x,
    y,
    z,
    phi = 0,
    bty = "g",
    type = "h",
    ticktype = "detailed",
    pch = 19,
    cex = 0.5,
    xlab = "NVIQ",
    ylab = "Age",
    zlab = "Motion",
    main = "BRIEF data before matching",
    clab = c("Motion")
  )
  # pvalue before matching ###
  pvalMotion <-
    (t.test(RMSD.PRE.censoring ~ group, data = Gdata))$p.value
  pvalAge <- (t.test(Age ~ group, data = Gdata))$p.value
  pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = Gdata))$p.value
  pvalHandedness <-
    chisq.test(Gdata$group, Gdata$Handedness)$p.value
  pvalGender <- chisq.test(Gdata$group, Gdata$Gender)$p.value
  
  A <- printPval(pvalMotion, "Gdata$RMSD.PRE.censoring")
  B <- printPval(pvalAge, "Gdata$Age")
  C <- printPval(pvalNVIQ, "Gdata$WASI.NVIQ")
  H <- printPval(pvalHandedness, "Gdata$Handedness")
  I <- printPval(pvalGender, "Gdata$Gender")
  
  smdMotion <- smd(Gdata, RMSD.PRE.censoring)
  smdAge <- smd(Gdata, "Age")
  smdNVIQ <- smd(Gdata, "WASI.NVIQ")
  smdGender <- smd(Gdata, "Gender")
  smdHandedness <- smd(Gdata, "Handedness")
  
  # Random Forest growing ###
  set.seed(184684)
  ntrees <- nTree
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
  # return(rfRF100)
  
  ### Reference group for matching to start with ########
  # which group has less subject. That group will be our reference for matching.
  ASD <- Gdata$group == 1
  if (NROW (Gdata[ASD, ]) < NROW (Gdata[!ASD,])) {
    (nMin = NROW (Gdata[ASD,]))
  } else {
    nMin <- NROW (Gdata[!ASD,])
  }
  
  
  ####New deifined matrix
  #terminal node for subjects across all trees.
  nodeResponse <-
    sapply(rfRF100, function(x) {
      predict(x$tree, newdata = Gdata , type = "node")
    })
  plot(rfRF100[[1]]$tree)
  
  # Terminal nodes
  terNode <- nodeids(rfRF100[[1]]$tree, terminal = TRUE)
  
  # How many subjects in each terminal node
  NObs = sapply(terNode, function(n) {
    nrow(rfRF100[[1]]$tree[n]$data)
  })
  
  # Table subjects based on group for each terminal node
  groupObs = t(sapply(terNode, function(n) {
    table(rfRF100[[1]]$tree[n]$data$group)
  }))
  
  #make data frame with terminal node labels, #ASD, #TD in that terminal node
  tabNode <- cbind(terNode, NObs, groupObs)
  
  # chisq test for all pairs of terminal nodes
  xPval <- matrix(NA, NROW(tabNode), (NROW(tabNode)))
  for (ter1 in 1:NROW(tabNode)) {
    for (ter2 in 1:NROW(tabNode)) {
      xPval[ter1, ter2] <- (chisq.test(tabNode[ter1:ter2, 3:4])$p.value)
    }
  }
  
  diag(xPval) <- NA
  colnames(xPval) <- terNode
  rownames(xPval) <- terNode
  
  #setting up chi-square distnca matrix.
  xDistMat <-
    matrix(0,
           nrow = NROW(nodeResponse),
           ncol = NROW(nodeResponse))
  colnames(xDistMat) <- rownames(Gdata)
  rownames(xDistMat) <- rownames(Gdata)
  
  #find subjects that are in the same terminal node
  ntrees = 1
  for (tree in 1:ntrees) {
    for (row in 1:NROW(nodeResponse)) {
      xDistMat[which(nodeResponse[row, tree] == nodeResponse[, tree]), row] <-
        0
      a <- as.character(nodeResponse[row, tree])
      b <- as.character(nodeResponse[(nodeResponse[row, tree] != nodeResponse[, tree]), tree])
      c <- which(nodeResponse[row, tree] != nodeResponse[, tree])
      xDistMat[c, row]  <- xPval[a, b]
    }
  }
  
  
  #ASD subjects in terminal node 3 in tree 1.
  (rfRF100[[1]]$tree[[3]]$data$group == 1)
  
  #terminal node for subject 022A
  b = nodeResponse[1, 1]
  
  #Terminal nodes for subjects that are not in the same terminal node as 022A
  d22A = nodeResponse[which(nodeResponse[1, 1] != nodeResponse[, 1]), 1]
  
  #subjects that are in the same terminal nodes as 022A
  s022A = which(nodeResponse[1, 1] == nodeResponse[, 1])
  
  #how many subjects are in terminal node 2- object b?
  table(nodeResponse[, 1] == b)[2]
  
  #How many subjects are in terminal node that a specific subject lays-for example 038A?
  table(nodeResponse[, 1] == d22A[1])[2]
  
  # what is their original group for subjects that are in terminal 2?
  Gdata[nodeResponse[, 1] == b, "group"]
  
  # what is their original group for subjects that are in terminal 11?
  Gdata[nodeResponse[, 1] == d22A[1], "group"]
  
  #put them together
  Gdata["node"] <- ifelse(nodeResponse[, 1] == b, 1, 0)
  
  #data that are in both terminal 2 and 11
  D = (Gdata[nodeResponse[, 1] == d22A[1] | nodeResponse[, 1] == b, ])
  
  #1-pvalue of contigency table for two different subjects classified by group.
  d = 1 - chisq.test(D$group, D$nodeI)$p.value
  
  #Setup chi-square distance matrix
  xdistanMat <- cematrix(0, nrow = nrow(Gdata), ncol = nrow(Gdata))
  
  
  
  
  #set same terminal node to zero
  ifelse(nodeResponse[1, 1] == nodeResponse[, 1], 0, NA)
  for (row in 1:NROW(nodeResponse)) {
    distMat[, row] <-
      distMat[, row] - (nodeResponse[1, 1] == nodeResponse[, 1])
  
    
  }
  
  
  #### Disstance matrix ####################
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
  mean(distForest)
  rownames(distForest) <- rownames(Gdata)
  colnames(distForest) <- rownames(Gdata)
  
  
  candidVal <-
    apply(distForest[1:nMin , (nMin + 1):ncol(distForest)], 1, function(x)
      return((sort(x))[1:3]))
  
  candid <-
    apply(distForest[1:nMin, (nMin + 1):ncol(distForest)], 1, function(x)
      return((order(x))[1:3] + nMin))
  #candid <- NULL
  # nMin is the nember of ASD subject. Choose the 3  closest candidate from distance Forest for each ASD subject
  #for (i in 1:nMin) {
  #  candid <-
  #    c(candid, (order(distForest[i,(nMin + 1):nrow(nodeResponse)])[1:3]) +
  #        nMin)
  #}
  distCandid <- (table(candid))
  names(distCandid) <-
    rownames(Gdata[as.numeric(names(table(candid))),])
  hist(distCandid, main = "Number of times a subject is selected")
  # heatmap of distance forest to see if there is any evident pattern between
  # ASD and TD.
  # (image(t(distForest[nrow(distForest):1,]), axes = FALSE,zlim = c(-4, 4),
  #       col = rainbow(5)
  #))
  
  # Check the mean and std of Forest distances in ASD vs TD
  #sel = distForest[1:table(Gdata$group)["1"],
  #                 table(Gdata$group)["1"] + 1:(table(Gdata$group)["0"]+
  #                 table(Gdata$group)["1"])
  
  #(image(t(sel[nrow(sel):1,]), axes = FALSE, zlim = c(-4, 4), col = rainbow(54)
  #))
  #mean(sel)
  # clusters <- kmeans(as.matrix(Gdata),10, nstart = 20)
  # clusters
  # names(clusters)
  # clusters$size
  #sd(sel)
  # Matched data using definde ditance
  dataNew <-
    rbind (Gdata[1:nMin,], (Gdata[rownames(distCandid),]))
  # PCA of the new dataset
  # p <- prcomp(dataNew[ ,c(-1, -2, -3)],scale. = T); p;
  # summary(p);biplot(p);plot(p,type="l");
  
  #### distance Matrix p-values after matching for  ####
  pvalHandedness <-
    chisq.test(dataNew$group, dataNew$Handedness)$p.value
  pvalGender <- chisq.test(dataNew$group, dataNew$Gender)$p.value
  pvalMotion <-
    (t.test(RMSD.PRE.censoring ~ group, data = dataNew))$p.value
  pvalAge <- (t.test(Age ~ group, data = dataNew))$p.value
  pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = dataNew))$p.value
  
  D <- printPval(pvalMotion, "dataNew$RMSD.PRE.censoring")
  E <- printPval(pvalAge, "dataNew$Age")
  G <- printPval(pvalNVIQ, "dataNew$WASI.NVIQ")
  J <- printPval(pvalHandedness, "dataNew$handedness")
  K <- printPval(pvalGender, "dataNew$Gender")
  
  smdMotion <- smd(dataNew, RMSD.PRE.censoring)
  smdAge <- smd(dataNew, Age)
  smdNVIQ <- smd(dataNew, WASI.NVIQ)
  smdGender <- smd(dataNew, Gender)
  smdhandedness <- smd(dataNew, Handedness)
  
  #### Disstance matrix deleting trouble makers##############
  
  distForSel <- distForest[1:nMin, (nMin + 1):ncol(distForest)]
  rownames(distForSel) <- rownames(Gdata[1:nMin, ])
  colnames(distForSel) <-
    rownames(Gdata[(nMin + 1):ncol(distForest), ])
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
  colnames(filtGen) <- rownames(Gdata[(nMin + 1):ncol(distForest), ])
  rownames(filtGen) <- rownames(Gdata[1:nMin, ])
  
  filtHand <- foreach (i = 1:nMin) %do%
    ifelse(Gdata[rownames(distForSel)[i],
                 "Handedness"] == Gdata[colnames(distForSel), "Handedness"] , 1, 0)
  filtHand <-
    matrix(unlist(filtHand),
           ncol = ncol(distForSel),
           byrow = T)
  colnames(filtHand) <- rownames(Gdata[(nMin + 1):ncol(distForest), ])
  rownames(filtHand) <- rownames(Gdata[1:nMin, ])
  
  filtGenHand <- filtGen + filtHand
  rownames(filtGenHand) <- rownames(Gdata[1:nMin, ])
  colnames(filtGenHand) <-
    rownames(Gdata[(nMin + 1):ncol(distForest),])
  
  distFiltVal <- 0.4
  filtGenHandDist <-
    matrix(
      ifelse(distForSel[,] <= distFiltVal & filtGenHand[,] == 2,
             distForSel, 100),
      nrow = nMin ,
      byrow = F
    )
  rownames(filtGenHandDist) <- rownames(Gdata[1:nMin,])
  colnames(filtGenHandDist) <-
    rownames(Gdata[(nMin + 1):ncol(distForest),])
  
  #outFileName = paste("output//Mat", distFiltVal, ".csv", sep="")
  #write.table(filtGenHandDist, outFileName, sep=",", col.names=TRUE)
  
  #ASD trouble makers
  #table((rowSums(filtGenHandDist) >= 5200))
  # Useless TD's
  #table(colSums(filtGenHandDist) >= 4700)
  
  remFilt = filtGenHandDist[!apply(filtGenHandDist, 1, function(x) {
    all(x == 100)
  }),!apply(filtGenHandDist, 2, function(x) {
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
  dataNewFilt <- rbind (Gdata[remRowName,], (Gdata[namesTD,]))
  
  ### distance Matrix p-values after matching Matched data deleting trouble ASD
  
  pvalHandedness <-
    chisq.test(dataNewFilt$group, dataNewFilt$Handedness)$p.value
  pvalGender <-
    chisq.test(dataNewFilt$group, dataNewFilt$Gender)$p.value
  pvalMotion <-
    (t.test(RMSD.PRE.censoring ~ group, data = dataNewFilt))$p.value
  pvalAge <- (t.test(Age ~ group, data = dataNewFilt))$p.value
  pvalNVIQ <-
    (t.test(WASI.NVIQ ~ group, data = dataNewFilt))$p.value
  
  DDF <- printPval(pvalMotion, "dataNewFilt$RMSD.PRE.censoring")
  EDF <- printPval(pvalAge, "dataNewFilt$Age")
  GDF <- printPval(pvalNVIQ, "dataNewFilt$WASI.NVIQ")
  JDF <- printPval(pvalHandedness, "dataNewFilt$handedness")
  KDF <- printPval(pvalGender, "dataNewFilt$Gender")
  
  smdMotion <- smd(dataNewFilt, RMSD.PRE.censoring)
  smdAge <- smd(dataNewFilt, Age)
  smdNVIQ <- smd(dataNewFilt, WASI.NVIQ)
  smdGender <- smd(dataNewFilt, Gender)
  smdhandedness <- smd(dataNewFilt, Handedness)
  
  
  
  
  ############## Output of the Gmatch
  print (list(
    table(data$group),
    table(dataNew$group),
    table(dataNewFilt$group)
  ))
  #cat("subjects with extreme motion are" , outliers)
  
}
