#==============================================================================
# Copyright statement comment
# Author: Afrooz Jahedi
# Goal: Create manual (1-1) matching.
# Inputs: data, variables to match on, number of trees, method of matching
# Outputs:pvalues and standardized mean difference for each variable.number 
#  of subjects in each group
# Description: Using the distance matrix from chi-squre definition, in which 
# rows are ASD and TD are columns. Any method is fine. manuall matching was done.
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
  
  x2distForSel <- xDistFor[1:nMin, (nMin + 1):ncol(xDistFor)]
  rownames(x2distForSel) <- rownames(Gdata[1:nMin,])
  colnames(x2distForSel) <-
    rownames(Gdata[(nMin + 1):ncol(xDistFor),])
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
    distFiltVal <- 0.84
    
    # Select TD subjects tha t qualify these 3 conditions.
    filtGenHandDist <- matrix(
      ifelse(
        x2distForSel[,] <= distFiltVal & filtGenHand[,] == 2,
        x2distForSel,
        NA
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
    
  # ASD Subjects that are selected manually 1-1 matching
  ASD <-c("022A","038A","045A","056A","071A","082A","091A","092A","095A","096A","111A","115A",
          "123A","131A","136A","138A","141A","144A","147A","148A","150A","153A","155A","157A",
           "160A","162A","167A","174A","175A","191A","192A","194A","198A_V2",
          "205A","209A","224A","225A","226A","405A","406A_V2", "409A","410A","413A",
          "420A")
  # ASD Subjects that are selected manually 1-1 matching
  TD <- c("048C","049C","064C","065C","067C","078C","086C","090C","097C","098C","100C","101C","104C","106C","108C",   
           "109C","110C","119C","121C","122C","128C","133C","154C","165C","166C","168C","172C",
            "173C","176C","178C","179C","182C","183C","185C","186C","195C","197C","211C","218C","221C",   
           "402C","407C","415C","418C_V2")
  
  # Create dataset for matching
  oneTonedat <- Gdata[c(ASD, TD), ]
  
  #==== Calculating pvalues and standardized mean difference ====
  pvalHandedness <-
    chisq.test(oneTonedat$group, oneTonedat$Handedness)$p.value
  pvalGender <-
    chisq.test(oneTonedat$group, oneTonedat$Gender)$p.value
  pvalMotion <-
    (t.test(RMSD.PRE.censoring ~ group, data = oneTonedat))$p.value
  pvalAge <- (t.test(Age ~ group, data = oneTonedat))$p.value
  pvalNVIQ <- (t.test(WASI.NVIQ ~ group, data = oneTonedat))$p.value
  
  # Tabling pvalues for output
  pvalsoneTone <-
    matrix(c(pvalMotion, pvalAge, pvalNVIQ, pvalGender, pvalHandedness),
           ncol = 5)
  colnames(pvalsoneTone) <-
    c("RMSD.PRE.censoring",
      "Age",
      "WASI.NVIQ",
      "Handedness",
      "Gender")
  rownames(pvalsoneTone) <- "oneTone"
  pvalsoneTone <- as.table(pvalsoneTone)
  print(pvalsoneTone)
  smdMotion <- smd(oneTonedat, RMSD.PRE.censoring)
  smdAge <- smd(oneTonedat, Age)
  smdNVIQ <- smd(oneTonedat, WASI.NVIQ)
  smdGender <- smd(oneTonedat, Gender)
  smdhandedness <- smd(oneTonedat, Handedness)
  
  # Tabling smd for output
  smdoneTone <-
    matrix(c(smdMotion, smdAge, smdNVIQ, smdGender, smdHandedness),
           ncol = 5)
  colnames(smdoneTone) <-
    c("RMSD.PRE.censoring",
      "Age",
      "WASI.NVIQ",
      "Handedness",
      "Gender")
  rownames(smdoneTone) <- "oneTone"
  smdoneTone <- as.table(smdoneTone)
  print(smdoneTone)
  }
}
