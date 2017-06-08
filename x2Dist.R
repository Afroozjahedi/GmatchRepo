#==============================================================================
# Copyright statement comment
# Author: Afrooz Jahedi
# Goal:Calculate chi-square pvalue of a two by two table of subjects in different 
#      nodes split by group variable.
# Inputs: info from a tree, terminal node index, of trees in the forest.
# Outputs: distance matrix based on p-values 
# Description: 
#==============================================================================

#x2Dist(rfRF100[[treeNum]]$tree,nodeResponse,treeNum)
x2Dist <- function(tree, terminalNode, treeNum) {
  pars <- as.list(match.call())
  
  # Identify terminal nodes in the tree
  terNode <- nodeids(tree, terminal = TRUE)
  
  # How many subjects in each terminal node
  NObs = sapply(terNode, function(n) {
    nrow(tree[n]$data)
  })
  
  # Table subjects based on group for each terminal node
  groupObs = t(sapply(terNode, function(n) {
    table(tree[n]$data$group)
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
  
  #take out diagonal elements of terminal nodes (e.g.2,2) we don't need them.
  diag(xPval) <- NA
  colnames(xPval) <- terNode
  rownames(xPval) <- terNode
  
  #setting up chi-square distnca matrix.
  xDistMat <-
    matrix(0,
           nrow = NROW(terminalNode),
           ncol = NROW(terminalNode))
  colnames(xDistMat) <- rownames(terminalNode)
  rownames(xDistMat) <- rownames(terminalNode)
  
  #find subjects that are in the same terminal node
    for (row in 1:NROW(terminalNode)) {
      
      #if same ter node, dist=0
      xDistMat[which(terminalNode[row, treeNum] == terminalNode[, treeNum]), row] <-
        0
      ##if not same ter node, what is the ter node for first fix subj
      a <- as.character(terminalNode[row, treeNum])
      #what is the ter node for second non-fix subj
      b <- as.character(terminalNode[(terminalNode[row, treeNum] != terminalNode[, treeNum]), treeNum])
      #what is the subj ID for subjects that are not in same ter node
      c <- which(terminalNode[row, treeNum] != terminalNode[, treeNum])
      #fill out the distance matrix accordingly
      xDistMat[c, row]  <- 1-xPval[a, b]
    }
  return(xDistMat)
  
 # cat("terminal node","=", terNode,"\n", 
      #NObs,"\n",
      #groupObs,"\n",tabNode,"\n",xDistMat)
}