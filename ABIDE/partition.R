#Main function to split data
partition<-function(vars, y, subset, search, method, split, nsplit,
                    minsplit, minbucket, a, scale.y, useSearch, useOptim, useRpart, allVars){
  
  if(sum(subset) < minsplit){return(NULL)}
  ###make an if condition that if variable is categorical use whatever is doing in findsplit and else
  ### do whatever in partition
  
  vars<-vars[subset,,drop=FALSE]
  y<-y[subset]
  
  
  if(NROW(vars)<2*minbucket){stop("Can't split tree to satisfy minbucket")}
  #If all y values are the same, rpart will give an error message.  Do not split tree
  if(length(unique(y))==1){return(NULL)}
  
  nVars<-NCOL(vars)
  stats<-rep(NA,nVars)
  cutoff<-rep(NA,nVars)
  breakLeft<-vector("list",nVars)
  
  
  
  for(v in 1:nVars){
    #If randomly picking a subset of categories, do not sort by mean.  Would be more likely to select variables when sorted
    if(search=="greedy" & !is.null(nsplit)){xTemp<-ordinalize(x=vars[,v],y,sortCat=FALSE)
    } else {xTemp<-ordinalize(x=vars[,v],y)}
    x<-xTemp$x
    #If all x values the same, do not check optimal split
    if(abs(max(x) - min(x)) > 1e-8){
      #The SSS partition deals with problems when there is a very small number of observations
      #Use greedy search in this case (or set minsplit >= 5)
      if (search=="greedy"){
        #Note: Rpart uses impurity gain
        if(useRpart==FALSE){
          nX<-length(x)
          #Due to ties, it's possible minbucket cannot be satisfied
          if(sort(x)[minbucket]==sort(x,decreasing=TRUE)[minbucket]){cutpts=NULL
          } else {
            #Cutpoints considered must satisfy minbucket
            cutpts<-unique(sort(x)[minbucket:(nX-minbucket+1)])
            #Take average of distinct points to determine cutoff (like rpart)
            if(length(cutpts)==1){stop(paste0("Only 1 cutpt:", cutpts, x))}
            cutpts<-(cutpts[1:(length(cutpts)-1)]+cutpts[2:length(cutpts)])/2
          }
          #It is possible (unlikely) no cutpoint can satisfy minbucket
          if(!is.null(cutpts)){
            mod<-splitrule(y=y,x=x,cutpts=cutpts,method=method,split=split)
            stats[v] <- mod$stat
            if(is.factor(vars[,v])){
              breakLeft[[v]] <- rep(NA, length(levels(vars[,v])))
              breakLeft[[v]][levels(vars[,v]) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl <= mod$cutoff]]=1L
              breakLeft[[v]][levels(vars[,v]) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl > mod$cutoff]]=2L
              if(all(is.na(breakLeft[[v]])) & length(unique(breakLeft[[v]]))<=1){stop("Did not find correct cutpoints")}
            } else {cutoff[v] <- mod$cutoff}
          }
        } else { #Can use rpart to do greedy search (can't do nsplit)
          x<-vars[,v]
          #Note: use cp=-1 to enforce exactly one split (otherwise, may only have root node)
          #Use minsplit=2 (already checked minsplit) and maxcompete=1, maxsurrogate=0 to speed up convergence
          mod<-rpart(y~x, method=method, parms=list(split=split),
                     control=rpart.control(minsplit=2, minbucket=minbucket, cp=-1, xval=0, maxdepth=1,
                                           maxcompete=1, maxsurrogate=0))
          #It is possible (very unlikely) that no single split can satisfy minbucket constraint
          if(!is.null(mod$splits)){
            stats[v] <- -mod$splits[3]
            if(is.factor(vars[,v])){
              #csplit - 1 left, 2 not present, 3 right.  Change to 1, NA, 2 
              breakLeft[[v]] <- as.vector(mod$csplit)
              breakLeft[[v]][mod$csplit==2]=NA
              breakLeft[[v]][mod$csplit==3]=2L
            } else {cutoff[v] <- mod$splits[4]}
          }
        }
      } else {stop("Unexpected search")}
    }
  }
  #If each candidate variable cannot be split (e.g. cannot satisfy minbucket), return null
  if(all(is.na(stats))){return(NULL)}
  if(is.na(cutoff[which.min(stats)])){
    #Index is used for categorical variable splits
    return(partysplit(varid=as.integer(colnames(vars)[which.min(stats)]),
                      index=breakLeft[[which.min(stats)]],
                      info=list(stats=stats)))
  } else {
    #Breaks is used for continuous variable splits
    return(partysplit(varid=as.integer(colnames(vars)[which.min(stats)]),
                      breaks=cutoff[which.min(stats)],
                      info=list(stats=stats)))
  }
}