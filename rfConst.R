#Grows a random forest
rfConst<-function(ntrees,formula,training,search,method,split=c("information","ttest","MSE","entropy", "gini"),
                  mtry,nsplit,minsplit,minbucket,maxdepth,bootstrap=TRUE,useRpart=FALSE){
  
  #Construct random forest
  randFor<-lapply(1:ntrees, function(b){
    if(b%%100==0){print(paste0("Tree Number: ",b))}
    if(bootstrap==FALSE){
      obs.b<-1:nrow(training)
      sample.b<-training
    } else{
      obs.b <- sample(1:nrow(training), size=nrow(training), replace=T)
      sample.b <- training[obs.b,]
    }
    tree.b<-grow(formula=formula, data=sample.b, search=search, method=method, split=split,
                 mtry=mtry, nsplit=nsplit, minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth, 
                 useRpart=useRpart)
    list(tree=tree.b,cases=sort(unique(obs.b)))
  })
  
  return(randFor)
}
