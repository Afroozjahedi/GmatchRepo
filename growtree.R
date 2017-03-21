#Grows tree by using partition() function several times
growtree<-function(id=1L, depth=1L, data, response, subset, search, method, split,
                   mtry, nsplit, minsplit, minbucket, maxdepth,
                   a, scale.y, useSearch, useOptim, useRpart, QCthreshold){
  #print(c(depth,id))
  if(depth > maxdepth){return(partynode(id=id))}
  
  y<-data[[response]]
  
  {
    #Select candidate variables
    #right here i need to insert the variable selection for both cat and cont
    #before it goes into partition or find split
    
    varSelected<-sort(sample(1:(ncol(data)-1),mtry,replace = F, prob = c(.2,.2,.2,.2,.2)))
    vars<-data[,!names(data) %in% response,drop=FALSE][,varSelected,drop=FALSE]
    colnames(vars)=varSelected #Have columns represent varid
    
    sp<-partition(vars=vars, y=y, subset=subset, search=search, method=method, split=split, nsplit=nsplit,
                  minsplit=minsplit, minbucket=minbucket, a=a, scale.y=scale.y,
                  useSearch=useSearch, useOptim=useOptim, useRpart=useRpart)
  }
  
  if(is.null(sp)){return(partynode(id=id))}
  
  # Split the data
  kidids<-kidids_split(sp, data=data)
  depth=depth+1
  
  kids<-vector(mode="list", length=max(kidids, na.rm=TRUE))
  for(kidid in 1:length(kids)){
    s<-subset
    s[kidids != kidid]<-FALSE
    # Node ID
    if(kidid > 1){myid<-max(nodeids(kids[[kidid-1]]))
    } else {myid<-id}
    # Start recursion on this daugther node
    kids[[kidid]]<-growtree(id=as.integer(myid+1), depth=depth, data=data, response=response, subset=s,
                            search=search, method=method, split=split, mtry=mtry, nsplit=nsplit,
                            minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth,
                            a=a, scale.y=scale.y, useSearch=useSearch, useOptim=useOptim, useRpart=useRpart,
                            QCthreshold=QCthreshold)
  }
  return(partynode(id=as.integer(id), split=sp, kids=kids,
                   info=list(stats=min(info_split(sp)$stats, na.rm=TRUE))))
}