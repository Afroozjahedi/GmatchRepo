smd = function(Data,Var){
        ASD = Data$group == 1 
  pars <- as.list(match.call())
  #formula <-group ~ RMSD.PRE.censoring + Age + WASI.NVIQ +Gender + Handedness 
    if (is.factor(Data[ASD,as.character(pars$Var)])){
    mASD = mean(as.numeric(Data[ASD,as.character(pars$Var)]),na.rm = TRUE)
    mTD = mean(as.numeric(Data[!ASD,as.character(pars$Var)]),na.rm = TRUE)
    varASD = var(as.numeric(Data[ASD,as.character(pars$Var)]),na.rm = TRUE)
    varTD = var(as.numeric(Data[!ASD,as.character(pars$Var)]),na.rm = TRUE)
    stdMeanDiff = (mASD - mTD)/(sqrt((varASD + varTD)/2))
    #cat("SMD",pars$Data,pars$Var,"=", stdMeanDiff,"\n")
    return(stdMeanDiff)
  }else {
  mASD = mean(Data[ASD,as.character(pars$Var)],na.rm = TRUE)
  mTD = mean(Data[!ASD,as.character(pars$Var)],na.rm = TRUE)
  varASD = var(Data[ASD,as.character(pars$Var)],na.rm = TRUE)
  varTD = var(Data[!ASD,as.character(pars$Var)],na.rm = TRUE)
  stdMeanDiff = (mASD - mTD)/(sqrt((varASD + varTD)/2))
  #cat("SMD",pars$Data,pars$Var,"=", stdMeanDiff,"\n")}
  return(stdMeanDiff)}
}
#debug(smd)
 #smd(nonOutlierLowADOS,RMSD.PRE.censoring)
#smdGender = smd(nonOutlierLowADOS,"Gender")
 
 #b <-  function(data,name) {
   
   ## match.call return a call containing the specified arguments 
   ## and the function name also 
   ## I convert it to a list , from which I remove the first element(-1)
   ## which is the function name
   
  # pars <- as.list(match.call()[-1])
  # data[,as.character(pars$name)]
   #cat(pars$name,"smd is",mean( data[,as.character(pars$name)]))
   #sd(data[,as.character(pars$name)])
   
# }

 