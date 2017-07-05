grow<-function(formula, data=NULL, subset=NULL, search="greedy", method="class", split="gini",
               mtry=NULL, nsplit=NULL, minsplit=20, minbucket=round(minsplit/3), maxdepth=30,
               useRpart=FALSE){
  search <- match.arg(search,c("greedy","sss","crar"))
  method <- match.arg(method,c("anova","class"))
  split <- match.arg(split,c("information","ttest","MSE","entropy", "gini"))
  if(is.null(mtry)){mtry<-length(all.vars(formula[[3]]))}
  
  response<-all.vars(formula)[1]
  # Complete cases only.  Rearrange data so that response comes last
  data<-data[complete.cases(data),c(all.vars(formula)[-1],response)]
  #if(is.factor(data[[response]])){data[[response]]=as.numeric(data[[response]]==levels(data[[response]])[1])}
  
  if(is.null(subset)){subset<-rep(TRUE, nrow(data))}
  
  # Grow tree
  nodes<-growtree(id=1L, data=data, response=response, subset=subset, search=search, method=method, split=split,
                  mtry=mtry, nsplit=nsplit, minsplit=minsplit, minbucket=minbucket,
                  maxdepth=maxdepth, useRpart=useRpart)
  
  # Compute terminal node number for each observation
  fitted <- fitted_node(nodes, data=data)
  # Return rich constparty object
  ret <- party(nodes, data = data,
               fitted = data.frame("(fitted)" = fitted,
                                   "(response)" = data[[response]],
                                   "(weights)" = subset+0,
                                   check.names = FALSE),
               terms = terms(formula))
  as.constparty(ret)
}
