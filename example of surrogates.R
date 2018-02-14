library(partykit)
# This is a node with specific split rule and coresponding kid nodes. No data in it yet.
pn <- partynode(id = 1L, split = partysplit(varid=2L, breaks = 3),
kids = lapply(2:3, partynode))

irisNAVar= as.data.frame(lapply(iris[-5],function(x) x[sample(c(TRUE,NA),prob=c(0.85,0.15), size=length(x), replace=TRUE)]))  
irisNA=data.frame(Species=iris[,5],irisNAVar)
# tree <- as.constparty(party(pn, data = iris,
#                             fitted = data.frame("(fitted)" = fitted_node(pn, data = iris),
#                                                 "(response)" = iris$Species,check.names = FALSE),
#                             terms = terms(Species ~ ., data = iris)))

set.seed(1)
treeNA <- as.constparty(party(pn, data = irisNA,
                            fitted = data.frame("(fitted)" = fitted_node(pn, data = irisNA),
                                                "(response)" = irisNA$Species,check.names = FALSE),
                            terms = terms(Species ~ ., data = irisNA)))
#plot(tree)
plot(treeNA)
set.seed(2)
treeNA <- as.constparty(party(pn, data = irisNA,
                              fitted = data.frame("(fitted)" = fitted_node(pn, data = irisNA),
                                                  "(response)" = irisNA$Species,check.names = FALSE),
                              terms = terms(Species ~ ., data = irisNA)))
plot(treeNA)
#Manually add a missing Sepal Width
#newdata <- data.frame(Sepal.Length=5, Sepal.Width=NA, Petal.Length=1.3, Petal.Width=0.2)


#Since Sepal Width is missing, randomly assign a node.  Suppress warning
set.seed(1)
a=suppressWarnings(predict(treeNA, newdata=irisNA, type="node"))
set.seed(4)
b=suppressWarnings(predict(treeNA, newdata=irisNA, type="node"))
