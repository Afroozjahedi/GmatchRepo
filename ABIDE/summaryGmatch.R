summaryGmatch <- function(DATA){
#=== Tabling pvalues and SMD for output ====
        GROUP = DATA$GROUP
        summary <- aggregate(.~GROUP,DATA,function(x) c(mean = mean(x), sd = sd(x),range = range(x)))
        table(DATA$GROUP,DATA$HANDEDNESS)
        pvalues <- t(sapply(DATA[c(-1,-5)], function(x) unlist(t.test(x~DATA$GROUP)["p.value"])))
        # pvalHandedness <-
        #         chisq.test(DATA$GROUP, DATA$HANDEDNESS)$p.value
        
#         apply(DATA[-1],2, smd)
#         #debug(smd)
#         #smd(nonOutlierLowADOS,RMSD.PRE.censoring)
#         #smdGender = smd(nonOutlierLowADOS,"Gender")
#         
#         #b <-  function(data,name) {
#         
#         ## match.call return a call containing the specified arguments 
#         ## and the function name also 
#         ## I convert it to a list , from which I remove the first element(-1)
#         ## which is the function name
#         
#         # pars <- as.list(match.call()[-1])
#         # data[,as.character(pars$name)]
#         #cat(pars$name,"smd is",mean( data[,as.character(pars$name)]))
#         #sd(data[,as.character(pars$name)])
#         
#         # }
#         
#         
# library(dplyr)  
#         DATA %>% 
#                 group_by(GROUP) %>%
#                 do(data.frame(val=smd(.)))
#         
# Calculate standardized mean difference
smdMotion <- smd(DATA, RMSD)
smdAge <- smd(DATA, "AGE")
#smdFIQ <- smd(DATA, "FIQ")
#smdHandedness <- smd(DATA, "Handedness")

#===Tabling smd for output ===
smdOpt <-
        matrix(c(smdMotion, smdAge),ncol = 2)
colnames(smdOpt) <-
        c("RMSD",
          "AGE")
rownames(smdOpt) <- "SMD"
smdOpt <- as.table(smdOpt)
#return(smdMotion )#= smdMotion)#,smdAge = smdAge, smdOpt = smdOpt))
return(list(smdMotion = smdMotion,smdAge = smdAge,smdFIQ = smdFIQ,
                       print(summary),print(smdOpt),print(table(DATA$GROUP))),print(pvalues))
#cat("pvalsOpt", "smdOpt","table(DATA$GROUP)","\n")
}
