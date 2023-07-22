source("functions.R")
load("pop.indices.RData")
load("srmr.ci.check.val.mod.orig.n150.RData")
load("srmr.ci.check.val.mod.orig.n200.RData")
load("srmr.ci.check.val.mod.orig.n300.RData")
load("srmr.ci.check.val.mod.orig.n500.RData")
load("srmr.ci.check.val.mod.orig.n800.RData")
load("srmr.ci.check.val.mod.orig.n1000.RData")
load("srmr.ci.check.val.mod.high.n150.RData")
load("srmr.ci.check.val.mod.high.n200.RData")
load("srmr.ci.check.val.mod.high.n300.RData")
load("srmr.ci.check.val.mod.high.n500.RData")
load("srmr.ci.check.val.mod.high.n800.RData")
load("srmr.ci.check.val.mod.high.n1000.RData")


#looking only at the k_s values
check.list <- srmr.ci.check.val.mod.orig.n1000
check.list <- check.list[!(sapply(check.list, is.null))]
test <- simplify2array(check.list, higher=T)
ks <- round(cbind(test[1,1,], test[1,2,], test[1,3,], 
            test[1,4,], test[1,5,], test[1,6,]), 3)
colnames(ks) <- paste("path.mod", 0:5, sep="")

#write.csv(ks, "ksn500.csv", row.names=F)

round(var(ks),3)
hist(ks[,1])
table(round(ks[,1],1))
0.01/0.6
0.01/0.8





srmr.check.val.orig.mod <- list("n150"=list.mean(srmr.ci.check.val.mod.orig.n150),
                                "n200"=list.mean(srmr.ci.check.val.mod.orig.n200),
                                "n300"=list.mean(srmr.ci.check.val.mod.orig.n300),
                                "n500"=list.mean(srmr.ci.check.val.mod.orig.n500),
                                "n800"=list.mean(srmr.ci.check.val.mod.orig.n800),
                                "n1000"=list.mean(srmr.ci.check.val.mod.orig.n1000))

save(srmr.check.val.orig.mod, 
     file="srmr.check.val.orig.mod.RData")



srmr.check.val.high.mod <- list("n150"=list.mean(srmr.ci.check.val.mod.high.n150),
                                "n200"=list.mean(srmr.ci.check.val.mod.high.n200),
                                "n300"=list.mean(srmr.ci.check.val.mod.high.n300),
                                "n500"=list.mean(srmr.ci.check.val.mod.high.n500),
                                "n800"=list.mean(srmr.ci.check.val.mod.high.n800),
                                "n1000"=list.mean(srmr.ci.check.val.mod.high.n1000))

save(srmr.check.val.high.mod, 
     file="srmr.check.val.high.mod.RData")
