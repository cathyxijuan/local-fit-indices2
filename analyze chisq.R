load("chisq.mod.orig.n150.RData")


colMeans(chisq.mod.orig.n150, na.rm=T)
colMeans(chisq.mod.orig.n150[,4:6]<=0.05, na.rm=T)

pval <- chisq.mod.orig.n150[complete.cases(chisq.mod.orig.n150),4:6]
colMeans(pval<0.05)