source("models original.R")
source("models high reliable.R")


### calculate true population fit indices ###

Data<-simulateData(pop.mod,sample.nobs = 100000)
fit1 <- sem(pop.mod, data=Data)
Sigma<-lavInspect(fit1,"cov.lv")

#lavInspect(fit1,"est") #it outputs all the parameters set in the pop.mod
mod.list <- list(path.mod1, path.mod2, path.mod3, path.mod4,path.mod5,path.mod6)

pop.fit.mod.orig <- matrix(NA, nrow=3, ncol=length(mod.list))
rownames(pop.fit.mod.orig) <- c("pop.rmsea", "pop.cfi", "pop.srmr")
colnames(pop.fit.mod.orig) <- paste("path.mod", 1:length(mod.list), sep="")

for(i in 1:length(mod.list)){
  fit2 <- sem(mod.list[[i]], sample.cov = Sigma, sample.nobs = 10000000, likelihood = "wishart")
  fmin <- lavInspect(fit2, "fit")["fmin"]*2
  df <- lavInspect(fit2, "fit")["df"]
  fminB <- lavInspect(fit2, "fit")["baseline.chisq"]/lavInspect(fit2, "fit")["ntotal"]
  rmsea <- sqrt(fmin/df)
  cfi <-1-fmin/fminB
  srmr <- lavInspect(fit2, "fit")[ "srmr"]
  pop.fit.mod.orig[,i] <- c(rmsea,cfi,srmr)
}
round(pop.fit.mod.orig, 5)

pop.rmsea <- pop.fit.mod.orig[1,]
pop.cfi <- pop.fit.mod.orig[2,]
pop.srmr <- pop.fit.mod.orig[3,]
pop.indices <- rbind(pop.rmsea, pop.rmsea, pop.rmsea, pop.rmsea, pop.rmsea,pop.rmsea,pop.rmsea,
                     pop.cfi, pop.cfi, pop.cfi, pop.cfi,pop.cfi,pop.cfi,pop.cfi,
                     pop.srmr, pop.srmr, pop.srmr, pop.srmr, pop.srmr, 
                     pop.srmr, pop.srmr,pop.srmr,pop.srmr)
pop.indices

save(pop.indices, file="pop.indices.RData")

#            path.mod1 path.mod2 path.mod3 path.mod4 path.mod5 path.mod6
#pop.rmsea         0   0.06851   0.15424   0.15787   0.18743   0.20661
#pop.cfi           1   0.99004   0.94954   0.93959   0.91484   0.88359
#pop.srmr          0   0.02407   0.04240   0.05346   0.06002   0.07077

