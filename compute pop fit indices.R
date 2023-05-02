source("models original.R")


### calculate true population SRMR ###

Data<-simulateData(pop.mod,sample.nobs = 10000)
fit1 <- sem(pop.mod, data=Data)
Sigma<-lavInspect(fit1,"cov.lv")

lavInspect(fit1,"est")$psi
mod.list <- list(path.mod0, path.mod1, path.mod2, path.mod3,path.mod4,path.mod5)

pop.fit.mod.orig <- matrix(NA, nrow=3, ncol=length(mod.list))
rownames(pop.fit.mod.orig) <- c("pop.rmsea", "pop.cfi", "pop.srmr")
colnames(pop.fit.mod.orig) <- paste("path.mod", 0:(length(mod.list)-1), sep="")

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


#          path.mod0  path.mod1  path.mod2 path.mod3 path.mod4  path.mod5
#pop.rmsea 1.720638e-08 0.06851110 0.15423663 0.1874318 0.1578682 0.20660795
#pop.cfi   1.000000e+00 0.99004460 0.94954412 0.9148440 0.9395887 0.88359401
#pop.srmr  2.085630e-08 0.02406553 0.04239781 0.0600216 0.0534605 0.07077444
