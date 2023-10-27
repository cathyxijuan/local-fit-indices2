source("models original.R")
source("functions.R")


n=400
#set.seed(777)
simuData<- simulateData(pop.mod, sample.nobs=n)
head(simuData)

#step 1
fit1 <- sem(sat.struct, data = simuData, estimator="ML")
lavInspect(fit1, "fit")
lavInspect(fit1, "converged")

#step 2 #I think stages 1 and 2 should really grouped into stage 1.
sat.struct.cov <- lavInspect(fit1, "cov.lv")[c(5:7, 1:4),c(5:7, 1:4) ] 
sat.struct.cov


fit2 <- sem(path.mod2, sample.cov = sat.struct.cov, sample.nobs = n, likelihood = "wishart")



struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                 "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                 "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                 "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                 "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) #the order of these names match the one for fit2 
fit1@Options$h1.information = "structured" 
W1.unstr.invert <- lavInspect(fit1, "inverted.information.observed")[struct.path, struct.path]


V1.unstr <- lavInspect(fit1, "information.first.order")[struct.path, struct.path]
#Note: W1.unstr.invert and V1.unstr should be the similar under normal data

#computing Gamma
Gamma.tri <- W1.unstr.invert%*%V1.unstr%*%W1.unstr.invert
Gamma <- W1.unstr.invert


c.adj.val(fit2, Gamma.tri, structured=T, expected=T)
c.adj.val(fit2, Gamma.tri, structured=T, expected=F)
c.adj.val(fit2, Gamma.tri, structured=F, expected=T)
c.adj.val(fit2, Gamma.tri, structured=F, expected=F)

rmsea.default <- lavInspect(fit2, "fit")["rmsea"]
rmsea.adj.unstr <- rmsea.adj(fit2, Gamma, structured = F)
rmsea.adj.str.exp <- rmsea.adj(fit2, Gamma, structured = T)
rmsea.adj.str.obs <- rmsea.adj(fit2, Gamma, structured = T, expected=F)
rmsea.adj.unstr.tri <-rmsea.adj(fit2, Gamma.tri, structured = F)
rmsea.adj.str.exp.tri <-rmsea.adj(fit2, Gamma.tri, structured = T)
rmsea.adj.str.obs.tri <-rmsea.adj(fit2, Gamma.tri, structured = T, expected=F)
rmsea.all <- c(rmsea.default,rmsea.adj.unstr, rmsea.adj.str.exp,rmsea.adj.str.obs,
               rmsea.adj.unstr.tri, rmsea.adj.str.exp.tri,rmsea.adj.str.obs.tri)


cfi.default <- lavInspect(fit2, "fit")["cfi"]
cfi.adj.unstr <- cfi.adj(fit2, Gamma, structured = F)
cfi.adj.str.exp <- cfi.adj(fit2, Gamma, structured = T)
cfi.adj.str.obs <- cfi.adj(fit2, Gamma, structured = F, expected = F)
cfi.adj.unstr.tri <- cfi.adj(fit2, Gamma.tri, structured = F)
cfi.adj.str.exp.tri <-cfi.adj(fit2, Gamma.tri, structured = T)
cfi.adj.str.obs.tri <-cfi.adj(fit2, Gamma.tri, structured = T, expected = F)
cfi.all <- c(cfi.default, cfi.adj.unstr,cfi.adj.str.exp,cfi.adj.str.obs,
             cfi.adj.unstr.tri,cfi.adj.str.exp.tri, cfi.adj.str.obs.tri)




srmr.unadj <- lavInspect(fit2, "fit")["srmr"]
srmr.default.adj.str <- srmr.adj.old(fit2, structured=T)
srmr.default.adj.unstr <- srmr.adj.old(fit2, structured=F)
srmr.adj.unstr <-srmr.adj(fit2, Gamma, structured = F)
srmr.adj.str.exp <- srmr.adj(fit2, Gamma, structured = T)
srmr.adj.str.obs <- srmr.adj(fit2, Gamma, structured = T, expected=F)
srmr.adj.unstr.tri <-srmr.adj(fit2, Gamma.tri, structured = F)
srmr.adj.str.exp.tri <-srmr.adj(fit2, Gamma.tri, structured = T)
srmr.adj.str.obs.tri <-srmr.adj(fit2, Gamma.tri, structured = T, expected = F)
srmr.all <- c(srmr.unadj, 
              srmr.default.adj.unstr,  srmr.default.adj.str, 
              srmr.adj.unstr, srmr.adj.str.exp, srmr.adj.str.obs,
              srmr.adj.unstr.tri, srmr.adj.str.exp.tri,srmr.adj.str.obs.tri)


 


fit.all <- c(rmsea.all, cfi.all, srmr.all)
fit.all
names(fit.all) <- c("rmsea.default",
                          "rmsea.adj.unstr", 
                          "rmsea.adj.str.exp",
                          "rmsea.adj.str.obs",
                          "rmsea.adj.unstr.tri", 
                          "rmsea.adj.str.exp.tri",
                          "rmsea.adj.str.obs.tri",
                          "cfi.default", 
                          "cfi.adj.unstr",
                          "cfi.adj.str.exp",
                          "cfi.adj.str.obs",
                          "cfi.adj.unstr.tri",
                          "cfi.adj.str.exp.tri", 
                          "cfi.adj.str.obs.tri",
                          "srmr.unadj",
                          "srmr.default.adj.unstr", 
                          "srmr.default.adj.str", 
                          "srmr.adj.unstr", 
                          "srmr.adj.str.exp", 
                          "srmr.adj.str.obs",
                          "srmr.adj.unstr.tri", 
                          "srmr.adj.str.exp.tri",
                          "srmr.adj.str.obs.tri")
length(fit.all)
fit.all
