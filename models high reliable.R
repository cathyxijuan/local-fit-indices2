library(lavaan)
library(matrixcalc)
library(OpenMx)
library(expm)

pop.mod <- '     
xi1 =~ .9*x1 + .9*x2 + .8*x3 + .8*x4 
xi2 =~ .9*x5 + .9*x6 + .8*x7 + .8*x8
xi3 =~ .8*x9 + .8*x10 + .7*x11 + .7*x12
xi4 =~ .8*x13 + .8*x14 + .7*x15 + .7*x16
eta1 =~ .9*y1 + .9*y2 + .8*y3 + .8*y4 
eta2 =~ .8*y5 + .8*y6 + .7*y7 + .7*y8
eta3 =~ .9*y9 + .9*y10 + .7*y11 + .7*y12
eta1 ~ 0.3*xi1 + 0.4*xi2 + 0.2*xi3 + 0.3*xi4
eta2 ~ 0.3*xi1 + 0.4*xi2  + 0.3*eta1 
eta3 ~ 0.5*eta2 + 0.2*eta1
xi1 ~~ 0.3*xi2 + 0.3*xi3 + 0.3*xi4
xi2 ~~ 0.3*xi3 + 0.3*xi4 
xi3 ~~ 0.3*xi4 
xi1 ~~ 1*xi1
xi2 ~~ 1*xi2
xi3 ~~ 1*xi3
xi4 ~~ 1*xi4
x1 ~~ .19*x1
x2 ~~ .19*x2
x3 ~~ .36*x3
x4 ~~ .36*x4
x5 ~~ .19*x5
x6 ~~ .19*x6
x7 ~~ .36*x7
x8 ~~ .36*x8
x9 ~~ .36*x9
x10 ~~ .36*x10
x11 ~~ .51*x11
x12 ~~ .51*x12
x13 ~~ .36*x13
x14 ~~ .36*x14
x15 ~~ .51*x15
x16 ~~ .51*x16
y1 ~~ .19*y1
y2 ~~ .19*y2
y3 ~~ .36*y3
y4 ~~ .36*y4
y5 ~~ .36*y5
y6 ~~ .36*y6
y7 ~~ .51*y7
y8 ~~ .51*y8
y9 ~~ .19*y9
y10 ~~ .19*y10
y11 ~~ .51*y11
y12 ~~ .51*y12
eta1 ~~ .302*eta1 
eta2 ~~ 0.3318*eta2
eta3 ~~ 0.5646*eta3
'



fit.mod <- "
xi1 =~ x1 + x2 + x3 + x4 
xi2 =~ x5 + x6 + x7 + x8
xi3 =~ x9 + x10 + x11 + x12
xi4 =~ x13 + x14 + x15 + x16
eta1 =~ y1 + y2 + y3 + y4 
eta2 =~ y5 + y6 + y7 + y8
eta3 =~ y9 + y10 + y11 + y12
eta1 ~ xi1 + xi2 + xi3
eta2 ~ xi1 + xi2 + eta1   
eta3 ~ eta2
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"


sat.struct <- " 
xi1 =~ x1 + x2 + x3 + x4 
xi2 =~ x5 + x6 + x7 + x8
xi3 =~ x9 + x10 + x11 + x12
xi4 =~ x13 + x14 + x15 + x16
eta1 =~ y1 + y2 + y3 + y4 
eta2 =~ y5 + y6 + y7 + y8
eta3 =~ y9 + y10 + y11 + y12 
eta1 ~~ eta2 + eta3 + xi1 +xi2 + xi3 +xi4
eta2 ~~ eta3 + xi1 +xi2 + xi3 +xi4
eta3 ~~ xi1 +xi2 + xi3 +xi4
xi1 ~~ xi2 + xi3 +xi4
xi2 ~~  xi3 +xi4
xi3 ~~ xi4
"

path.mod1 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4
eta2 ~ xi1 + xi2 + eta1   
eta3 ~ eta2 + eta1
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"


path.mod2 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4 
eta2 ~ xi1 + xi2 + eta1   
eta3 ~ eta2              #delete eta1 here
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"

path.mod3 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4
eta2 ~ xi2 + eta1    #delete xi1 here
eta3 ~ eta2 + eta1
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"


path.mod4 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4
eta2 ~ xi2  + eta1     #delete xi1 here
eta3 ~ eta2            #delete eta1 here
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"

path.mod5 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4
eta2 ~ xi1 + eta1    #delete xi2 here
eta3 ~ eta2          #delete eta1 here
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"



path.mod6 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4
eta2 ~ eta1    #delete xi1 and xi2 here
eta3 ~ eta2    #delete eta1 here
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"




path.mod.HW <- " #Hao's original model (problem: xi4 is not used in the regression equations..)
eta1 ~ xi1 + xi2 + xi3 
eta2 ~ xi1 + xi2 + eta1   
eta3 ~ eta2
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"


path.mod.list <- list(path.mod1, path.mod2, path.mod3, path.mod4, path.mod5, path.mod6)


### calculate true population SRMR ###
# Data<-simulateData(pop.mod,sample.nobs = 500)
# fit1 <- sem(pop.mod, data=Data)
# Sigma<-lavInspect(fit1,"cov.lv")
# 
# lavInspect(fit1,"est")$psi
# mod.list <- list(path.mod0, path.mod1, path.mod2, path.mod3,path.mod4,path.mod5)
# 
# pop.fit.mod.orig <- matrix(NA, nrow=3, ncol=length(mod.list))
# rownames(pop.fit.mod.orig) <- c("pop.rmsea", "pop.cfi", "pop.srmr")
# colnames(pop.fit.mod.orig) <- paste("path.mod", 0:(length(mod.list)-1), sep="")
# 
# for(i in 1:length(mod.list)){
#   fit2 <- sem(mod.list[[i]], sample.cov = Sigma, sample.nobs = 500, likelihood = "wishart")
#   pop.fit.mod.orig[,i] <- lavInspect(fit2, "fit")[c("rmsea", "cfi", "srmr")]
# }
# round(pop.fit.mod.orig, 5)
# 
# 
#round(pop.fit.mod.orig,4)



