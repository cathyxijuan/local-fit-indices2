library(lavaan)
library(matrixcalc)
library(OpenMx)
library(expm)

pop.mod <- '     
xi1 =~ .5*x1 + .5*x2 + .4*x3 + .4*x4 
xi2 =~ .5*x5 + .5*x6 + .4*x7 + .4*x8
xi3 =~ .4*x9 + .4*x10 + .3*x11 + .3*x12
xi4 =~ .4*x13 + .4*x14 + .3*x15 + .3*x16
eta1 =~ .5*y1 + .5*y2 + .4*y3 + .4*y4 
eta2 =~ .4*y5 + .4*y6 + .3*y7 + .3*y8
eta3 =~ .5*y9 + .5*y10 + .3*y11 + .3*y12
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
x1 ~~ .75*x1
x2 ~~ .75*x2
x3 ~~ .84*x3
x4 ~~ .84*x4
x5 ~~ .75*x5
x6 ~~ .75*x6
x7 ~~ .84*x7
x8 ~~ .84*x8
x9 ~~ .84*x9
x10 ~~ .84*x10
x11 ~~ .91*x11
x12 ~~ .91*x12
x13 ~~ .84*x13
x14 ~~ .84*x14
x15 ~~ .91*x15
x16 ~~ .91*x16
y1 ~~ .75*y1
y2 ~~ .75*y2
y3 ~~ .84*y3
y4 ~~ .84*y4
y5 ~~ .84*y5
y6 ~~ .84*y6
y7 ~~ .91*y7
y8 ~~ .91*y8
y9 ~~ .75*y9
y10 ~~ .75*y10
y11 ~~ .91*y11
y12 ~~ .91*y12
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

path.mod0 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4
eta2 ~ xi1 + xi2 + eta1   
eta3 ~ eta2 + eta1
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"


path.mod1 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4 
eta2 ~ xi1 + xi2 + eta1   
eta3 ~ eta2              #delete eta1 here
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"

path.mod2 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4
eta2 ~ xi2 + eta1    #delete xi1 here
eta3 ~ eta2 + eta1
xi1 ~~ xi2 + xi3 + xi4
xi2 ~~ xi3  + xi4
xi3  ~~ xi4
"



path.mod3 <- "
eta1 ~ xi1 + xi2 + xi3 + xi4
eta2 ~ xi1 + eta1    #delete xi2 here
eta3 ~ eta2          #delete eta1 here
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


path.mod.list <- list(path.mod0, path.mod1, path.mod2, path.mod3, path.mod4, path.mod5)
