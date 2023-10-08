library(lavaan)
library(matrixcalc)
library(OpenMx)
library(expm)

pop.mod <- '     
xi1 =~ .7*x1 + .7*x2 + .6*x3 + .6*x4 
xi2 =~ .7*x5 + .7*x6 + .6*x7 + .6*x8
xi3 =~ .6*x9 + .6*x10 + .5*x11 + .5*x12
xi4 =~ .6*x13 + .6*x14 + .5*x15 + .5*x16
eta1 =~ .7*y1 + .7*y2 + .6*y3 + .6*y4 
eta2 =~ .6*y5 + .6*y6 + .5*y7 + .5*y8
eta3 =~ .7*y9 + .7*y10 + .5*y11 + .5*y12
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
x1 ~~ .51*x1
x2 ~~ .51*x2
x3 ~~ .64*x3
x4 ~~ .64*x4
x5 ~~ .51*x5
x6 ~~ .51*x6
x7 ~~ .64*x7
x8 ~~ .64*x8
x9 ~~ .64*x9
x10 ~~ .64*x10
x11 ~~ .75*x11
x12 ~~ .75*x12
x13 ~~ .64*x13
x14 ~~ .64*x14
x15 ~~ .75*x15
x16 ~~ .75*x16
y1 ~~ .51*y1
y2 ~~ .51*y2
y3 ~~ .64*y3
y4 ~~ .64*y4
y5 ~~ .64*y5
y6 ~~ .64*y6
y7 ~~ .75*y7
y8 ~~ .75*y8
y9 ~~ .51*y9
y10 ~~ .51*y10
y11 ~~ .75*y11
y12 ~~ .75*y12
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

