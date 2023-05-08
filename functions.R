####Usage: find a matrix of means based on a list of matrices
#Argument: lis: a list of matrices
list.mean <- function(lis ){
  if(sum(sapply(lis, is.null))==0){
    apply(simplify2array(lis), 1:2, mean, na.rm = T)
  } else{
    lis <- lis[!(sapply(lis, is.null))]
    apply(simplify2array(lis), 1:2, mean, na.rm = T)
    }
}


####Usage: find a matrix of sds based on a list of matrices
#Argument: lis: a list of matrices
list.sd <- function(lis ){
  if(sum(sapply(lis, is.null))==0){
    round(apply(simplify2array(lis), 1:2, sd, na.rm = T),4)
  } else{
    lis <- lis[!(sapply(lis, is.null))]
    round(apply(simplify2array(lis), 1:2, sd, na.rm = T),4)
  }
}




####Usage: root mean square error for our simulation results 
#Argument: simu.list: a list of matrices from our small sample simulation studies
########## pop.matrix: a matrix of population values from our previous paper. 


rmse <- function(simu.list, pop.matrix){
  if(sum(sapply(simu.list, is.null))==0){
    to.array <- simplify2array(simu.list)
    subtr.pop <- sweep(to.array, 1:2, pop.matrix)
    sqr.dif  <- subtr.pop ^2
    exp.val <- apply( sqr.dif, 1:2, mean, na.rm=T)
    sqrt.val  <- sqrt(exp.val)
    round(sqrt.val,4)
  } else{
    simu.list <- simu.list[!(sapply(simu.list, is.null))]
    to.array <- simplify2array(simu.list)
    subtr.pop <- sweep(to.array, 1:2, pop.matrix)
    sqr.dif  <- subtr.pop ^2
    exp.val <- apply( sqr.dif, 1:2, mean, na.rm=T)
    sqrt.val  <- sqrt(exp.val)
    round(sqrt.val ,4)
  }
}
  







##Purpose: compute the new adjusted SRMR
##Argument:
## 1) fit: the SEM fit from step 2
## 2) gamma: a gamma matrix computed using fit1 in step 1; can either be triple product or not
## 3) structured: logical; whether the weight matrix should be structured.
## 4) expected: logical; whether the weight matrix should be expected or observed; only differ if structured=T
srmr.adj <- function(fit, gamma, structured=T, expected=T ){
  if (structured==T){
    fit@Options$h1.information = "structured" 
    if(expected==T){
      W <- lavInspect(fit, "h1.information.expected")
    } else{
      W <- lavInspect(fit, "h1.information.observed")
    }
  } else {
    fit@Options$h1.information = "unstructured" 
    W <- lavInspect(fit, "wls.v")
  }
  
  n <- lavInspect(fit, "nobs")
  residual.cov <- lavInspect(fit, "resid")$cov
  S <- lavInspect(fit, "sampstat")$cov 
  
  samp.var <- diag(S)
  resid.cov.vech <- OpenMx::vech(residual.cov)
  p <-ncol(residual.cov)
  t <- p*(p+1)/2
  
  s1 <- OpenMx::vech(samp.var%*%t(samp.var))
  G <- diag(s1) #G matrix in Maydeu-Olivares 2017's paper; used later in the computation of unbiased SRMR
  G.inv <- diag(1/s1)
  T_s  <- t(resid.cov.vech)%*%G.inv%*%resid.cov.vech
  
  
  delta <- lavInspect(fit, "delta")
  # W.R<-chol(W)
  # H<-delta%*%qr.solve(W.R%*%delta, W.R) #has problems when W is not positive definite, which happens when W is structured and observed
  H <- delta%*%solve(t(delta)%*%W%*%delta)%*%t(delta)%*%W
  G_sqrt_inv <- diag(1/sqrt(s1)) #inverse of the square root of G 
  e_s <- G_sqrt_inv%*%resid.cov.vech #before equation 15
  
  Xi_u <- (diag(nrow(H)) - H)%*%gamma%*%t((diag(nrow(H)) - H))/n #before Equation 16
  Xi_s <- G_sqrt_inv%*%Xi_u%*%G_sqrt_inv
  k_s <- 1 -((sum(Xi_s^2) + 2*t(e_s)%*%Xi_s%*%e_s)/(4*T_s^2))
  
  tr.Xi<-sum(diag(Xi_s))
  SRMR.sq<-(T_s-tr.Xi)/t
  
  sqrt(max(0,SRMR.sq))/k_s
}








##Purpose: compute the old/default adjusted SRMR in step 2 of the two-stage procedure
##Argument:
## 1) fit: the SEM fit from step 2
## 3) structured: logical; whether the weight matrix should be structured. 
srmr.adj.old <- function(fit, structured){
  if (structured==T){
    fit@Options$h1.information = "structured" 
  } else {
    fit@Options$h1.information = "unstructured" 
  }
  lavResiduals(fit)$summary$cov[5]
}



##Purpose: compute the new adjusted SRMR confidence interval (CI)
##Argument:
## 1) fit: the SEM fit from step 2
## 2) gamma: a gamma matrix computed using fit1 in step 1; can either be triple product or not
## 3) structured: logical; whether the weight matrix should be structured.
## 4) expected: logical; whether the weight matrix should be expected or observed; only differ if structured=T
srmr.adj.ci <- function(fit, gamma, structured=T, expected=T){
  if (structured==T){
    fit@Options$h1.information = "structured" 
    if(expected==T){
      W <- lavInspect(fit, "h1.information.expected")
    } else{
      W <- lavInspect(fit, "h1.information.observed")
    }
  } else {
    fit@Options$h1.information = "unstructured" 
    W <- lavInspect(fit, "wls.v")
  }
  
  n <- lavInspect(fit, "nobs")
  residual.cov <- lavInspect(fit, "resid")$cov
  S <- lavInspect(fit, "sampstat")$cov 
  
  samp.var <- diag(S)
  resid.cov.vech <- OpenMx::vech(residual.cov)
  p <-ncol(residual.cov)
  t <- p*(p+1)/2
  
  s1 <- OpenMx::vech(samp.var%*%t(samp.var))
  G <- diag(s1) #G matrix in Maydeu-Olivares 2017's paper; used later in the computation of unbiased SRMR
  G.inv <- diag(1/s1)
  T_s  <- t(resid.cov.vech)%*%G.inv%*%resid.cov.vech
  
  
  delta <- lavInspect(fit, "delta")
  # W.R<-chol(W)
  # H<-delta%*%qr.solve(W.R%*%delta, W.R) #has problems when W is not positive definite, which happens when W is structured and observed
  H <- delta%*%solve(t(delta)%*%W%*%delta)%*%t(delta)%*%W
  G_sqrt_inv <- diag(1/sqrt(s1)) #inverse of the square root of G 
  e_s <- G_sqrt_inv%*%resid.cov.vech #before equation 15
  
  Xi_u <- (diag(nrow(H)) - H)%*%gamma%*%t((diag(nrow(H)) - H))/n #before Equation 16
  Xi_s <- G_sqrt_inv%*%Xi_u%*%G_sqrt_inv
  tr.Xissqr<-sum(Xi_s^2)
  k_s <- 1 -((tr.Xissqr + 2*t(e_s)%*%Xi_s%*%e_s)/(4*T_s^2))
  
  tr.Xi<-sum(diag(Xi_s))
  SRMR.sq<-(T_s-tr.Xi)/t
  srmr_unbias.new <-  sqrt(max(0,SRMR.sq))/k_s
  
  se.new <- sqrt(k_s^(-2)*(tr.Xissqr+2*t(e_s)%*%Xi_s%*%e_s)/(2*t*T_s))
  alpha=0.1
  c<-1-alpha/2
  zc<-qnorm(c) 
  srmr_unbias.new.ci.lower <- max(0, srmr_unbias.new-zc*se.new )
  srmr_unbias.new.ci.upper <- srmr_unbias.new+zc*se.new
  c(srmr_unbias.new.ci.lower, srmr_unbias.new.ci.upper)
}






##Purpose: compute the old/default adjusted SRMR confidence interval (CI) in step 2 of the two-stage procedure
##Argument:
## 1) fit: the SEM fit from step 2
## 3) structured: logical; whether the weight matrix should be structured. 
srmr.adj.old.ci <- function(fit, structured=T){
  if (structured==T){
    fit@Options$h1.information = "structured" 
  } else {
    fit@Options$h1.information = "unstructured" 
  }
  c(lavResiduals(fit)$summary$cov[7],lavResiduals(fit)$summary$cov[8])
}





##Purpose: compute the adjust constant for the hypothesized model for computing RMSEA and CFI
##Argument:
## 1) fit: the SEM fit from step 2
## 2) gamma: a gamma matrix computed using fit1 in step 1; can either be triple product or not
## 3) structured: logical; whether the weight matrix should be structured. 
## 4) expected: logical; whether the expected weight matrix should be used. It only varies under structured=T
c.adj.val <- function(fit, gamma, structured, expected=T){
  if (structured==T){
    fit@Options$h1.information = "structured" 
    if(expected==T){
      W2 <- lavInspect(fit, "h1.information.expected")
    } else{
      W2 <- lavInspect(fit, "h1.information.observed")
    }
  } else {
    fit@Options$h1.information = "unstructured" 
    W2 <- lavInspect(fit, "wls.v")
  }
  
  deltabreve <- lavInspect(fit, "delta") #jacobian matrix
  Uc <- 
    W2-W2%*%
    deltabreve%*%
    solve(t(deltabreve)%*%
            W2%*%deltabreve)%*%
    t(deltabreve)%*%W2
  
  lav_matrix_trace(Uc%*%gamma) 
}



##Purpose: compute the new adjusted RMSEA 
##Argument:
## 1) fit: the SEM fit from step 2
## 2) gamma: a gamma matrix computed using fit1 in step 1; can either be triple product or not
## 3) structured: logical; whether the weight matrix should be structured. 
## 4) expected: logical; whether the expected weight matrix should be used. 
rmsea.adj <- function(fit, gamma, structured=T, expected=T){
  n <- lavInspect(fit, "nobs")
  c.adj <- c.adj.val(fit, gamma, structured, expected)
  Fmin <- lavInspect(fit, "fit")["fmin"]*2
  df <- lavInspect(fit, "fit")["df"]
  sqrt(max((Fmin- c.adj/n)/df, 0))
}






##Purpose: compute adjusted value for the baseline model for computing RMSEA and CFI
##Argument:
## 1) fit: the SEM fit from step 2
## 2) gamma: a gamma matrix computed using fit1 in step 1; can either be triple product or not
## 3) structured: logical; whether the weight matrix should be structured. 
## 4) expected: logical; whether the expected weight matrix should be used. 
cB.adj.val <- function(fit, gamma, structured=T, expected=T){
  fit2B <-  lavaan:::lav_object_independence(fit, se=T) 
  if (structured==T){
    fit@Options$h1.information = "structured" 
    if (expected==T){
      W2B <-lavInspect(fit2B, "h1.information.expected")
    } else {
      W2B <-lavInspect(fit2B, "h1.information.observed")
    }
  } else {
    fit@Options$h1.information = "unstructured" 
    W2B <-lavInspect(fit2B, "wls.v")
  }
  
  #should we vary it under the structured case, which is the same as "wls.v"
  deltabreveB<- lavInspect(fit2B, "delta") 
  
  Uc.B <- 
    W2B-W2B%*%
    deltabreveB%*%
    solve(t(deltabreveB)%*%
            W2B%*%deltabreveB)%*%
    t(deltabreveB)%*%W2B
  
  lav_matrix_trace(Uc.B%*%gamma) 
}




##Purpose: compute the new adjusted CFI. 
##Argument:
## 1) fit: the SEM fit from step 2
## 2) gamma: a gamma matrix computed using fit1 in step 1; can either be triple product or not
## 3) structured: logical; whether the weight matrix should be structured. 
## 4) expected: logical; whether the expected weight matrix should be used. 
cfi.adj <- function(fit, gamma, structured=T, expected=T){
  if (structured==T){
    fit@Options$h1.information = "structured" 
  } else {
    fit@Options$h1.information = "unstructured" 
  }
  n <- lavInspect(fit, "nobs")
  fit2B <-  lavaan:::lav_object_independence(fit, se=T) 
  c.adj <- c.adj.val(fit, gamma, structured, expected)
  cB.adj <- cB.adj.val(fit, gamma, structured, expected)
  
  dfB <- lavInspect(fit2B, "fit")
  FminB <- lavInspect(fit2B, "fit")["fmin"]*2
  Fmin <- lavInspect(fit, "fit")["fmin"]*2
  df <- lavInspect(fit, "fit")["df"]
  
  if(max(FminB-cB.adj/n, Fmin-c.adj/n, 0) ==0 ){
    1
  } else{
    1- max(Fmin-c.adj/n, 0) / max(FminB-cB.adj/n, Fmin-c.adj/n)}  
}











##Purpose: generate a list of matrices with fit indices for a given condition. 
### Each simulated dataset has a matrix of fit indices with rows being different kinds of fit and columns being different path models. 
### These matrices are combined as a list across repetitions. 
####Argument:
#pop.model: the population model that is used for generating data. 
#sat.model: the model used in step 1 of the two-stage procedure; it is saturated in the latent variables. 
#path.model.list: a list of path model used in step 2 of the two-stage procedure
#sample.size: sample size in the condition
#rep.num: number of repetitions for the stimulation study. 
simu.fit <- function(pop.model, sat.model, path.model.list, sample.size, rep.num){
  num.fit.indices <- 23
  num.path.mod <- length(path.model.list)
  
  fit.list <- vector(mode="list", length=rep.num)
  struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                   "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                   "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                   "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                   "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) 
  
  for(j in 1:rep.num){
    fit.matrix <- matrix(nrow=num.fit.indices, ncol=num.path.mod)
    colnames(fit.matrix) <-paste("path.mod", 0:(num.path.mod-1), sep="")
    rownames(fit.matrix) <- c("rmsea.default",
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
                              "srmr.default.adj.str", 
                              "srmr.default.adj.unstr", 
                              "srmr.adj.unstr", 
                              "srmr.adj.str.exp", 
                              "srmr.adj.str.obs",
                              "srmr.adj.unstr.tri", 
                              "srmr.adj.str.exp.tri",
                              "srmr.adj.str.obs.tri")
    
    simuData<- simulateData(pop.model, sample.nobs=sample.size)
    fit1 <- sem(sat.model, data = simuData, estimator="ML", likelihood="wishart")
    if (lavInspect(fit1, "converged")==F){ 
      print(j)
    } else {
      sat.struct.cov <- lavInspect(fit1, "cov.lv")[c(5:7, 1:4),c(5:7, 1:4) ]
      if (is.positive.definite(sat.struct.cov)==F){
        print(j)
      } else{
        #compute gamma
        fit1@Options$h1.information = "structured" 
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.observed")[struct.path, struct.path]
        V1.unstr <- lavInspect(fit1, "information.first.order")[struct.path, struct.path]
        Gamma.tri <- W1.unstr.invert%*%V1.unstr%*%W1.unstr.invert
        Gamma <- W1.unstr.invert
        
        
        
        for(i in 1:num.path.mod){
          fit2 <- sem(path.model.list[[i]], sample.cov = sat.struct.cov, sample.nobs = sample.size, likelihood = "wishart")
          
          #compute fit indices
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
          srmr.all <- c(srmr.unadj, srmr.default.adj.str, srmr.default.adj.unstr, 
                        srmr.adj.unstr, srmr.adj.str.exp, srmr.adj.str.obs,
                        srmr.adj.unstr.tri, srmr.adj.str.exp.tri,srmr.adj.str.obs.tri)
          
          fit.all <- c(rmsea.all, cfi.all, srmr.all)
          fit.matrix[,i] <- fit.all 
          
        }
        
        fit.matrix <- round(fit.matrix, 8)
        fit.list[[j]] <- fit.matrix
        print(j)
        
      }
    }
  }
  fit.list
  
}





##Purpose: count the nonconvergence number in a simulation 
####Argument:
#pop.model: the population model that is used for generating data. 
#sat.model: the model used in step 1 of the two-stage procedure; it is saturated in the latent variables. 
#sample.size: sample size in the condition
#rep.num: number of repetitions for the stimulation study. 
##VALUE: numeric number for nonconvergence count. 
nonconverge.count <- function(pop.model, sat.model, sample.size, rep.num){
  nonconverge <- 0
  for(j in 1:rep.num){
    simuData<- simulateData(pop.model, sample.nobs=sample.size)
    fit1 <- sem(sat.model, data = simuData, estimator="ML")
    if (lavInspect(fit1, "converged")==F){ 
      nonconverge <- nonconverge + 1
      print(j)
    } else { 
      print(j)
    }
  }
  nonconverge
}



##Purpose: count the number of nulls in a list
##Argument: 
#matrix.list: a list of matrix that may contain null values
count.null <- function(matrix.list){
  sum(sapply(matrix.list, is.null))
}


##Purpose: remove the null matrices in a list of matrices
##Argument: 
#matrix.list: a list of matrix that may contain null values
remove.null <- function(matrix.list){
  if(sum(sapply(matrix.list, is.null))==0){
    matrix.list
  } else{
    matrix.list[!(sapply(matrix.list, is.null))]
  }
}




##Purpose: generate a list of matrices with components that compute indices for a given condition. 
### Each simulated dataset has a matrix of fit indices with rows being different kinds of fit and columns being different path models. 
### These matrices are combined as a list across repetitions. 
####Argument:
#pop.model: the population model that is used for generating data. 
#sat.model: the model used in step 1 of the two-stage procedure; it is saturated in the latent variables. 
#path.model.list: a list of path model used in step 2 of the two-stage procedure
#sample.size: sample size in the condition
#rep.num: number of repetitions for the stimulation study. 

simu.component <- function(pop.model, sat.model, path.model.list, sample.size, rep.num){
  num.fit.indices <- 17 #need to change
  num.path.mod <- length(path.model.list)
  
  component.list <- vector(mode="list", length=rep.num)
  struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                   "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                   "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                   "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                   "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) 
  
  for(j in 1:rep.num){
    component.matrix <- matrix(nrow=num.fit.indices, ncol=num.path.mod)
    colnames(component.matrix) <-paste("path.mod", 0:(num.path.mod-1), sep="")
    rownames(component.matrix) <- c("rmsea.default","rmsea.adj.unstr", "rmsea.adj.str",
                                    "rmsea.adj.unstr.tri",   "rmsea.adj.str.tri", 
                                    "cfi.default", "cfi.adj.unstr", "cfi.adj.str",
                                    "cfi.adj.unstr.tri","cfi.adj.str.tri", 
                                    "srmr.unadj", "srmr.default.adj.str", "srmr.default.adj.unstr", 
                                    "srmr.adj.unstr", "srmr.adj.str", "srmr.adj.unstr.tri", 
                                    "srmr.adj.str.tri")
    
    simuData<- simulateData(pop.model, sample.nobs=sample.size)
    fit1 <- sem(sat.model, data = simuData, estimator="ML", likelihood="wishart")
    if (lavInspect(fit1, "converged")==F){ 
      print(j)
    } else {
      sat.struct.cov <- lavInspect(fit1, "cov.lv")[c(5:7, 1:4),c(5:7, 1:4) ]
      if (is.positive.definite(sat.struct.cov)==F){
        print(j)
      } else{
        #compute gamma
        fit1@Options$h1.information = "unstructured" 
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.observed")[struct.path, struct.path]
        V1.unstr <- lavInspect(fit1, "information.first.order")[struct.path, struct.path]
        Gamma.tri <- W1.unstr.invert%*%V1.unstr%*%W1.unstr.invert
        Gamma <- W1.unstr.invert
        
        
        
        for(i in 1:num.path.mod){
          fit2 <- sem(path.model.list[[i]], sample.cov = sat.struct.cov, sample.nobs = sample.size, likelihood = "wishart")
          
          #compute fit indices
          fmin <- lavInspect(fit2, "fit")["fmin"]*2
          df <- lavInspect(fit2, "fit")["df"]
          c.adj.unstr <- c.adj.val(fit2, Gamma, structured = F)
          c.adj.str <- c.adj.val(fit2, Gamma, structured = T)
          c.adj.unstr.tri <-c.adj(fit2, Gamma.tri, structured = F)
          c.adj.str.tri <-c.adj(fit2, Gamma.tri, structured = T)
          c.all <- c(fmin, df,c.adj.unstr, c.adj.str, c.adj.unstr.tri, c.adj.str.tri)
          
          
          fminB <- lavInspect(fit2, "fit")["baseline.chisq"]/lavInspect(fit2, "fit")["ntotal"]
          dfB <- lavInspect(fit2, "fit")["baseline.df"]
          cB.adj.unstr <- cB.adj.val(fit2, Gamma, structured = F)
          cB.adj.str <- cB.adj.val(fit2, Gamma, structured = T)
          cB.adj.unstr.tri <- cB.adj.val(fit2, Gamma.tri, structured = F)
          cB.adj.str.tri <-cB.adj.val(fit2, Gamma.tri, structured = T)
          cB.all <- c(fminB, dfB, cB.adj.unstr, cB.adj.str,
                      cB.adj.unstr.tri,  cB.adj.str.tri )
          
          
          
          component.all <- c( c.all, cB.all)
          component.matrix[,i] <- component.all 
          
        }
        
        component.matrix <- round(component.matrix, 8)
        component.list[[j]] <- component.matrix
        print(j)
        
      }
    }
  }
  component.list
  
}









##Purpose: generate a list of matrices with confidence interval for RMSEA for a given condition. 
### Each simulated dataset has a matrix of fit indices with rows being different kinds of fit and columns being different path models. 
### These matrices are combined as a list across repetitions. 
####Argument:
#pop.model: the population model that is used for generating data. 
#sat.model: the model used in step 1 of the two-stage procedure; it is saturated in the latent variables. 
#path.model.list: a list of path model used in step 2 of the two-stage procedure
#sample.size: sample size in the condition
#rep.num: number of repetitions for the stimulation study. 
simu.rmsea.ci <- function(pop.model, sat.model, path.model.list, sample.size, rep.num){
  num.fit.indices <- 6
  num.path.mod <- length(path.model.list)
  
  ci.list <- vector(mode="list", length=rep.num)
  struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                   "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                   "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                   "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                   "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) 
  
  for(j in 1:rep.num){
    ci.matrix <- matrix(nrow=num.fit.indices, ncol=num.path.mod)
    colnames(ci.matrix) <-paste("path.mod", 0:(num.path.mod-1), sep="")
    rownames(ci.matrix) <- c("rmsea.ci.lower.default",
                             "rmsea.ci.upper.default", 
                             "rmsea.ci.lower.adj.str.exp",
                             "rmsea.ci.upper.adj.str.exp",
                             "rmsea.ci.lower.adj.str.exp.tri",
                             "rmsea.ci.upper.adj.str.exp.tri")
    
    simuData<- simulateData(pop.model, sample.nobs=sample.size)
    fit1 <- sem(sat.model, data = simuData, estimator="ML", likelihood="wishart")
    if (lavInspect(fit1, "converged")==F){ 
      print(j)
    } else {
      sat.struct.cov <- lavInspect(fit1, "cov.lv")[c(5:7, 1:4),c(5:7, 1:4) ]
      if (is.positive.definite(sat.struct.cov)==F){
        print(j)
      } else{
        #compute gamma
        fit1@Options$h1.information = "structured" 
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.observed")[struct.path, struct.path]
        V1.unstr <- lavInspect(fit1, "information.first.order")[struct.path, struct.path]
        Gamma.tri <- W1.unstr.invert%*%V1.unstr%*%W1.unstr.invert
        Gamma <- W1.unstr.invert
        
        
        
        for(i in 1:num.path.mod){
          fit2 <- sem(path.model.list[[i]], 
                      sample.cov = sat.struct.cov, 
                      sample.nobs = sample.size, 
                      likelihood = "wishart")
          df <- lavInspect(fit2, "fit")["df"]
          
          
          
          
          #compute fit indices
          rmsea.ci.default <- lavInspect(fit2, "fit")[c("rmsea.ci.lower",
                                                        "rmsea.ci.upper")]
          
          c.adj <- c.adj.val(fit2, gamma=Gamma, structured=T, expected=T)/df
          fit.sample.adj <- sem(path.model.list[[i]], 
                                sample.cov = sat.struct.cov, 
                                sample.nobs = sample.size/c.adj, 
                                likelihood = "wishart", start=fit2)
          rmsea.ci.adj.str.exp <- lavInspect(fit.sample.adj, "fit")[c("rmsea.ci.lower",
                                                                      "rmsea.ci.upper")]
          
          
          
          c.adj.tri <- c.adj.val(fit2, gamma=Gamma.tri, structured=T, expected=T)/df
          fit.sample.adj.tri <- sem(path.model.list[[i]], 
                                sample.cov = sat.struct.cov, 
                                sample.nobs = sample.size/c.adj.tri, 
                                likelihood = "wishart", start=fit2)
          rmsea.ci.adj.str.exp.tri <- lavInspect(fit.sample.adj.tri, 
                                                 "fit")[c("rmsea.ci.lower", 
                                                          "rmsea.ci.upper")]
          rmsea.ci.all <- c(rmsea.ci.default, 
                            rmsea.ci.adj.str.exp,
                            rmsea.ci.adj.str.exp.tri)

          ci.matrix[,i] <-   rmsea.ci.all
          
        }
        
        ci.matrix <- round(ci.matrix, 8)
        ci.list[[j]] <- ci.matrix
        print(j)
        
      }
    }
  }
  ci.list
}



##Purpose: generate a list of matrices with confidence interval for SRMR for a given condition. 
### Each simulated dataset has a matrix of fit indices with rows being different kinds of fit and columns being different path models. 
### These matrices are combined as a list across repetitions. 
####Argument:
#pop.model: the population model that is used for generating data. 
#sat.model: the model used in step 1 of the two-stage procedure; it is saturated in the latent variables. 
#path.model.list: a list of path model used in step 2 of the two-stage procedure
#sample.size: sample size in the condition
#rep.num: number of repetitions for the stimulation study. 
simu.srmr.ci <- function(pop.model, sat.model, path.model.list, sample.size, rep.num){
  num.fit.indices <- 6
  num.path.mod <- length(path.model.list)
  
  ci.list <- vector(mode="list", length=rep.num)
  struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                   "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                   "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                   "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                   "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) 
  
  for(j in 1:rep.num){
    ci.matrix <- matrix(nrow=num.fit.indices, ncol=num.path.mod)
    colnames(ci.matrix) <-paste("path.mod", 0:(num.path.mod-1), sep="")
    rownames(ci.matrix) <- c("srmr.ci.lower.default",
                             "srmr.ci.upper.default", 
                             "srmr.ci.lower.adj.str.exp",
                             "srmr.ci.upper.adj.str.exp",
                             "srmr.ci.lower.adj.str.exp.tri",
                             "srmr.ci.upper.adj.str.exp.tri")
    
    simuData<- simulateData(pop.model, sample.nobs=sample.size)
    fit1 <- sem(sat.model, data = simuData, estimator="ML", likelihood="wishart")
    if (lavInspect(fit1, "converged")==F){ 
      print(j)
    } else {
      sat.struct.cov <- lavInspect(fit1, "cov.lv")[c(5:7, 1:4),c(5:7, 1:4) ]
      if (is.positive.definite(sat.struct.cov)==F){
        print(j)
      } else{
        #compute gamma
        fit1@Options$h1.information = "structured" 
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.observed")[struct.path, struct.path]
        V1.unstr <- lavInspect(fit1, "information.first.order")[struct.path, struct.path]
        Gamma.tri <- W1.unstr.invert%*%V1.unstr%*%W1.unstr.invert
        Gamma <- W1.unstr.invert
        
        
        
        for(i in 1:num.path.mod){
          fit2 <- sem(path.model.list[[i]], 
                      sample.cov = sat.struct.cov, 
                      sample.nobs = sample.size, 
                      likelihood = "wishart")
          df <- lavInspect(fit2, "fit")["df"]
          
          
          #compute fit indices
          srmr.ci.default <- srmr.adj.old.ci(fit2)
            
            srmr.ci.adj.str.exp <- srmr.adj.ci(fit2, Gamma)
            
            srmr.ci.adj.str.exp.tri <- srmr.adj.ci(fit2, Gamma.tri)
            
            srmr.ci.all <- c(srmr.ci.default, 
                             srmr.ci.adj.str.exp,
                             srmr.ci.adj.str.exp.tri)
          
          ci.matrix[,i] <-   srmr.ci.all
          
        }
        
        ci.matrix <- round(ci.matrix, 8)
        ci.list[[j]] <- ci.matrix
        print(j)
        
      }
    }
  }
  ci.list
}









##Purpose: compute the coverage rate of the confidence intervals from the simulation study. 
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
ci.coverage <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.default <- ci[1,] <= pop.indices &pop.indices <=ci[2,]
    ci.adj.str.exp <- ci[3,] <= pop.indices &pop.indices <=ci[4,]
    ci.adj.str.exp.tri <- ci[5,] <= pop.indices &pop.indices <=ci[6,]
    ci.check[[i]] <- rbind(ci.default,
                           ci.adj.str.exp,
                           ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}
