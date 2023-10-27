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
    round(apply(simplify2array(lis), 1:2, sd, na.rm = T),3)
  } else{
    lis <- lis[!(sapply(lis, is.null))]
    round(apply(simplify2array(lis), 1:2, sd, na.rm = T),3)
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
    round(sqrt.val,9)
  } else{
    simu.list <- simu.list[!(sapply(simu.list, is.null))]
    to.array <- simplify2array(simu.list)
    subtr.pop <- sweep(to.array, 1:2, pop.matrix)
    sqr.dif  <- subtr.pop ^2
    exp.val <- apply( sqr.dif, 1:2, mean, na.rm=T)
    sqrt.val  <- sqrt(exp.val)
    round(sqrt.val ,9)
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



##Purpose: check SRMR CI computation; prints out the k_s, adjusted SRMR, standard error and the CI
##Argument:
## 1) fit: the SEM fit from step 2
## 2) gamma: a gamma matrix computed using fit1 in step 1; can either be triple product or not
## 3) structured: logical; whether the weight matrix should be structured.
## 4) expected: logical; whether the weight matrix should be expected or observed; only differ if structured=T
srmr.adj.ci.check.val <- function(fit, gamma, structured=T, expected=T){
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
  
  out <- c(k_s, srmr_unbias.new, se.new, srmr_unbias.new.ci.lower, srmr_unbias.new.ci.upper)
  names(out)<- c("k_s", "srmr", "se", "srmr.ci.lower", "srmr.ci.upper")
  out
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





##Purpose: compute the adjust constant for the hypothesized model for computing RMSEA and CFI.
## PLEASE NOTE THAT IN THE PAPER, the adjustment constant is c.adj.val/df.
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
  sqrt(max((Fmin- c.adj/(n-1))/df, 0))
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
  
  if(max(FminB-cB.adj/(n-1), Fmin-c.adj/(n-1), 0) ==0 ){
    1
  } else{
    1- max(Fmin-c.adj/(n-1), 0) / max(FminB-cB.adj/(n-1), Fmin-c.adj/(n-1))}  
}


##Purpose: compute the new adjusted confidence interval for CFI. 
##Argument:
## 1) fit: the SEM fit from step 2
## 2) gamma: a gamma matrix computed using fit1 in step 1; can either be triple product or not
## 3) structured: logical; whether the weight matrix should be structured. 
## 4) expected: logical; whether the expected weight matrix should be used. 
cfi.adj.ci <- function(fit2, gamma, structured=T, expected=T){
  if (structured==T){
    fit2@Options$h1.information = "structured" 
  } else {
    fit2@Options$h1.information = "unstructured" 
  }
  n <- lavInspect(fit2, "nobs")
  fit2B <-  lavaan:::lav_object_independence(fit2, se=T) 
  c.adj <- c.adj.val(fit2, gamma, structured, expected)
  cB.adj <- cB.adj.val(fit2, gamma, structured, expected)
  
  dfB <- lavInspect(fit2B, "fit")
  FminB <- lavInspect(fit2B, "fit")["fmin"]*2
  Fmin <- lavInspect(fit2, "fit")["fmin"]*2
  df <- lavInspect(fit2, "fit")["df"]
  
  if(max(FminB-cB.adj/(n-1), Fmin-c.adj/(n-1), 0) ==0 ){
    cfi.adjust <- 1
  } else{
    cfi.adjust <-1- max(Fmin-c.adj/(n-1), 0) / max(FminB-cB.adj/(n-1), Fmin-c.adj/(n-1))
  }  
  ####Computing CI using Lai(2019)'s method
  #CFI CI based on 
  #Keke Lai (2019) A Simple Analytic Confidence Interval for CFI Given
  #Nonnormal Data, Structural Equation Modeling: A Multidisciplinary Journal, 26:5, 757-777
  #Code is from Appendix B
  #Note: .A indicates things related to the hypothesized model
  #      .Z indicates things related to the baseline model. These are used in Lai's original code
  
  H.A <- inspect(fit2, "hessian")*2
  H.A.inv <- try(chol2inv(chol(H.A)), TRUE)
  if(class(H.A.inv)[1]!="matrix") stop("Model A did not converge to a minimizer")
  H.Z <- inspect(fit2B, "hessian")*2
  H.Z.inv <- try(chol2inv(chol(H.Z)), TRUE)
  N <- inspect(fit2, "nobs") #same as n
  ## Rearrange some matrices ##
  sigma.dev1.A <- lavaan:::computeDelta(fit2@Model)[[1]] #A: hypothesized model
  sigma.dev1.Z0 <- lavaan:::computeDelta(fit2B@Model)[[1]]#Z: baseline model
  S <- inspect(fit2, "sampstat")$'cov'
  p <- dim(S)[1]
  S <- as.matrix(S, p, p)
  Sigma.theta.A <- inspect(fit2, "cov.ov") #note: this is the model-implied variance-covariance matrix of the observed variables. Aliases "sigma", "sigma.hat"
  Sigma.theta.Z0 <- inspect(fit2B, "cov.ov")
  p.star <- p*(p+1)/2
  target.var.names <- rownames(S)
  current.var.names <- rownames(Sigma.theta.Z0)
  Sigma.theta.Z <- matrix(NA, p, p)
  rownames(Sigma.theta.Z) <- colnames(Sigma.theta.Z) <- target.var.names
  for (i.row in 1:p){
    for(i.col in 1:p){
      row.name <- target.var.names[i.row]
      col.name <- target.var.names[i.col]
      pick.row <- which(current.var.names==row.name)
      pick.col <- which(current.var.names==col.name)
      Sigma.theta.Z[i.row, i.col] <- Sigma.theta.Z0[pick.row, pick.col]
    }
  }
  current.matrix <- matrix(NA, p, p)
  current.matrix[lower.tri(current.matrix,
                           diag=TRUE)] <- 1:p.star
  pick.vech <- rep(NA, p.star)
  j <- 1
  for(i.col in 1:p){
    for(i.row in i.col:p){
      row.name <- target.var.names[i.row]
      col.name <- target.var.names[i.col]
      pick.row <- which(current.var.names==row.name)
      pick.col <- which(current.var.names==col.name)
      if(pick.row >= pick.col) pick.vech[j] <- current.matrix[pick.row, pick.col]
      if(pick.row < pick.col) pick.vech[j] <- current.matrix[pick.col, pick.row]
      j <- j+1
    }
  }
  
  q.Z <- dim(sigma.dev1.Z0)[2]
  sigma.dev1.Z <- matrix(NA, p.star, q.Z)
  for(i in 1:p.star){
    pick <- pick.vech[i]
    sigma.dev1.Z[i,] <- sigma.dev1.Z0[pick,]
  }
  ## Finish rearranging matrices ##
  S.inv <- chol2inv(chol(S))
  SS <- S %x% S
  s <- lav_matrix_vech(S)
  D <- lav_matrix_duplication(p)
  Sigma.theta.A.inv <- chol2inv(chol(Sigma.theta.A))
  Sigma.theta.Z.inv <- chol2inv(chol(Sigma.theta.Z))
  W.A <- Sigma.theta.A.inv %x% Sigma.theta.A.inv
  W.Z <- Sigma.theta.Z.inv %x% Sigma.theta.Z.inv
  DWD.A <- lav_matrix_duplication_pre_post(W.A)
  DWD.Z <- lav_matrix_duplication_pre_post(W.Z)
  #Derivatives of F wrt S and theta ('th' in code)
  #The leading 'J' and 'H' denote 1st and 2nd derivatives
  vec.S.inv <- lav_matrix_vec(S.inv)
  S.inv.x.S.inv <- S.inv %x% S.inv
  H.Fs.A <- H.Fs.Z <- lav_matrix_duplication_pre_post(S.inv.x.S.inv)
  vec.Sig.theta.A.inv <- lav_matrix_vec(Sigma.theta.A.inv)
  vec.Sig.theta.Z.inv <- lav_matrix_vec(Sigma.theta.Z.inv)
  J.Fs.A <- t(vec.Sig.theta.A.inv)%*%D - t(vec.S.inv)%*%D #Appendix A2 for hypothesized model
  J.Fs.Z <- t(vec.Sig.theta.Z.inv)%*%D - t(vec.S.inv)%*%D #Appendix A2 for baseline model
  #H.F.ths.A <- (-1)*t(sigma.dev1.A) %*% DWD.A
  # H.F.ths.Z <- (-1)*t(sigma.dev1.Z) %*% DWD.Z
  # H.F.sth.A <- t(H.F.ths.A)
  
  # H.F.sth.Z <- t(H.F.ths.Z)
  # H.F.thth.A <- H.A
  # H.F.thth.Z <- H.Z
  # J.ths.A <- (-1)*H.A.inv%*% H.F.ths.A
  # J.ths.Z <- (-1)*H.Z.inv%*% H.F.ths.Z
  #Full derivative d2.F(theta,s)/ds.ds
  #H.Fs.Full.1.A <- H.Fs.A
  #H.Fs.Full.2.A <- H.F.sth.A %*% J.ths.A
  #H.Fs.Full.3.A <- t(J.ths.A) %*% (H.F.ths.A + H.A %*% J.ths.A)
  # H.Fs.full.A <- H.Fs.Full.1.A + H.Fs.Full.2.A + H.Fs.Full.3.A
  # H.Fs.Full.1.Z <- H.Fs.Z
  # H.Fs.Full.2.Z <- H.F.sth.Z %*% J.ths.Z
  # H.Fs.Full.3.Z <- t(J.ths.Z) %*% (H.F.ths.Z + H.Z %*% J.ths.Z)
  # H.Fs.full.Z <- H.Fs.Full.1.Z + H.Fs.Full.2.Z + H.Fs.Full.3.Z
  F.A.ML <- log(det(Sigma.theta.A)) -log(det(S))+sum(diag(S %*% Sigma.theta.A.inv)) -p
  F.Z.ML <- log(det(Sigma.theta.Z)) -log(det(S))+sum(diag(S %*% Sigma.theta.Z.inv)) -p
  
  J.fs <- rbind(J.Fs.A, J.Fs.Z)
  #H.fs <- matrix(NA, nrow=2*p.star, ncol=p.star)
  # ind.FA <- seq(from=1, to=2*p.star, by=2)
  # ind.FZ <- ind.FA + 1
  # H.fs[ind.FA,] <- H.Fs.full.A
  # H.fs[ind.FZ,] <- H.Fs.full.Z
  #Derivatives of CFI wrt S
  J.CFI.f <- cbind((-1)/F.Z.ML, F.A.ML*F.Z.ML^(-2)) #appendix A15
  # H.CFI.f <- diag(0, 2)
  #H.CFI.f[1,2] <- H.CFI.f[2,1] <- F.Z.ML^(-2)
  #H.CFI.f[2,2] <- (-2)*F.A.ML * F.Z.ML^(-3)
  J.CFI.s <- J.CFI.f %*% J.fs #appendix A14
  #H.CFI.s.1 <- t(J.fs) %*% H.CFI.f %*% J.fs
  # H.CFI.s.2 <- diag(1, p.star) %x% J.CFI.f
  # H.CFI.s.3 <- H.CFI.s.2 %*% H.fs
  #H.CFI.s <- H.CFI.s.1 + H.CFI.s.3
  #SE1 and SE2 for CFI
  V1.CFI <- J.CFI.s %*% gamma %*% t(J.CFI.s)
  # HG.CFI <- H.CFI.s %*% gamma
  #HG2.CFI <- HG.CFI %*% HG.CFI
  #V2.CFI <- V1.CFI + sum(diag(HG2.CFI))/(2*N)
  SE1 <- sqrt(V1.CFI/N) #According to Lai (2019), this one should work in most situations. I will just use this one.
  #SE2 <- sqrt(V2.CFI/N)
  
  
  
  
  alpha <- 0.1
  c<-1-alpha/2
  zc<-qnorm(c) 
  
  cfi.adjust.ci.lower <- min(max(cfi.adjust-zc*SE1,0),1)
  cfi.adjust.ci.upper <- min(cfi.adjust+zc*SE1,1)
  
  
  c(cfi.adjust.ci.lower, cfi.adjust.ci.upper)
  
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
    colnames(fit.matrix) <-paste("path.mod", 1:(num.path.mod), sep="")
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
                              "srmr.default.adj.unstr",
                              "srmr.default.adj.str", 
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
        fit1@Options$h1.information = "unstructured" 
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.expected")[struct.path, struct.path]
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
          srmr.default.adj.unstr <- srmr.adj.old(fit2, structured=F)
          srmr.default.adj.str <- srmr.adj.old(fit2, structured=T)
          srmr.adj.unstr <-srmr.adj(fit2, Gamma, structured = F)
          srmr.adj.str.exp <- srmr.adj(fit2, Gamma, structured = T)
          srmr.adj.str.obs <- srmr.adj(fit2, Gamma, structured = T, expected=F)
          srmr.adj.unstr.tri <-srmr.adj(fit2, Gamma.tri, structured = F)
          srmr.adj.str.exp.tri <-srmr.adj(fit2, Gamma.tri, structured = T)
          srmr.adj.str.obs.tri <-srmr.adj(fit2, Gamma.tri, structured = T, expected = F)
          srmr.all <- c(srmr.unadj,
                        srmr.default.adj.unstr,srmr.default.adj.str, 
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
    colnames(component.matrix) <-paste("path.mod", 1:(num.path.mod), sep="")
    rownames(component.matrix) <- c("rmsea.default",
                                    "rmsea.adj.unstr", "rmsea.adj.str",
                                    "rmsea.adj.unstr.tri",   "rmsea.adj.str.tri", 
                                    "cfi.default", 
                                    "cfi.adj.unstr", "cfi.adj.str",
                                    "cfi.adj.unstr.tri","cfi.adj.str.tri", 
                                    "srmr.unadj", 
                                    "srmr.default.adj.unstr",   "srmr.default.adj.str", 
                                    "srmr.adj.unstr", "srmr.adj.str", 
                                    "srmr.adj.unstr.tri",  "srmr.adj.str.tri")
    
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
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.expected")[struct.path, struct.path]
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
  num.fit.indices <- 10
  num.path.mod <- length(path.model.list)
  
  ci.list <- vector(mode="list", length=rep.num)
  struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                   "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                   "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                   "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                   "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) 
  
  for(j in 1:rep.num){
    ci.matrix <- matrix(nrow=num.fit.indices, ncol=num.path.mod)
    colnames(ci.matrix) <-paste("path.mod", 1:(num.path.mod), sep="")
    rownames(ci.matrix) <- c("rmsea.ci.lower.default",
                             "rmsea.ci.upper.default", 
                             "rmsea.ci.lower.adj.unstr",
                             "rmsea.ci.upper.adj.unstr",
                             "rmsea.ci.lower.adj.str.exp",
                             "rmsea.ci.upper.adj.str.exp",
                             "rmsea.ci.lower.adj.unstr.tri",
                             "rmsea.ci.upper.adj.unstr.tri",
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
        fit1@Options$h1.information = "unstructured" 
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.expected")[struct.path, struct.path]
        V1.unstr <- lavInspect(fit1, "information.first.order")[struct.path, struct.path]
        Gamma.tri <- W1.unstr.invert%*%V1.unstr%*%W1.unstr.invert
        Gamma <- W1.unstr.invert
        
        
        
        for(i in 1:num.path.mod){
          fit2 <- sem(path.model.list[[i]], 
                      sample.cov = sat.struct.cov, 
                      sample.nobs = sample.size, 
                      likelihood = "wishart",check.vcov = FALSE)
          df <- lavInspect(fit2, "fit")["df"]
          
          
          
          
          #compute fit indices
          rmsea.ci.default <- lavInspect(fit2, "fit")[c("rmsea.ci.lower",
                                                        "rmsea.ci.upper")]
          
          c.adj.unstr <- c.adj.val(fit2, Gamma, structured = F)/df #note that in the paper, c=c.adj/df
          
          fit.sample.adj.unstr <- sem(path.model.list[[i]], 
                                      sample.cov = sat.struct.cov, 
                                      sample.nobs = (sample.size-1)/c.adj.unstr+1, 
                                      likelihood = "wishart", start=fit2,check.vcov = FALSE)
          rmsea.ci.adj.unstr <- lavInspect(fit.sample.adj.unstr, "fit")[c("rmsea.ci.lower",
                                                                          "rmsea.ci.upper")]
          
          c.adj.str.exp <- c.adj.val(fit2, gamma=Gamma, structured=T, expected=T)/df
          fit.sample.adj.str.exp <- sem(path.model.list[[i]], 
                                        sample.cov = sat.struct.cov, 
                                        sample.nobs = (sample.size-1)/c.adj.str.exp+1, 
                                        likelihood = "wishart", start=fit2,check.vcov = FALSE)
          rmsea.ci.adj.str.exp <- lavInspect( fit.sample.adj.str.exp, "fit")[c("rmsea.ci.lower",
                                                                               "rmsea.ci.upper")]
          
          c.adj.unstr.tri <- c.adj.val(fit2, gamma=Gamma.tri, structured=F)/df
          fit.sample.adj.unstr.tri <- sem(path.model.list[[i]], 
                                          sample.cov = sat.struct.cov, 
                                          sample.nobs = (sample.size-1)/c.adj.unstr.tri+1, 
                                          likelihood = "wishart", start=fit2,check.vcov = FALSE)
          rmsea.ci.adj.unstr.tri <- lavInspect(fit.sample.adj.unstr.tri, "fit")[c("rmsea.ci.lower",
                                                                                  "rmsea.ci.upper")]
          
          c.adj.str.exp.tri <- c.adj.val(fit2, gamma=Gamma.tri, structured=T, expected=T)/df
          fit.sample.adj.str.exp.tri <- sem(path.model.list[[i]], 
                                            sample.cov = sat.struct.cov, 
                                            sample.nobs = (sample.size-1)/c.adj.str.exp.tri+1,
                                            likelihood = "wishart", start=fit2,check.vcov = FALSE)
          rmsea.ci.adj.str.exp.tri <- lavInspect(fit.sample.adj.str.exp.tri, 
                                                 "fit")[c("rmsea.ci.lower", 
                                                          "rmsea.ci.upper")]
          
          rmsea.ci.all <- c(rmsea.ci.default, 
                            rmsea.ci.adj.unstr,
                            rmsea.ci.adj.str.exp,
                            rmsea.ci.adj.unstr.tri,
                            rmsea.ci.adj.str.exp.tri)
          
          ci.matrix[,i] <-   rmsea.ci.all
          
        }
        
        ci.matrix <- round(ci.matrix, 9)
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
  num.fit.indices <- 12
  num.path.mod <- length(path.model.list)
  
  ci.list <- vector(mode="list", length=rep.num)
  struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                   "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                   "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                   "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                   "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) 
  
  for(j in 1:rep.num){
    ci.matrix <- matrix(nrow=num.fit.indices, ncol=num.path.mod)
    colnames(ci.matrix) <-paste("path.mod", 1:(num.path.mod), sep="")
    rownames(ci.matrix) <- c("srmr.ci.lower.default.adj.unstr",
                             "srmr.ci.upper.default.adj.unstr",
                             "srmr.ci.lower.default.adj.str",
                             "srmr.ci.upper.default.adj.str", 
                             "srmr.ci.lower.adj.unstr",
                             "srmr.ci.upper.adj.unstr",
                             "srmr.ci.lower.adj.str.exp",
                             "srmr.ci.upper.adj.str.exp", 
                             "srmr.ci.lower.adj.unstr.tri",
                             "srmr.ci.upper.adj.unstr.tri",
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
        fit1@Options$h1.information = "unstructured" 
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.expected")[struct.path, struct.path]
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
          srmr.ci.default.adj.unstr <- srmr.adj.old.ci(fit2,structured = F)
          srmr.ci.default.adj.str <- srmr.adj.old.ci(fit2)
          srmr.ci.adj.unstr <- srmr.adj.ci(fit2, gamma=Gamma, structured = F)
          srmr.ci.adj.str.exp <- srmr.adj.ci(fit2, gamma=Gamma)
          srmr.ci.adj.unstr.tri <- srmr.adj.ci(fit2, gamma=Gamma.tri, structured = F)
          srmr.ci.adj.str.exp.tri <- srmr.adj.ci(fit2, gamma=Gamma.tri)
          
          srmr.ci.all <- c(srmr.ci.default.adj.unstr,
                           srmr.ci.default.adj.str,
                           srmr.ci.adj.unstr,
                           srmr.ci.adj.str.exp,
                           srmr.ci.adj.unstr.tri,
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


##Purpose: generate a list of matrices that can check the computations of SRMR CIs. 
### Each simulated dataset has a matrix of fit indices with rows being different kinds of fit and columns being different path models. 
### These matrices are combined as a list across repetitions. 
####Argument:
#pop.model: the population model that is used for generating data. 
#sat.model: the model used in step 1 of the two-stage procedure; it is saturated in the latent variables. 
#path.model.list: a list of path model used in step 2 of the two-stage procedure
#sample.size: sample size in the condition
#rep.num: number of repetitions for the stimulation study. 
simu.srmr.ci.check.val <- function(pop.model, sat.model, path.model.list, sample.size, rep.num){
  num.fit.indices <- 10
  num.path.mod <- length(path.model.list)
  
  ci.list <- vector(mode="list", length=rep.num)
  struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                   "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                   "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                   "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                   "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) 
  
  for(j in 1:rep.num){
    ci.matrix <- matrix(nrow=num.fit.indices, ncol=num.path.mod)
    colnames(ci.matrix) <-paste("path.mod", 1:(num.path.mod), sep="")
    rownames(ci.matrix) <- c("k_s.adj.str.exp", 
                             "srmr.adj.str.exp", 
                             "se.adj.str.exp", 
                             "srmr.ci.lower.adj.str.exp", 
                             "srmr.ci.upper.adj.str.exp",
                             "k_s.adj.str.exp.tri", 
                             "srmr.adj.str.exp.tri", 
                             "se.adj.str.exp.tri", 
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
        fit1@Options$h1.information = "unstructured" 
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.expected")[struct.path, struct.path]
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
          
          
          srmr.ci.adj.str.exp <- srmr.adj.ci.check.val(fit2, Gamma)
          
          srmr.ci.adj.str.exp.tri <- srmr.adj.ci.check.val(fit2, Gamma.tri)
          
          srmr.ci.all <- c(srmr.ci.adj.str.exp,
                           srmr.ci.adj.str.exp.tri)
          
          ci.matrix[,i] <-   srmr.ci.all
          
        }
        
        ci.matrix <- round(ci.matrix, 3)
        ci.list[[j]] <- ci.matrix
        print(j)
        
      }
    }
  }
  ci.list
}










##Purpose: generate a list of matrices with confidence interval for CFI for a given condition. 
### Each simulated dataset has a matrix of fit indices with rows being different kinds of fit and columns being different path models. 
### These matrices are combined as a list across repetitions. 
####Argument:
#pop.model: the population model that is used for generating data. 
#sat.model: the model used in step 1 of the two-stage procedure; it is saturated in the latent variables. 
#path.model.list: a list of path model used in step 2 of the two-stage procedure
#sample.size: sample size in the condition
#rep.num: number of repetitions for the stimulation study. 
simu.cfi.ci <- function(pop.model, sat.model, path.model.list, sample.size, rep.num){
  num.fit.indices <- 8
  num.path.mod <- length(path.model.list)
  
  ci.list <- vector(mode="list", length=rep.num)
  struct.path <- c("eta1~~eta1","eta1~~eta2", "eta1~~eta3", "xi1~~eta1", "xi2~~eta1",  "xi3~~eta1", 
                   "xi4~~eta1", "eta2~~eta2", "eta2~~eta3", "xi1~~eta2" , "xi2~~eta2" , "xi3~~eta2", 
                   "xi4~~eta2" , "eta3~~eta3" ,"xi1~~eta3",  "xi2~~eta3" , "xi3~~eta3" , "xi4~~eta3", 
                   "xi1~~xi1" ,  "xi1~~xi2" ,  "xi1~~xi3",   "xi1~~xi4" ,  "xi2~~xi2" ,  "xi2~~xi3"  ,
                   "xi2~~xi4" ,  "xi3~~xi3" ,  "xi3~~xi4"  , "xi4~~xi4"  ) 
  
  for(j in 1:rep.num){
    ci.matrix <- matrix(nrow=num.fit.indices, ncol=num.path.mod)
    colnames(ci.matrix) <-paste("path.mod", 1:(num.path.mod), sep="")
    rownames(ci.matrix) <- c("cfi.ci.lower.adj.unstr",
                             "cfi.ci.upper.adj.unstr",
                             "cfi.ci.lower.adj.str.exp",
                             "cfi.ci.upper.adj.str.exp",
                             "cfi.ci.lower.adj.unstr.tri",
                             "cfi.ci.upper.adj.unstr.tri",
                             "cfi.ci.lower.adj.str.exp.tri",
                             "cfi.ci.upper.adj.str.exp.tri")
    
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
        W1.unstr.invert <- lavInspect(fit1, "inverted.information.expected")[struct.path, struct.path]
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
          cfi.ci.adj.unstr <- cfi.adj.ci(fit2, Gamma, structured=T)
          
          cfi.ci.adj.str.exp <- cfi.adj.ci(fit2, Gamma)
          cfi.ci.adj.unstr.tri <- cfi.adj.ci(fit2, Gamma.tri)
          
          cfi.ci.adj.str.exp.tri <- cfi.adj.ci(fit2, Gamma.tri)
          
          cfi.ci.all <- c(   cfi.ci.adj.unstr,
                             cfi.ci.adj.str.exp, 
                             cfi.ci.adj.unstr.tri,
                             cfi.ci.adj.str.exp.tri )
          
          ci.matrix[,i] <-   cfi.ci.all
          
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
ci.coverage.rmsea <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.default <- ci[1,] <= pop.indices &pop.indices <=ci[2,]
    ci.adj.unstr <- ci[3,] <= pop.indices &pop.indices <=ci[4,]
    ci.adj.str.exp <- ci[5,] <= pop.indices &pop.indices <=ci[6,]
    ci.adj.unstr.tri <- ci[7,] <= pop.indices &pop.indices <=ci[8,]
    ci.adj.str.exp.tri <- ci[9,] <= pop.indices &pop.indices <=ci[10,]
    ci.check[[i]] <- rbind(ci.default,
                           ci.adj.unstr,
                           ci.adj.str.exp,
                           ci.adj.unstr.tri,
                           ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}





##Purpose: compute the coverage rate of the confidence intervals from the simulation study. 
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
ci.coverage.srmr <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.default.unstr <- ci[1,] <= pop.indices &pop.indices <=ci[2,]
    ci.default.str <- ci[3,] <= pop.indices &pop.indices <=ci[4,]
    ci.adj.unstr <- ci[5,] <= pop.indices &pop.indices <=ci[6,]
    ci.adj.str.exp <- ci[7,] <= pop.indices &pop.indices <=ci[8,]
    ci.adj.unstr.tri <- ci[9,] <= pop.indices &pop.indices <=ci[10,]
    ci.adj.str.exp.tri <- ci[11,] <= pop.indices &pop.indices <=ci[12,]
    ci.check[[i]] <- rbind(ci.default.unstr,
                           ci.default.str,
                           ci.adj.unstr,
                           ci.adj.str.exp,
                           ci.adj.unstr.tri,
                           ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}







##Purpose: compute the coverage rate of the confidence intervals from the simulation study. 
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
ci.coverage.cfi <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.adj.unstr <- ci[1,] <= pop.indices &pop.indices <=ci[2,]
    ci.adj.str.exp <- ci[3,] <= pop.indices &pop.indices <=ci[4,]
    ci.adj.unstr.tri <- ci[5,] <= pop.indices &pop.indices <=ci[6,]
    ci.adj.str.exp.tri <- ci[7,] <= pop.indices &pop.indices <=ci[8,]
    ci.check[[i]] <- rbind(ci.adj.unstr,
                           ci.adj.str.exp,
                           ci.adj.unstr.tri,
                           ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}




##Purpose: compute the confidence interval's missing rate from below (i.e., the population value is lower than the lower bound of the CI)
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
missing.rate.below.rmsea <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.default <- ci[1,] >= pop.indices
    ci.adj.unstr <- ci[3,] >= pop.indices 
    ci.adj.str.exp <- ci[5,] >= pop.indices 
    ci.adj.unstr.tri <- ci[7,] >= pop.indices 
    ci.adj.str.exp.tri <- ci[9,] >= pop.indices 
    ci.check[[i]] <- rbind(ci.default,
                           ci.adj.unstr,
                           ci.adj.str.exp,
                           ci.adj.unstr.tri,
                           ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}




##Purpose: compute the confidence interval's missing rate from above (i.e., the population value is higher than the upper bound of the CI)
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
missing.rate.above.rmsea <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.default <- pop.indices >=ci[2,]
    ci.adj.unstr <- pop.indices >=ci[4,]
    ci.adj.str.exp <- pop.indices >=ci[6,]
    ci.adj.unstr.tri <- pop.indices >=ci[8,]
    ci.adj.str.exp.tri <- pop.indices >=ci[10,]
    ci.check[[i]] <- rbind(ci.default,
                           ci.adj.unstr,
                           ci.adj.str.exp,
                           ci.adj.unstr.tri,
                           ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}





##Purpose: compute the confidence interval's missing rate from below for CFI (i.e., the population value is lower than the lower bound of the CI)
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
missing.rate.below.cfi <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.adj.unstr <- ci[1,] >= pop.indices
    ci.adj.str.exp <- ci[3,] >= pop.indices 
    ci.adj.unstr.tri <- ci[5,] >= pop.indices 
    ci.adj.str.exp.tri <- ci[7,] >= pop.indices 
    ci.check[[i]] <- 
      rbind(ci.adj.unstr,
            ci.adj.str.exp, 
            ci.adj.unstr.tri,
            ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}




##Purpose: compute the confidence interval's missing rate from above for CFI (i.e., the population value is higher than the upper bound of the CI)
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
missing.rate.above.cfi <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.adj.unstr <- pop.indices >=ci[2,]
    ci.adj.str.exp <- pop.indices >=ci[4,]
    ci.adj.unstr.tri <- pop.indices >=ci[6,]
    ci.adj.str.exp.tri <- pop.indices >=ci[8,]
    ci.check[[i]] <- 
      rbind(ci.adj.unstr,
            ci.adj.str.exp, 
            ci.adj.unstr.tri,
            ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}




##Purpose: compute the confidence interval's missing rate from below (i.e., the population value is lower than the lower bound of the CI)
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
missing.rate.below.srmr <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.default.unstr <- ci[1,] >= pop.indices
    ci.default.str <- ci[3,] >= pop.indices
    ci.adj.unstr <- ci[5,] >= pop.indices 
    ci.adj.str.exp <- ci[7,] >= pop.indices 
    ci.adj.unstr.tri <- ci[9,] >= pop.indices 
    ci.adj.str.exp.tri <- ci[11,] >= pop.indices 
    ci.check[[i]] <- rbind(ci.default.unstr,
                           ci.default.str,
                           ci.adj.unstr,
                           ci.adj.str.exp,
                           ci.adj.unstr.tri,
                           ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}




##Purpose: compute the confidence interval's missing rate from above (i.e., the population value is higher than the upper bound of the CI)
##Arguments:
##### 1) pop.indices: a numeric vector that contains the population values of the fit indices.
##### 2) ci.data.list: a list of confidence intervals from the simulation. 
missing.rate.above.srmr <- function(pop.indices, ci.data.list){
  ci.check <- list()
  ci.data <- remove.null(ci.data.list)
  simu.num <- length(ci.data)
  for(i in 1:simu.num){
    ci <- ci.data[[i]]
    ci.default.str <- pop.indices >=ci[2,]
    ci.default.unstr <- pop.indices >=ci[4,]
    ci.adj.unstr <- pop.indices >=ci[6,]
    ci.adj.str.exp <- pop.indices >=ci[8,]
    ci.adj.unstr.tri <- pop.indices >=ci[10,]
    ci.adj.str.exp.tri <- pop.indices >=ci[12,]
    ci.check[[i]] <- rbind(ci.default.unstr,
                           ci.default.str,
                           ci.adj.unstr,
                           ci.adj.str.exp,
                           ci.adj.unstr.tri,
                           ci.adj.str.exp.tri)
  }
  
  list.mean(ci.check)
  
}
