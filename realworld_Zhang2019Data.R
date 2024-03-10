library(tidyverse)
library(lavaan)
dat <- read.csv("zhang2019data.csv")




full.mod <- "
conscientiousness =~ CO_I1+CO_I2+CO_I3+CO_I4+CO_I5+CO_I6+CO_I7+CO_I8+CO_I9
optimism =~ LOT_I1+LOT_I2+LOT_I3+LOT_I4+LOT_I5+LOT_I6+LOT_I7+LOT_I8+LOT_I9+LOT_I10
selfEsteem =~ RS_I1+RS_I2+RS_I3+RS_I4+RS_I5+RS_I6+RS_I7+RS_I8+RS_I9+RS_I10
grade=~1*GPA

#regression
grade ~ conscientiousness 
selfEsteem ~ grade
optimism ~ selfEsteem
"

sat.struct <- "
conscientiousness =~ CO_I1+CO_I2+CO_I3+CO_I4+CO_I5+CO_I6+CO_I7+CO_I8+CO_I9
optimism =~ LOT_I1+LOT_I2+LOT_I3+LOT_I4+LOT_I5+LOT_I6+LOT_I7+LOT_I8+LOT_I9+LOT_I10
selfEsteem =~ RS_I1+RS_I2+RS_I3+RS_I4+RS_I5+RS_I6+RS_I7+RS_I8+RS_I9+RS_I10
grade=~1*GPA 

#regression
grade ~~ optimism + selfEsteem + conscientiousness 
optimism ~~ selfEsteem + conscientiousness
selfEsteem ~~ conscientiousness 
"

path.mod <- "
grade ~ conscientiousness 
selfEsteem ~ grade
optimism ~ selfEsteem
"

struct.path <- c(
  "conscientiousness~~conscientiousness","conscientiousness~~selfEsteem", "conscientiousness~~grade" ,  "conscientiousness~~optimism" ,
  "selfEsteem~~selfEsteem", "selfEsteem~~grade", "optimism~~selfEsteem" , 
  "grade~~grade", "optimism~~grade",
  "optimism~~optimism"
)
struct.var <- c("conscientiousness", "selfEsteem", "grade", "optimism")


fit <- sem(full.mod, data = dat)
lavInspect(fit, "fit")[c("rmsea", "cfi", "srmr")]
summary(fit)
n <- nrow(dat)

# fit the model using the two-stage procedure
fit1 <- sem(sat.struct, data = dat, estimator = "ML", likelihood = "wishart")

sat.struct.cov <- lavInspect(fit1, "cov.lv")[struct.var, struct.var]

fit2 <- sem(path.mod,
            sample.cov = sat.struct.cov, sample.nobs = n,
            likelihood = "wishart")

summary(fit2)
lavInspect(fit2, "fit")[c("rmsea", "cfi", "srmr")]


##################### Estimating Structural Fit Indices ###############################

## computing structural RMSEA
fit1@Options$h1.information <- "structured"
W1.str.invert <- lavInspect(fit1, "inverted.information.expected")
V1.str <- lavInspect(fit1, "information.first.order")
rownames(W1.str.invert)
# Gamma for V2
Gamma <- W1.str.invert[struct.path, struct.path]
round(Gamma, 2)
# Gamma for V5
Gamma.tri <- (W1.str.invert %*% V1.str %*% W1.str.invert)[struct.path, struct.path]



fit2@Options$h1.information <- "structured"
W2 <- lavInspect(fit2, "h1.information.expected")
deltabreve <- lavInspect(fit2, "delta") # jacobian matrix
Uc <-
  W2 - W2 %*%
  deltabreve %*%
  solve(t(deltabreve) %*%
          W2 %*% deltabreve) %*%
  t(deltabreve) %*% W2

Fmin <- lavInspect(fit2, "fit")["fmin"] * 2
df <- lavInspect(fit2, "fit")["df"]

# small-sample corrections for the hypothesized structural model
c.adj.v2 <- lav_matrix_trace(Uc %*% Gamma) / df
c.adj.v5 <- lav_matrix_trace(Uc %*% Gamma.tri) / df
rmsea.v2 <- sqrt(max(Fmin / df - c.adj.v2 / (n - 1), 0))
rmsea.v5 <- sqrt(max(Fmin / df - c.adj.v5 / (n - 1), 0))
rmsea.v2
rmsea.v5


## computing structural CFI
fit2B <- lavaan:::lav_object_independence(fit2, se = T)
deltabreveB <- lavInspect(fit2B, "delta")
W2B <- lavInspect(fit2B, "h1.information.expected")

Uc.B <-
  W2B - W2B %*%
  deltabreveB %*%
  solve(t(deltabreveB) %*%
          W2B %*% deltabreveB) %*%
  t(deltabreveB) %*% W2B


dfB <- lavInspect(fit2B, "fit")["df"]
FminB <- lavInspect(fit2B, "fit")["fmin"] * 2

# small-sample corrections for the baseline model
cB.adj.v2 <- lav_matrix_trace(Uc.B %*% Gamma) / dfB
cB.adj.v5 <- lav_matrix_trace(Uc.B %*% Gamma.tri) / dfB


cfi.v2 <- 1 - max(Fmin - c.adj.v2 * df / (n - 1), 0) / max(FminB - cB.adj.v2 * dfB / (n - 1), Fmin - c.adj.v2 * df / (n - 1))
cfi.v5 <- 1 - max(Fmin - c.adj.v5 * df / (n - 1), 0) / max(FminB - cB.adj.v5 * dfB / (n - 1), Fmin - c.adj.v5 * df / (n - 1))
cfi.v2
cfi.v5


## computing structural SRMR
residual.cov <- lavInspect(fit2, "resid")$cov

samp.var <- diag(sat.struct.cov)
resid.cov.vech <- OpenMx::vech(residual.cov)
p <- ncol(residual.cov)
t <- p * (p + 1) / 2

s1 <- OpenMx::vech(samp.var %*% t(samp.var))
G <- diag(s1) # G matrix in Maydeu-Olivares 2017's paper
G.inv <- diag(1 / s1)
T_s <- t(resid.cov.vech) %*% G.inv %*% resid.cov.vech


H <- deltabreve %*% solve(t(deltabreve) %*% W2 %*% deltabreve) %*% t(deltabreve) %*% W2
G_sqrt_inv <- diag(1 / sqrt(s1))
e_s <- G_sqrt_inv %*% resid.cov.vech


Xi_u.v2 <- (diag(nrow(H)) - H) %*% Gamma %*% t((diag(nrow(H)) - H)) / n
Xi_s.v2 <- G_sqrt_inv %*% Xi_u.v2 %*% G_sqrt_inv
Xi_u.v5 <- (diag(nrow(H)) - H) %*% Gamma.tri %*% t((diag(nrow(H)) - H)) / n
Xi_s.v5 <- G_sqrt_inv %*% Xi_u.v5 %*% G_sqrt_inv

# small-sample corrections for SRMR
k_s.v2 <- 1 - ((sum(Xi_s.v2^2) + 2 * t(e_s) %*% Xi_s.v2 %*% e_s) / (4 * T_s^2))
tr.Xi.v2 <- sum(diag(Xi_s.v2))
k_s.v5 <- 1 - ((sum(Xi_s.v5^2) + 2 * t(e_s) %*% Xi_s.v5 %*% e_s) / (4 * T_s^2))
tr.Xi.v5 <- sum(diag(Xi_s.v5))


srmr.v2 <- sqrt(max(0, (T_s - tr.Xi.v2) / t)) / k_s.v2
srmr.v5 <- sqrt(max(0, (T_s - tr.Xi.v5) / t)) / k_s.v5

srmr.v2
srmr.v5


all.structural.fit <- c(
  rmsea.v2, rmsea.v5,
  cfi.v2, cfi.v5,
  srmr.v2, srmr.v5
)
names(all.structural.fit) <- c(
  "rmsea.v2", "rmsea.v5",
  "cfi.v2", "cfi.v5",
  "srmr.v2", "srmr.v5"
)
all.structural.fit




############# Computing the Confidence Intervals (CIs) for the Structural Fit Indices ##############



## computing CI for structural RMSEA

# adjust the sample.nobs argument so that the estimated noncentrality parameter takes into account the small-sample corrections (c.adj.v2 or c.adj.v5)
fit.rmsea.v2 <- sem(path.mod,
                    sample.cov = sat.struct.cov,
                    sample.nobs = (n - 1) / c.adj.v2 + 1,
                    likelihood = "wishart", start = fit2
)
fit.rmsea.v5 <- sem(path.mod,
                    sample.cov = sat.struct.cov,
                    sample.nobs = (n - 1) / c.adj.v5 + 1,
                    likelihood = "wishart", start = fit2
)

rmsea.ci.v2 <- lavInspect(fit.rmsea.v2, "fit")[c("rmsea.ci.lower", "rmsea.ci.upper")]
rmsea.ci.v5 <- lavInspect(fit.rmsea.v5, "fit")[c("rmsea.ci.lower", "rmsea.ci.upper")]
rmsea.ci.v2
rmsea.ci.v5


## compute CI for structural CFI

H <- inspect(fit2, "hessian") * 2
H.inv <- solve(H)
H.B <- inspect(fit2B, "hessian") * 2
H.B.inv <- solve(H.B)

# rearranging some matrices
S <- sat.struct.cov
p <- dim(S)[1]
S <- as.matrix(S, p, p)
Sigma.theta <- inspect(fit2, "cov.ov")
Sigma.theta.B0 <- inspect(fit2B, "cov.ov")

target.var.names <- rownames(S)
current.var.names <- rownames(Sigma.theta.B0)
Sigma.theta.B <- matrix(NA, p, p)
rownames(Sigma.theta.B) <- colnames(Sigma.theta.B) <- target.var.names
for (i.row in 1:p) {
  for (i.col in 1:p) {
    row.name <- target.var.names[i.row]
    col.name <- target.var.names[i.col]
    pick.row <- which(current.var.names == row.name)
    pick.col <- which(current.var.names == col.name)
    Sigma.theta.B[i.row, i.col] <- Sigma.theta.B0[pick.row, pick.col]
  }
}
current.matrix <- matrix(NA, p, p)
current.matrix[lower.tri(current.matrix,
                         diag = TRUE
)] <- 1:t
pick.vech <- rep(NA, t)
j <- 1
for (i.col in 1:p) {
  for (i.row in i.col:p) {
    row.name <- target.var.names[i.row]
    col.name <- target.var.names[i.col]
    pick.row <- which(current.var.names == row.name)
    pick.col <- which(current.var.names == col.name)
    if (pick.row >= pick.col) pick.vech[j] <- current.matrix[pick.row, pick.col]
    if (pick.row < pick.col) pick.vech[j] <- current.matrix[pick.col, pick.row]
    j <- j + 1
  }
}

q.B <- dim(deltabreveB)[2]
deltabreve.B <- matrix(NA, t, q.B)
for (i in 1:t) {
  pick <- pick.vech[i]
  deltabreve.B[i, ] <- deltabreveB[pick, ]
}


S.inv <- solve(S)
SS <- S %x% S
s <- lav_matrix_vech(S)
D <- lav_matrix_duplication(p)
Sigma.theta.inv <- solve(Sigma.theta)
Sigma.theta.B.inv <- solve(Sigma.theta.B)
W <- Sigma.theta.inv %x% Sigma.theta.inv
W.B <- Sigma.theta.B.inv %x% Sigma.theta.B.inv
DWD <- lav_matrix_duplication_pre_post(W)
DWD.B <- lav_matrix_duplication_pre_post(W.B)

vec.S.inv <- lav_matrix_vec(S.inv)
S.inv.x.S.inv <- S.inv %x% S.inv
H.Fs <- H.Fs.B <- lav_matrix_duplication_pre_post(S.inv.x.S.inv)
vec.Sig.theta.inv <- lav_matrix_vec(Sigma.theta.inv)
vec.Sig.theta.B.inv <- lav_matrix_vec(Sigma.theta.B.inv)
J.Fs <- t(vec.Sig.theta.inv) %*% D - t(vec.S.inv) %*% D
J.Fs.B <- t(vec.Sig.theta.B.inv) %*% D - t(vec.S.inv) %*% D

F.ML <- log(det(Sigma.theta)) - log(det(S)) + sum(diag(S %*% Sigma.theta.inv)) - p
F.B.ML <- log(det(Sigma.theta.B)) - log(det(S)) + sum(diag(S %*% Sigma.theta.B.inv)) - p
J.fs <- rbind(J.Fs, J.Fs.B)
J.CFI.f <- cbind((-1) / F.B.ML, F.ML * F.B.ML^(-2))
J.CFI.s <- J.CFI.f %*% J.fs

# SE for CFI
se.cfi.v2 <- sqrt((J.CFI.s %*% Gamma %*% t(J.CFI.s)) / n)
se.cfi.v5 <- sqrt((J.CFI.s %*% Gamma.tri %*% t(J.CFI.s)) / n)


alpha <- 0.1
c <- 1 - alpha / 2
zc <- qnorm(c)
zc

cfi.ci.v2 <- c(
  min(max(cfi.v2 - zc * se.cfi.v2, 0), 1),
  min(cfi.v2 + zc * se.cfi.v2, 1)
)
cfi.ci.v5 <- c(
  min(max(cfi.v5 - zc * se.cfi.v5, 0), 1),
  min(cfi.v5 + zc * se.cfi.v5, 1)
)
cfi.ci.v2
cfi.ci.v5


## compute CI for structural SRMR
tr.Xissqr.v2 <- sum(Xi_s.v2^2)
tr.Xissqr.v5 <- sum(Xi_s.v5^2)

se.srmr.v2 <- sqrt(k_s.v2^(-2) * (tr.Xissqr.v2 + 2 * t(e_s) %*% Xi_s.v2 %*% e_s) / (2 * t * T_s))
se.srmr.v5 <- sqrt(k_s.v5^(-2) * (tr.Xissqr.v5 + 2 * t(e_s) %*% Xi_s.v5 %*% e_s) / (2 * t * T_s))

srmr.ci.v2 <- c(
  max(0, srmr.v2 - zc * se.srmr.v2),
  srmr.v2 + zc * se.srmr.v2
)
srmr.ci.v5 <- c(
  max(0, srmr.v5 - zc * se.srmr.v5),
  srmr.v5 + zc * se.srmr.v5
)
srmr.ci.v2
srmr.ci.v5

all.ci <- c(
  rmsea.ci.v2, rmsea.ci.v5,
  cfi.ci.v2, cfi.ci.v5,
  srmr.ci.v2, srmr.ci.v5
)
names(all.ci) <- c(
  "rmsea.ci.v2.lower", "rmsea.ci.v2.upper",
  "rmsea.ci.v5.lower", "rmsea.ci.v5.upper",
  "cfi.ci.v2.lower", "cfi.ci.v2.upper",
  "cfi.ci.v5.lower", "cfi.ci.v5.upper",
  "srmr.ci.v2.lower", "srmr.ci.v2.upper",
  "srmr.ci.v5.lower", "srmr.ci.v5.upper"
)
all.ci
all.structural.fit

