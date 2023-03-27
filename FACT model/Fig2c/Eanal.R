# Get rid of the outlier
Egood <- subset(Eresults,Sample != "F8")

# Set up the first function
dk <- 1
rhocost <- 0
Efitfun4 <- function(x) {
 rhoR <- x[1]
 rhoS <- x[2]
 aR <- x[3]
 aS <- x[4]
 with(Egood,(1+rhocost*abs(rhoS)*(rhoS<0))*
      sum((Estradiol-(rhoR*R+rhoS*S)/(dk+aR*R+aS*S))^2))
}
x4 <- c(100,0,1,1)
rhocost <- 0
Efit4 <- optim(x4,Efitfun4)
#### With LLOD = 0 #####
# It wants a negative value of rhoS which does give a better fit 
# 98.71836 -1.12578  1.00715  0.29052
# 10679
# rhocost <- 1
# Efit4 <- optim(x4,Efitfun4)
# 1.0094e+02 1.2545e-07 1.0303e+00 3.2193e-01
# 10986

#### With LLOD = 35 #####
# 85.72565  4.27545  0.88123  0.29326
# 8463.6

#### With LLOD = 62.5 #####
# 83.58364 10.11159  0.86354  0.34527
# 7321.1

Efitfun3 <- function(x) {
 rhoR <- x[1]
 aR <- x[2]
 aS <- x[3]
 with(Egood,sum((Estradiol-rhoR*R/(dk+aR*R+aS*S))^2))
}
x3 <- c(100,1,1)
Efit3 <- optim(x3,Efitfun3)
# 98.69654  1.00744  0.31475
# 10986
#### With LLOD = 35 ####
# 91.75880  0.94395  0.21866
# 12752
# Much worse
#### With LLOD = 62.5 ###
# 98.40648  1.02700  0.16222
# 24466
# Way worse

# Keep rhoS
Efitfun2 <- function(x) {
 rhoR <- x[1]
 rhoS <- x[2]
 aS <- x[3]
 with(Egood,sum((Estradiol-(rhoR*R+rhoS*S)/(R+aS*S))^2))
}
x2 <- c(100,5,1)
Efit2 <- optim(x2,Efitfun2)
# 97.86384  0.29972
# 11776
#### With LLOD = 35 ####
# 97.27822  4.85080  0.33273
# 8463.5
# Nice.
#### With LLOD = 62.5 ###
# 96.7872 11.7068  0.3998
# 7321.1
# Nice.

Efun2 <- function(S,R) {
 rhoR <- Efit2$par[1]
 rhoS <- Efit2$par[2]
 aS <- Efit2$par[3]
 (rhoR*R+rhoS*S)/(R+aS*S)
}

Efun3 <- function(S,R) {
 rhoR <- Efit3$par[1]
 aR <- Efit3$par[2]
 aS <- Efit3$par[3]
 rhoR*R/(dk+aR*R+aS*S)
}

Efun4 <- function(S,R) {
 rhoR <- Efit4$par[1]
 rhoS <- Efit4$par[2]
 aR <- Efit4$par[3]
 aS <- Efit4$par[4]
 (rhoR*R+rhoS*S)/(dk+aR*R+aS*S)
}

Eresults <- within(Eresults,Epred2 <- Efun2(S,R))
# Eresults <- within(Eresults,Epred3 <- Efun3(S,R))
# Eresults <- within(Eresults,Epred4 <- Efun4(S,R))
Egood <- subset(Eresults,Sample != "F8")

# Find estimate of sigma
sigma2 <- Efit2$value
x2 <- Efit2$par
## rhoR     rhoS      aS
## 96.7872 11.7068  0.3998
# rhoR is per cell production by R cells
# rhoS is per cell production by S cells
# aS is per cell uptake by S cells relative to R cells
xc <- x2
chkfun <- function(xc) nrow(Egood)*(Efitfun2(xc)/(2*sigma2)-0.5)
chkfun(xc)

# Find confidence limits.  cdiff is the critical log likelihood
# difference, here set to the traditional 2
xc <- x2
cdiff <- 2
f1 <- function(x) {
       xc[1] <- x
       chkfun(xc)-cdiff
}
lims1 <- c(uniroot(f1,c(x2[1],0.6*x2[1]))$root,
           uniroot(f1,c(x2[1],1.4*x2[1]))$root)
# 91.492 102.080

xc <- x2
f2 <- function(x) {
       xc[2] <- x
       chkfun(xc)-cdiff
}
lims2 <- c(uniroot(f2,c(x2[2],0.6*x2[2]))$root,
           uniroot(f2,c(x2[2],1.4*x2[2]))$root)
#  9.8781 13.5364

xc <- x2
f3 <- function(x) {
       xc[3] <- x
       chkfun(xc)-cdiff
}
lims3 <- c(uniroot(f3,c(x2[3],0.6*x2[3]))$root,
           uniroot(f3,c(x2[3],1.4*x2[3]))$root)
# 0.35470 0.46016


