# Functions for constrained logistic fits
rmax <- 0.5
Kmax <- 3e5

lfit <- function(x) {
 K <- x[1]
 r <- x[2]
 firstc <- mean(tmp$cells[tmp$day==4])
 sum((K*firstc/((K-firstc)*exp(-r*(tmp$day-4))+firstc)-tmp$cells)^2)*
       (1+abs(r)*(abs(r) > rmax))*  #Pay a price for extreme values
       (1+log(abs(K))*(abs(K) > Kmax))  #Pay a price for extreme values
}

lfun <- function(firstc,K,r,day) {
        K*firstc/((K-firstc)*exp(-r*(day-4))+firstc)
      }

LVmodel <- function(t, x, parms) {
# Extract state variable from x
  S <- x["S"]
  R <- x["R"]
  chkS <- with(as.list(parms),(1-(S+alphaSR*R)/KS))
  dS <- with(as.list(parms),rS*S*chkS*(1-2*(rS<0 & chkS<0)))
  chkR <- with(as.list(parms),(1-(alphaRS*S+R)/KR))
  dR <- with(as.list(parms),rR*R*chkR*(1-2*(rR<0 & chkR<0)))
  res<- c(dS,dR)
  names(res) <- paste0("d",names(x))
  return(list(res))
}

# Estimate LV coefficients only, use tmpSM, tmpRM and parms
alphamax <- 10
LVfunalpha <- function(x) {
  parms["alphaRS"] <- x[1]
  parms["alphaSR"] <- x[2]
  LVout <- as.data.frame(ode(y=init, parms=parms,
                              times=days, func=LVmodel))
  tmpSM$pred <- LVout$S
  tmpRM$pred <- LVout$R
  RMS <- with(tmpSM,sqrt(mean((cells-pred)^2)))+
              with(tmpRM,sqrt(mean((cells-pred)^2)))*
#              (1+with(as.list(parms),abs(alphaRS)*(abs(alphaRS) > alphamax))+
#                 with(as.list(parms),abs(alphaSR)*(abs(alphaSR) > alphamax)))
              (1+with(as.list(parms),alphamax*(alphaRS < 0))+
                 with(as.list(parms),alphamax*(alphaSR < 0))+
                 with(as.list(parms),abs(alphaRS)*(abs(alphaRS) > alphamax))+
                 with(as.list(parms),abs(alphaSR)*(abs(alphaSR) > alphamax)))
  return(RMS)
}

# Estimate r and K only, use tmpSM, tmpRM and parms, and don't let the
# values go too crazy
rmin <- -0.1
Kmin <- 1e3
LVfunrK <- function(x) {
  parms["rS"] <- x[1]
  parms["KS"] <- x[2]
  parms["rR"] <- x[3]
  parms["KR"] <- x[4]
  LVout <- as.data.frame(ode(y=init, parms=parms,
                              times=days, func=LVmodel))
  tmpSM$pred <- LVout$S
  tmpRM$pred <- LVout$R
  RMS <- (with(tmpSM,sqrt(mean((cells-pred)^2)))+
              with(tmpRM,sqrt(mean((cells-pred)^2))))*
              (1+with(as.list(parms),abs(rS)*(rS < rmin | rS > rmax))+
                 with(as.list(parms),abs(rR)*(rR < rmin | rR > rmax))+
                 with(as.list(parms),abs(KS)*(KS < Kmin | KS > Kmax))+
                 with(as.list(parms),abs(KR)*(KR < Kmin | KR > Kmax)))
  return(RMS)
}


# Alternative models with facilitation
facilmodel <- function(t, x, fparms) {
  S <- x["S"]
  R <- x["R"]
  g <- with(as.list(fparms),(rS+beta*rS*R/kE)/(1+R/kE))
  chkS <- with(as.list(fparms),(1-(S+R)/KS))
  dS <- with(as.list(fparms),g*S*chkS*(1-2*(rS<0 & chkS<0)))
  chkR <- with(as.list(fparms),(1-(S+R)/KR))
  dR <- with(as.list(fparms),rR*R*chkR*(1-2*(rR<0 & chkR<0)))
  res<- c(dS,dR)
  names(res) <- paste0("d",names(x))
  return(list(res))
}

betamax <- 5.0
kEmax <- 2*Kmax
kEmin <- 2e3
facilfunmax <- function(x) {
  fparms["beta"] <- x[1]
  fparms["kE"] <- x[2]
  facilout <- as.data.frame(ode(y=init, parms=fparms,
                              times=days, func=facilmodel))
  tmpSM$pred <- facilout$S
  tmpRM$pred <- facilout$R
  RMS <- with(tmpSM,sqrt(mean((cells-pred)^2)))+
              with(tmpRM,sqrt(mean((cells-pred)^2)))*
              (1+with(as.list(fparms),betamax*(beta < 0))+
                 with(as.list(fparms),betamax*(kE < 0))+
                 with(as.list(fparms),abs(beta)*(abs(beta) > betamax))+
                 with(as.list(fparms),abs(beta)*(abs(kE) > kEmax))+
                 with(as.list(fparms),abs(beta)*(abs(kE) < kEmin)))
 return(RMS)
}

