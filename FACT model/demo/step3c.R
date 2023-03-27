# Estimate the alphas
alphaRSval <- 1.4  #Consensus values
alphaSRval <- 0.3
# alphaRSval <- 1.88  #Consensus values for Letrozole
# alphaSRval <- 0.28
tparms <- NULL
for (nm in usefiles) {
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  tmpRM.rK <- subset(tfile.rK,type=="R" & env=="M" & level1==0 & level2==0)
  tmpSM.rK <- subset(tfile.rK,type=="S" & env=="M" & level1==0 & level2==0)

# Pull out the data from the M files
  tmpRM0 <- subset(tfile,type =="R" & env=="M" & level1==0 & level2==0)
  tmpSM0 <- subset(tfile,type =="S" & env=="M" & level1==0 & level2==0)
  tmpRM <- subset(tmpRM0,day > 0)
  tmpSM <- subset(tmpSM0,day > 0)
  days <- unique(tmpRM$day)
 
# Pull out the parameters from the S files
  parms <- c(rS=tmpSM.rK$r,
             KS=tmpSM.rK$K,
             rR=tmpRM.rK$r,
             KR=tmpRM.rK$K,
             alphaRS=tmpSM.rK$alpha,
             alphaSR=tmpRM.rK$alpha)

# Estimate alphas
  init <- c(mean(tmpSM$cells[tmpSM$day==4]),mean(tmpRM$cells[tmpRM$day==4]))
  names(init) <- c("S","R")
  LV.opt <- optim(c(alphaRSval,alphaSRval),LVfunalpha)
  parms["alphaRS"] <- LV.opt$par[1]
  parms["alphaSR"] <- LV.opt$par[2]
  print(parms)
  tparms <- rbind(tparms,parms)

  if (plotm==1) {
    LVout <- as.data.frame(ode(y=init, parms=parms,
                               times=days, func=LVmodel))
    tmpSM$pred <- LVout$S
    tmpRM$pred <- LVout$R

    testparms <- parms
    testparms["alphaRS"] <- alphaRSval
    testparms["alphaSR"] <- alphaSRval
    LVout <- as.data.frame(ode(y=init, parms=testparms,
                               times=days, func=LVmodel))
    tmpSM$testpred <- LVout$S
    tmpRM$testpred <- LVout$R

    plot(cells ~ day,tmpSM0,pch=19,col=camacols["S"])
    points(cells ~ day,tmpRM0,pch=19,col=camacols["R"])
    lines(pred ~ day,tmpSM,subset=rep==1,pch=19,lwd=3,col=camacols["S"])
    lines(pred ~ day,tmpRM,subset=rep==1,pch=19,lwd=3,col=camacols["R"])
    lines(testpred ~ day,tmpSM,subset=rep==1,pch=19,lwd=3,lty=2,col=camacols["S"])
    lines(testpred ~ day,tmpRM,subset=rep==1,pch=19,lwd=3,lty=2,col=camacols["R"])
    title(main=nm)
    readline('hit return for next plot> ')
#      Sys.sleep(0.1)
  }
  tfile.rK$alpha[tfile.rK$env=="M" & tfile.rK$type=="R"] <- parms["alphaRS"]
  tfile.rK$alpha[tfile.rK$env=="M" & tfile.rK$type=="S"] <- parms["alphaSR"]
  tfile.rK$alphaval[tfile.rK$env=="M" & tfile.rK$type=="R"] <- alphaRSval
  tfile.rK$alphaval[tfile.rK$env=="M" & tfile.rK$type=="S"] <- alphaSRval

# Assign to its name and just overwrite cama.rK
  assign(paste0(nm,".rK"),tfile.rK)
}

