# Step three for LY2 cells
# pdf('Figs/step3l.pdf')
tparms <- NULL
nm <- "LY2"
Mnms <- paste0("M",1:2) 
# Consensus values
alphaRSval <- 0.2
alphaSRval <- 1.5
for (envm in Mnms) {
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  tmpRS.rK <- subset(tfile.rK,type=="R" & env=="R0" & level1==0 & level2==0)
  tmpSS.rK <- subset(tfile.rK,type=="S" & env=="S0" & level1==0 & level2==0)

# Pull out the data from the M files
  tmpRM0 <- subset(tfile,type =="R" & env==envm & level1==0 & level2==0)
  tmpSM0 <- subset(tfile,type =="S" & env==envm & level1==0 & level2==0)
  tmpRM <- subset(tmpRM0,day > 0)
  tmpSM <- subset(tmpSM0,day > 0)
  days <- unique(tmpRM$day)

# Pull out the parameters from the S files
  parms <- c(rS=tmpSS.rK$r,
             KS=tmpSS.rK$K,
             rR=tmpRS.rK$r,
             KR=tmpRS.rK$K,
             alphaRS=tmpSS.rK$alpha,
             alphaSR=tmpRS.rK$alpha)

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

    plot(cells ~ day,tmpSM0,pch=19,ylim=range(tmpSM0$cells,tmpRM0$cells),
         col=LY2cols["S"],xlab="Day",ylab="Cells",cex.lab=1.3,cex.axis=1.3)
    points(cells ~ day,tmpRM0,pch=19,col=LY2cols["R"])
    lines(pred ~ day,tmpSM,subset=rep==1,pch=19,lwd=3,col=LY2cols["S"])
    lines(pred ~ day,tmpRM,subset=rep==1,pch=19,lwd=3,col=LY2cols["R"])
    lines(testpred ~ day,tmpSM,subset=rep==1,pch=19,lwd=3,lty=2,col=LY2cols["S"])
    lines(testpred ~ day,tmpRM,subset=rep==1,pch=19,lwd=3,lty=2,col=LY2cols["R"])
    title(main=envm,cex=1.3)
    readline('hit return for next plot> ')
  }
  tfile.rK$r[tfile.rK$env==envm & tfile.rK$type=="S"] <- parms["rS"]
  tfile.rK$K[tfile.rK$env==envm & tfile.rK$type=="S"] <- parms["KS"]
  tfile.rK$r[tfile.rK$env==envm & tfile.rK$type=="R"] <- parms["rR"]
  tfile.rK$K[tfile.rK$env==envm & tfile.rK$type=="R"] <- parms["KR"]
  tfile.rK$alpha[tfile.rK$env==envm & tfile.rK$type=="R"] <- parms["alphaRS"]
  tfile.rK$alpha[tfile.rK$env==envm & tfile.rK$type=="S"] <- parms["alphaSR"]

  tfile.rK$alphaval[tfile.rK$env==envm & tfile.rK$type=="R"] <- alphaRSval
  tfile.rK$alphaval[tfile.rK$env==envm & tfile.rK$type=="S"] <- alphaSRval

# Assign to its name and just overwrite LY2.rK
  assign(paste0(nm,".rK"),tfile.rK)
  LY2.rK[LY2.rK$file==nm,] <- tfile.rK
}
# dev.off()

tparms <- as.data.frame(tparms)
row.names(tparms) <- Mnms 
#pdf("Figs/alphas.pdf")
# plot(alphaSR ~ alphaRS,tparms,pch=19,col="blue")
# abline(0,1,col="gray",lty=2)
# abline(v=alphaRSval,col="red",lty=2)
# abline(h=alphaSRval,col="red",lty=2)
# dev.off()

