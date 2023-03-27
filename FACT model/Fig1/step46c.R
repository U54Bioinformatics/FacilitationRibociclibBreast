# For steps 4 and 6, use the estimated alphas to estimate rR and KR 
all.rK <- NULL
# for (nm in camafiles) {
for (nm in usefiles) {
  print(nm)
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  tmpSM.rK <- subset(tfile.rK,type=="S" & env=="M" & (level1>0 | level2>0))
  tmpRM.rK <- subset(tfile.rK,type=="R" & env=="M" & (level1>0 | level2>0))
  tmpSS.rK <- subset(tfile.rK,type=="S" & env=="S" & (level1>0 | level2>0))
  tmpRS.rK <- subset(tfile.rK,type=="R" & env=="R" & (level1>0 | level2>0))
  for (irow in 1:nrow(tmpRM.rK)) {
# Pull out the data from the M files
    tmpSM0 <- subset(tfile,type =="S" & env=="M" & 
                    level1==tmpRM.rK$level1[irow] & level2==tmpRM.rK$level2[irow])
    tmpRM0 <- subset(tfile,type =="R" & env=="M" & 
                    level1==tmpRM.rK$level1[irow] & level2==tmpRM.rK$level2[irow])
    tmpSM <- subset(tmpSM0,day > 0)
    tmpRM <- subset(tmpRM0,day > 0)
    days <- unique(tmpRM$day)
 
# Pull out the parameters 
    parms <- c(rS=tmpSM.rK$r[irow],
               KS=tmpRM.rK$K[irow],
               rR=tmpRM.rK$r[irow],
               KR=tmpRM.rK$K[irow],
               alphaRS=tmpRM.rK$alpha[irow],
               alphaSR=tmpSM.rK$alpha[irow])

#  Estimate r and K given alpha with LVfunrK
    init <- c(mean(tmpSM$cells[tmpSM$day==4]),mean(tmpRM$cells[tmpRM$day==4]))
    names(init) <- c("S","R")
    guess <- c(0.2,2e5,0.2,2e5)
    if (nm==camafiles[2] & irow %in% 2:4) guess <- c(4.0719e-01,2.0459e+05,3.1106e-01,2.1971e+05)
    LV.opt <- optim(guess,LVfunrK)
    parms["rS"] <- LV.opt$par[1]
    parms["KS"] <- LV.opt$par[2]
    parms["rR"] <- LV.opt$par[3]
    parms["KR"] <- LV.opt$par[4]
    
    print(c(parms,level1=tmpRM.rK$level1[irow],level2=tmpRM.rK$level2[irow]))
    tmpSM.rK$r[irow] <- parms["rS"]
    tmpSM.rK$K[irow] <- parms["KS"]
    tmpRM.rK$r[irow] <- parms["rR"]
    tmpRM.rK$K[irow] <- parms["KR"]

    if (plotm==1) {
      LVout <- as.data.frame(ode(y=init, parms=parms,
                                 times=days, func=LVmodel))
      tmpSM$pred <- LVout$S
      tmpRM$pred <- LVout$R

      plot(cells ~ day,tmpSM0,ylim=range(c(tmpSM0$cells,tmpRM0$cells)),pch=19,col="green")
      points(cells ~ day,tmpRM0,pch=19,col="cyan")
      lines(pred ~ day,tmpSM,subset=rep==1,pch=19,lwd=3,col="green")
      lines(pred ~ day,tmpRM,subset=rep==1,pch=19,lwd=3,col="cyan")
      title(main=nm)
      readline('hit return for next plot> ')
#      Sys.sleep(0.1)
    }
  }
# Replace values in tfile with those in tmpRM for suitable rows
  tfile.rK$r[tfile.rK$env=="M" & 
             tfile.rK$type=="S"& 
             (tfile.rK$level1>0 | tfile.rK$level2>0)] <- tmpSM.rK$r
  tfile.rK$K[tfile.rK$env=="M" & 
             tfile.rK$type=="S"& 
             (tfile.rK$level1>0 | tfile.rK$level2>0)] <- tmpSM.rK$K
  tfile.rK$r[tfile.rK$env=="M" & 
             tfile.rK$type=="R"& 
             (tfile.rK$level1>0 | tfile.rK$level2>0)] <- tmpRM.rK$r
  tfile.rK$K[tfile.rK$env=="M" & 
             tfile.rK$type=="R"& 
             (tfile.rK$level1>0 | tfile.rK$level2>0)] <- tmpRM.rK$K
  assign(paste0(nm,".rK"),tfile.rK)
  all.rK <- rbind(all.rK,tfile.rK)
} 

# Save it in case you stupidly start the program again 
# cama.rK <- all.rK
# all.rK.save <- all.rK

