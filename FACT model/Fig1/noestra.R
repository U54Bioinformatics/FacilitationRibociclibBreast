inm <- 1
  nm <- camafiles[inm]
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  lvl1 <- 2
  lvl2  <- 1

# Pull out the S files
  controlS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                   ilvl1==1 & ilvl2==lvl2)
  riboS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                  ilvl1==lvl1 & ilvl2==lvl2)
  competeS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                   ilvl1==1 & ilvl2==lvl2)
  facilS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                  ilvl1==lvl1 & ilvl2==lvl2)


# Pull out the R files
  controlR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                   ilvl1==1 & ilvl2==lvl2)
  riboR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                  ilvl1==lvl1 & ilvl2==lvl2)
  competeR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                   ilvl1==1 & ilvl2==lvl2)
  facilR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                  ilvl1==lvl1 & ilvl2==lvl2)

# Pull out the fits
  controlS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                   ilvl1==1 & ilvl2==lvl2)
  riboS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                  ilvl1==lvl1 & ilvl2==lvl2)
  competeS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                   ilvl1==1 & ilvl2==lvl2)
  facilS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                  ilvl1==lvl1 & ilvl2==lvl2)

  controlR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                   ilvl1==1 & ilvl2==lvl2)
  riboR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                  ilvl1==lvl1 & ilvl2==lvl2)
  competeR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                   ilvl1==1 & ilvl2==lvl2)
  facilR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                  ilvl1==lvl1 & ilvl2==lvl2)

# Make the predictions
  controlpredS <- data.frame(day=days,
                cells=lfun(mean(controlS$cells[controlS$day==4]),
                                controlS.rK$K,controlS.rK$r,days))

  ribopredS <- data.frame(day=days,
                cells=lfun(mean(riboS$cells[riboS$day==4]),
                                riboS.rK$K,riboS.rK$r,days))

  controlpredR <- data.frame(day=days,
                cells=lfun(mean(controlR$cells[controlR$day==4]),
                                controlR.rK$K,controlR.rK$r,days))

  ribopredR <- data.frame(day=days,
                cells=lfun(mean(riboR$cells[riboR$day==4]),
                                riboR.rK$K,riboR.rK$r,days))

# Pull out the parameters 
  competeparms <- c(rS=competeS.rK$r,
                  KS=competeS.rK$K,
                  rR=competeR.rK$r,
                  KR=competeR.rK$K,
                  alphaRS=competeR.rK$alphaval,
                  alphaSR=competeS.rK$alphaval)

  facilparms <- c(rS=facilS.rK$r,
                KS=facilS.rK$K,
                rR=facilR.rK$r,
                KR=facilR.rK$K,
                alphaRS=facilR.rK$alphaval,
                alphaSR=facilS.rK$alphaval)

  init <- c(mean(competeS$cells[competeS$day==4]),
          mean(competeR$cells[competeR$day==4]))
  names(init) <- c("S","R")
  LVout <- as.data.frame(ode(y=init,parms=competeparms,times=days,func=LVmodel))
  competepredS <- LVout$S
  competepredR <- LVout$R

  init <- c(mean(facilS$cells[facilS$day==4]),mean(facilR$cells[facilR$day==4]))
  names(init) <- c("S","R")
  LVout <- as.data.frame(ode(y=init, parms=facilparms,times=days, func=LVmodel))
  facilpredS <- LVout$S
  facilpredR <- LVout$R

# Make the null models to compare
  nullparms <- facilparms
  nullparms["rS"] <- (riboS.rK$r)*(competeS.rK$r)/controlS.rK$r
  nullparms["KS"] <- (riboS.rK$K)*(competeS.rK$K)/controlS.rK$K
  nullparms["rR"] <- (riboR.rK$r)*(competeR.rK$r)/controlR.rK$r
  nullparms["KR"] <- (riboR.rK$K)*(competeR.rK$K)/controlR.rK$K
  LVout <- as.data.frame(ode(y=init, parms=nullparms,times=days, func=LVmodel))
  nullpredS <- LVout$S
  nullpredR <- LVout$R

facilS$nullS <- nullpredS
facilS$logratS <- with(facilS,log(cells/nullS))

