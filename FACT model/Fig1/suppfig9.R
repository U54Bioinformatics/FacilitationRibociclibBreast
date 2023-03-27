########################################################################
# Illustrate facilitation with other treatments (allfacil.R)
# Letrozole1 with lvl2=2: Point, no facilitation of S cells

pdf("suppfig9.pdf")
par(mar=c(5,5,4,2),font=2, font.axis=2, font.lab=2,
#    font.main=2,family="Arial")
    font.main=2)
days <- c(4,7,11,14,18)
lvl2s <- c(2,2,3)

########################################################################
# Illustrate facilitation modification (allfacilmod.R)
# Letrozole1 with lvl2=2: Point, facilitation of S cells reduced

inm <- 1
  nm <- "Letrozole4"
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  lvl1 <- 2
  lvl2  <- 2

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
                  alphaRS=competeR.rK$alpha,
                  alphaSR=competeS.rK$alpha)

  facilparms <- c(rS=facilS.rK$r,
                KS=facilS.rK$K,
                rR=facilR.rK$r,
                KR=facilR.rK$K,
                alphaRS=facilR.rK$alpha,
                alphaSR=facilS.rK$alpha)

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

# Plot up the data for S
  ymax <- max(c(controlS$cells,riboS$cells,competeS$cells,facilS$cells))
  plot(cells ~ day,controlS,pch="",col=camacols["S"],
       xlim=c(4,18.5),ylim=c(0,ymax),
       xlab="Day",ylab="Cell number",cex.lab=1.7,cex.axis=1.3)
  points(cells ~ day,controlS,pch=1,col=camacols["S"])
  points(cells ~ day,riboS,pch=19,col=camacols["S"])
  points(cells ~ day,competeS,pch=1,col="purple")
  points(cells ~ day,facilS,pch=19,col="purple")
  legend("topleft",c("Monoculture: Letrozole",
                     "Monoculture: Ribociclib and letrozole",
                     "Coculture: Letrozole",
                     "Coculture: Ribociclib and letrozole",
                     "Expected combined effect"),lwd=3,seg.len=4,
         col=c(camacols["S"],camacols["S"],"purple","purple","black"),
         lty=c(2,1,2,1,1),cex=1.0)

# Add predicted to the graph
  lines(cells ~ day,untreatedpredS,lwd=2,lty=1,col="gray")
  lines(cells ~ day,controlpredS,lwd=3,lty=2,col=camacols["S"])
  lines(cells ~ day,ribopredS,lwd=3,col=camacols["S"])
  lines(days,competepredS,lwd=3,lty=2,col="purple")
  lines(days,facilpredS,lwd=3,col="purple")
  lines(days,nullpredS,lwd=3,col="black")

# Add arrow
    lday <- 5
    fudge <- -0.2
    fp <- facilpredS[lday]
    lp <- nullpredS[lday]
    dy <- max(days)+0.7
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple4",code=1,lty=1,lwd=3,length=0.1)

# Do the stats
facilS$nullS <- nullpredS
facilS$logratS <- with(facilS,log(cells/nullS))
# print(summary(lm(logratS ~ day,facilS)))
## (Intercept)  -0.4352     0.1913   -2.28  0.04045 *  
## day           0.0720     0.0161    4.47  0.00063 ***
# There is still facilitation, but this has to be compared with the
# facilitation in the absence of a modifier.
letrofacilS <- facilS
# Just remake facilS without Letrozole
source("noestra.R")
noestrafacilS <- facilS

letrofacilS$flag <- "E"
noestrafacilS$flag <- "N"
allletro <- rbind(letrofacilS,noestrafacilS)
# print(summary(lm(logratS ~ day*flag,allletro)))
## (Intercept)  -0.4352     0.1709   -2.55  0.01716 *  
## day           0.0720     0.0144    5.00  3.3e-05 ***
## flagN         0.0275     0.2417    0.11  0.91043    
## day:flagN     0.0759     0.0203    3.73  0.00093 ***

########################################################################
# FOHTamoxifen with level2=1: Point, facilitation of S cells cancelled

inm <- 1
  nm <- "FOHTamoxifen"
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))

# Pull out the S files
  controlS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                   level1==0 & level2==1)
  riboS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                  level1==200 & level2==1)
  competeS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                   level1==0 & level2==1)
  facilS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                  level1==200 & level2==1)

# Pull out the R files
  controlR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                   level1==0 & level2==1)
  riboR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                  level1==200 & level2==1)
  competeR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                   level1==0 & level2==1)
  facilR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                  level1==200 & level2==1)

# Pull out the fits
  controlS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                   level1==0 & level2==1)
  riboS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                  level1==200 & level2==1)
  competeS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                   level1==0 & level2==1)
  facilS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                  level1==200 & level2==1)

  controlR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                   level1==0 & level2==1)
  riboR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                  level1==200 & level2==1)
  competeR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                   level1==0 & level2==1)
  facilR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                  level1==200 & level2==1)

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
                  alphaRS=competeR.rK$alpha,
                  alphaSR=competeS.rK$alpha)

  facilparms <- c(rS=facilS.rK$r,
                KS=facilS.rK$K,
                rR=facilR.rK$r,
                KR=facilR.rK$K,
                alphaRS=facilR.rK$alpha,
                alphaSR=facilS.rK$alpha)

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

# Plot up the data for S
  ymax <- max(c(controlS$cells,riboS$cells,competeS$cells,facilS$cells))
  plot(cells ~ day,controlS,pch="",col=camacols["S"],
       xlim=c(4,18.5),ylim=c(0,ymax),
       xlab="Day",ylab="Cell number",cex.lab=1.7,cex.axis=1.3)
  points(cells ~ day,controlS,pch=1,col=camacols["S"])
  points(cells ~ day,riboS,pch=19,col=camacols["S"])
  points(cells ~ day,competeS,pch=1,col="purple")
  points(cells ~ day,facilS,pch=19,col="purple")
  legend("topleft",c("Monoculture: Tamoxifen",
                     "Monoculture: Ribociclib and tamoxifen",
                     "Coculture: Tamoxifen",
                     "Coculture: Ribociclib and tamoxifen",
                     "Expected combined effect"),lwd=3,seg.len=4,
         col=c(camacols["S"],camacols["S"],"purple","purple","black"),
         lty=c(2,1,2,1,1),cex=1.0)

# Add predicted to the graph
  lines(cells ~ day,untreatedpredS,lwd=2,lty=1,col="gray")
  lines(cells ~ day,controlpredS,lwd=3,lty=2,col=camacols["S"])
  lines(cells ~ day,ribopredS,lwd=3,col=camacols["S"])
  lines(days,competepredS,lwd=3,lty=2,col="purple")
  lines(days,facilpredS,lwd=3,col="purple")
  lines(days,nullpredS,lwd=3,col="black")

# Add arrow
    lday <- 5
    fudge <- -0.2
    fp <- facilpredS[lday]
    lp <- nullpredS[lday]
    dy <- max(days)+0.7
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple4",code=1,lty=1,lwd=3,length=0.1)

# Do the stats
facilS$nullS <- nullpredS
facilS$logratS <- with(facilS,log(cells/nullS))
# print(summary(lm(logratS ~ day,facilS)))
## (Intercept)  0.32204    0.10612    3.03   0.0096 **
## day         -0.02854    0.00893   -3.20   0.0070 **
# There is negative facilitation, but this has to be compared with the
# facilitation in the absence of a modifier.
tamoxfacilS <- facilS
# Just remake facilS without Tamoxifen
tamoxfacilS$flag <- "E"
tmp <- noestrafacilS
# tmp$ilvl1 <- NULL
# tmp$ilvl2 <- NULL
alltamox <- rbind(tamoxfacilS,tmp)
# print(summary(lm(logratS ~ day*flag,alltamox)))
## (Intercept)   0.3220     0.1287    2.50  0.01895 *  
## day          -0.0285     0.0108   -2.64  0.01398 *  
## flagN        -0.7298     0.1820   -4.01  0.00046 ***
## day:flagN     0.1765     0.0153   11.52    1e-11 ***

########################################################################
# Raloxifene with level2=75: Point, facilitation of S cells cancelled

inm <- 1
  nm <- "Raloxifene3"
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))

# Pull out the S files
  controlS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                   level1==0 & level2==75)
  riboS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                  level1==200 & level2==75)
  competeS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                   level1==0 & level2==75)
  facilS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                  level1==200 & level2==75)

# Pull out the R files
  controlR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                   level1==0 & level2==75)
  riboR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                  level1==200 & level2==75)
  competeR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                   level1==0 & level2==75)
  facilR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                  level1==200 & level2==75)

# Pull out the fits
  controlS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                   level1==0 & level2==75)
  riboS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                  level1==200 & level2==75)
  competeS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                   level1==0 & level2==75)
  facilS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                  level1==200 & level2==75)

  controlR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                   level1==0 & level2==75)
  riboR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                  level1==200 & level2==75)
  competeR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                   level1==0 & level2==75)
  facilR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                  level1==200 & level2==75)

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
                  alphaRS=competeR.rK$alpha,
                  alphaSR=competeS.rK$alpha)

  facilparms <- c(rS=facilS.rK$r,
                KS=facilS.rK$K,
                rR=facilR.rK$r,
                KR=facilR.rK$K,
                alphaRS=facilR.rK$alpha,
                alphaSR=facilS.rK$alpha)

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

# Plot up the data for S
  ymax <- max(c(controlS$cells,riboS$cells,competeS$cells,facilS$cells))
  plot(cells ~ day,controlS,pch="",col=camacols["S"],
       xlim=c(4,18.5),ylim=c(0,ymax),
       xlab="Day",ylab="Cell number",cex.lab=1.7,cex.axis=1.3)
  points(cells ~ day,controlS,pch=1,col=camacols["S"])
  points(cells ~ day,riboS,pch=19,col=camacols["S"])
  points(cells ~ day,competeS,pch=1,col="purple")
  points(cells ~ day,facilS,pch=19,col="purple")
  legend("topleft",c("Monoculture: Raloxifene",
                     "Monoculture: Ribociclib and raloxifene",
                     "Coculture: Raloxifene",
                     "Coculture: Ribociclib and raloxifene",
                     "Expected combined effect"),lwd=3,seg.len=4,
         col=c(camacols["S"],camacols["S"],"purple","purple","black"),
         lty=c(2,1,2,1,1),cex=1.0)

# Add predicted to the graph
  lines(cells ~ day,untreatedpredS,lwd=2,lty=1,col="gray")
  lines(cells ~ day,controlpredS,lwd=3,lty=2,col=camacols["S"])
  lines(cells ~ day,ribopredS,lwd=3,col=camacols["S"])
  lines(days,competepredS,lwd=3,lty=2,col="purple")
  lines(days,facilpredS,lwd=3,col="purple")
  lines(days,nullpredS,lwd=3,col="black")

# Add arrow
    lday <- 5
    fudge <- -0.2
    fp <- facilpredS[lday]
    lp <- nullpredS[lday]
    dy <- max(days)+0.7
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple4",code=1,lty=1,lwd=3,length=0.1)

# Do the stats
facilS$nullS <- nullpredS
facilS$logratS <- with(facilS,log(cells/nullS))
# print(summary(lm(logratS ~ day,facilS)))
## (Intercept) -0.01810    0.06005   -0.30     0.77    
## day          0.03065    0.00505    6.07    4e-05 ***
# There is still facilitation, but this has to be compared with the
# facilitation in the absence of a modifier.
raloxfacilS <- facilS
# Just remake facilS without Raloxifene
raloxfacilS$flag <- "E"
tmp <- noestrafacilS
# tmp$ilvl1 <- NULL
# tmp$ilvl2 <- NULL
allralox <- rbind(raloxfacilS,tmp)
# print(summary(lm(logratS ~ day*flag,allralox)))
## (Intercept)  -0.0181     0.1128   -0.16   0.8738    
## day           0.0307     0.0095    3.23   0.0034 ** 
## flagN        -0.3896     0.1596   -2.44   0.0217 *  
## day:flagN     0.1173     0.0134    8.73  3.3e-09 ***

dev.off()
