########################################################################
# Illustrate facilitation with other treatments (allfacil.R)
# Estradiol1 with lvl2=2: Point, no facilitation of S cells

pdf("suppfig8.pdf")

#########################################################
############ Supp Fig 8a, the pedagogical figure ########
#########################################################

par(mar=c(5.1,5.1,8.1,2.1),xpd=TRUE)
nm <- camafiles[2]
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  cases <- subset(tfile.rK,level1 > 0 & level2 == 0 & env=="M" & type=="S")
  controlS <- subset(tfile,env=="S" & type=="S" & day > 0 & day < 21 &
                    level1==0 & level2==0)
  riboS <- subset(tfile,env=="S" & type=="S" & day > 0 & day < 21 &
                    level1==cases$level1 & level2==0)
  competeS <- subset(tfile,env=="M" & type=="S" & day > 0 & day < 21 &
                    level1==0 & level2==0)
  facilS <- subset(tfile,env=="M" & type=="S" & day > 0 & day < 21 &
                    level1==cases$level1 & level2==0)
  days <- unique(controlS$day)

  controlR <- subset(tfile,env=="R" & type=="R" & day > 0 & day < 21 &
                    level1==0 & level2==0)
  riboR <- subset(tfile,env=="R" & type=="R" & day > 0 & day < 21 &
                    level1==cases$level1 & level2==0)
  competeR <- subset(tfile,env=="M" & type=="R" & day > 0 & day < 21 &
                    level1==0 & level2==0)
  facilR <- subset(tfile,env=="M" & type=="R" & day > 0 & day < 21 &
                    level1==cases$level1 & level2==0)

# Pull out the fits
  controlS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                    level1==0 & level2==0)
  riboS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                    level1==cases$level1 & level2==0)
  competeS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                    level1==0 & level2==0)
  facilS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                    level1==cases$level1 & level2==0)

  controlR.rK <- subset(tfile.rK,env=="R" & type=="R" &
                    level1==0 & level2==0)
  riboR.rK <- subset(tfile.rK,env=="R" & type=="R" &
                    level1==cases$level1 & level2==0)
  competeR.rK <- subset(tfile.rK,env=="M" & type=="R" &
                    level1==0 & level2==0)
  facilR.rK <- subset(tfile.rK,env=="M" & type=="R" &
                    level1==cases$level1 & level2==0)

# Make the predictions
  days <- c(4,7,11,14,18)
  controlpredS <- data.frame(day=days,
                  cells=lfun(mean(controlS$cells[controlS$day==4]),
                                  controlS.rK$K,controlS.rK$r,days))
# Save for later
  untreatedpredS <- controlpredS

  ribopredS <- data.frame(day=days,
                  cells=lfun(mean(riboS$cells[riboS$day==4]),
                                  riboS.rK$K,riboS.rK$r,days))

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

  init <- c(mean(facilS$cells[facilS$day==4]),mean(facilR$cells[facilR$day==4]))
  names(init) <- c("S","R")
  LVout <- as.data.frame(ode(y=init, parms=facilparms,times=days, func=LVmodel))
  facilpredS <- LVout$S

# Make the null models to compare
  nullparms <- facilparms
  nullparms["rS"] <- (riboS.rK$r)*(competeS.rK$r)/controlS.rK$r
  nullparms["KS"] <- (riboS.rK$K)*(competeS.rK$K)/controlS.rK$K
  nullparms["rR"] <- (riboR.rK$r)*(competeR.rK$r)/controlR.rK$r
  nullparms["KR"] <- (riboR.rK$K)*(competeR.rK$K)/controlR.rK$K
  LVout <- as.data.frame(ode(y=init, parms=nullparms,times=days, func=LVmodel))
  nullpredS <- LVout$S

# Plot up the data
    plot(cells ~ day,controlS,pch="",col=camacols["S"],
         xlim = c(4,max(days)+9.0),
         ylim = c(0,max(cells)*1.02),
         xlab="Day",ylab="Cell number",cex.lab=1.7,cex.axis=1.3)
    points(cells ~ day,controlS,pch=1,col=camacols["S"])
    points(cells ~ day,riboS,pch=19,col=camacols["S"])
    points(cells ~ day,competeS,pch=1,col="purple")
    points(cells ~ day,facilS,pch=19,col="purple")
    legend("topleft",inset=c(0,-0.3),
                     c("Untreated Monoculture",
                       "Treated Monoculture",
                       "Untreated Coculture",
                       "Treated Coculture",
                       "Expected: Treated Coculture"),lwd=1,seg.len=4,
           col=c(camacols["S"],camacols["S"],"purple","purple","black"),
           pch=c(1,19,1,19,NA),lty=c(2,1,2,1,1),cex=1.0)

    legend("topright",inset=c(0,-0.3),
                      c("a. Treatment cost",
                       "b. Competition cost",
                       "c. Observed combined cost",
                       "d. Expected combined cost",
                       "e. Facilitation"),
                        lwd=3,pch=rep(NA,5),
           col=c(camacols["S"],"purple","blue","black","purple4"),
           lty=rep(NA,5),seg.len=2)
  par(font = 5) #change font to get arrows
    legend("topright", inset=c(0.35,-0.3),legend = rep(NA,5), pch = rep(174,5),
       lwd = 3, col=c(camacols["S"],"purple","blue","black","purple4"),
       lty = rep(NA,5), bty = "n",seg.len=2) 
par(font = 1) #back to default
# Add predicted to the graph

    lines(cells ~ day,controlpredS,lwd=3,lty=2,col=camacols["S"])
    lines(cells ~ day,ribopredS,lwd=3,col=camacols["S"])
    lines(days,competepredS,lwd=3,lty=2,col="purple")
    lines(days,facilpredS,lwd=3,col="purple")
    lines(days,nullpredS,lwd=3,col="black")

# Here's the cost of competition arrow
    lday <- length(days)
    fudge <- 0.2
    fp <-  (controlpredS$cells[3]+controlpredS$cells[4])/2
    lp <- (competepredS[3]+competepredS[4])/2
    dy <- (days[3]+days[4])/2
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple",code=2,lty=1,lwd=2,length=0.1)
    tsize <- 1.2
    text(dy-0.5,(lp+fp)/2,"b. Competition cost",col="purple",pos=2,cex=tsize)

# Here's the cost of treatment arrow
    fudge <- 0.0
    fp <- controlpredS$cells[lday]
    lp <- ribopredS$cells[lday]
    dy <- max(days)+1
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col=camacols["S"],code=2,lty=1,lwd=3,length=0.1)
    text(dy+15*fudge,0.3*(lp+fp),"a. Treatment \n cost",
         col=camacols["S"],pos=2,cex=tsize)

# Here's the observed combined cost arrow
    fudge <- 0.0
    fp <- controlpredS$cells[lday]
    lp <- facilpredS[lday]
    dy <- max(days)+2
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="blue",code=2,lty=1,lwd=3,length=0.1)
    text(dy,(lp+fp)/2,"c. Observed \ncombined \ncost",
         col="blue",pos=4,cex=tsize)

# Here's the predicted combined cost arrow (has to jump)
    fudge <- 0.02
    dy <- max(days)+3
# Bottom bit
    fp <- 0.53*controlpredS$cells[lday]
    lp <- nullpredS[lday]
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col=1,code=2,lty=1,lwd=3,length=0.1)
# Top bit
    fudge <- 0.0
    fp <- controlpredS$cells[lday]
    lp <- 0.7*controlpredS$cells[lday]
    segments(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col=1,lty=1,lwd=3)

    fudge <- 0.02
    fp <- controlpredS$cells[lday]
    lp <- nullpredS[lday]
    text(dy,(lp+fp)/3,"d. Expected \ncombined \ncost",col=1,pos=4,cex=tsize)

# Here's the facilitation combined cost arrow
    fudge <- 0.02
    fp <- facilpredS[lday]
    lp <- nullpredS[lday]
    dy <- max(days)+4
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple4",code=3,lty=1,lwd=3,length=0.1)
    text(dy,(lp+fp)/2,"e. Facilitation",col="purple4",pos=4,cex=tsize)

par(mar=c(5.1,4.1,4.1,2.1))
# End of pedagogical figure

par(mar=c(5,5,4,2),font=2, font.axis=2, font.lab=2,
#    font.main=2,family="Arial")
    font.main=2)
days <- c(4,7,11,14,18)
lvl2s <- c(2,2,3)

xmax <- 19.0
dplus <- 0.8

# First, interaction of estradiol with competition
inm <- 1
  nm <- camafiles[inm]
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  lvl1 <- 1
  lvl2 <- lvl2s[inm]

# Pull out the S files
  controlS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                     ilvl1==1 & ilvl2==1)
  treatS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  competeS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                     ilvl1==1 & ilvl2==1)
  facilS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))

# Pull out the R files
  controlR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                     ilvl1==1 & ilvl2==1)
  treatR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  competeR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                     ilvl1==1 & ilvl2==1)
  facilR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))

# Pull out the fits
  controlS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                        ilvl1==1 & ilvl2==1)
  treatS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                      ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  competeS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                        ilvl1==1 & ilvl2==1)
  facilS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                      ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  
  controlR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                        ilvl1==1 & ilvl2==1)
  treatR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                      ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  competeR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                        ilvl1==1 & ilvl2==1)
  facilR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))

# Make the predictions
  controlpredS <- data.frame(day=days,
                  cells=lfun(mean(controlS$cells[controlS$day==4]),
                                  controlS.rK$K,controlS.rK$r,days))

  treatpredS <- data.frame(day=days,
                  cells=lfun(mean(treatS$cells[treatS$day==4]),
                                  treatS.rK$K,treatS.rK$r,days))

  controlpredR <- data.frame(day=days,
                  cells=lfun(mean(controlR$cells[controlR$day==4]),
                                  controlR.rK$K,controlR.rK$r,days))

  treatpredR <- data.frame(day=days,
                  cells=lfun(mean(treatR$cells[treatR$day==4]),
                                  treatR.rK$K,treatR.rK$r,days))

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
  nullparms["rS"] <- (treatS.rK$r)*(competeS.rK$r)/controlS.rK$r
  nullparms["KS"] <- (treatS.rK$K)*(competeS.rK$K)/controlS.rK$K
  nullparms["rR"] <- (treatR.rK$r)*(competeR.rK$r)/controlR.rK$r
  nullparms["KR"] <- (treatR.rK$K)*(competeR.rK$K)/controlR.rK$K
  LVout <- as.data.frame(ode(y=init, parms=nullparms,times=days, func=LVmodel))
  nullpredS <- LVout$S
  nullpredR <- LVout$R

# Plot up the data for S
  ymax <- max(c(controlS$cells,treatS$cells,competeS$cells,facilS$cells))
  plot(cells ~ day,controlS,pch="",col=camacols["S"],
       xlim=c(4,xmax),ylim=c(0,ymax),
       xlab="Day",ylab="Cell number",cex.lab=1.7,cex.axis=1.3)
  points(cells ~ day,controlS,pch=1,col=camacols["S"])
  points(cells ~ day,treatS,pch=19,col=camacols["S"])
  points(cells ~ day,competeS,pch=1,col="purple")
  points(cells ~ day,facilS,pch=19,col="purple")
  legend("topleft",c("Monoculture: Untreated",
                     "Monoculture: Estradiol",
                     "Coculture: Untreated",
                     "Coculture: Estradiol",
                     "Expected combined effect"),lwd=3,seg.len=4,
         col=c(camacols["S"],camacols["S"],"purple","purple","black"),
         lty=c(2,1,2,1,1),cex=1.0)

# Add predicted to the graph
  lines(cells ~ day,controlpredS,lwd=3,lty=2,col=camacols["S"])
  lines(cells ~ day,treatpredS,lwd=3,col=camacols["S"])
  lines(days,competepredS,lwd=3,lty=2,col="purple")
  lines(days,facilpredS,lwd=3,col="purple")
  lines(days,nullpredS,lwd=3,col="black")

# Add arrow
    lday <- 5
    fudge <- 0.02
    fp <- facilpredS[lday]
    lp <- nullpredS[lday]
    dy <- max(days)+dplus
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple4",code=1,lty=1,lwd=3,length=0.1)

# Do the stats
facilS$nullS <- nullpredS
facilS$logratS <- with(facilS,log(cells/nullS))
summary(lm(logratS ~ day,facilS))
## (Intercept)  0.06387    0.04053    1.58   0.1390   
## day         -0.01253    0.00341   -3.67   0.0028 **

########################################################################
# Illustrate synergy (allsynergy.R)
# Estradiol1 with lvl2=2: Point, E replaces the facilitation by R cells

nm  <- camafiles[1]
tfile <- get(nm)
tfile.rK <- get(paste0(nm,".rK"))
lvl1 <- 2
lvl2 <- 2

# Pull out the S files
controlS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                 ilvl1==1 & ilvl2==1)
riboS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                ilvl1==lvl1 & ilvl2==1)
treatS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                 ilvl1==1 & ilvl2==lvl2)
bothS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                ilvl1==lvl1 & ilvl2==lvl2)

# Pull out the R files
controlR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                 ilvl1==1 & ilvl2==1)
riboR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                ilvl1==lvl1 & ilvl2==1)
treatR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                 ilvl1==1 & ilvl2==lvl2)
bothR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                ilvl1==lvl1 & ilvl2==lvl2)


# Pull out the fits
controlS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                    ilvl1==1 & ilvl2==1)
riboS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                ilvl1==lvl1 & ilvl2==1)
treatS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                 ilvl1==1 & ilvl2==lvl2)
bothS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                 ilvl1==lvl1 & ilvl2==lvl2)

controlR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                    ilvl1==1 & ilvl2==1)
riboR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                ilvl1==lvl1 & ilvl2==1)
treatR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                 ilvl1==1 & ilvl2==lvl2)
bothR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                 ilvl1==lvl1 & ilvl2==lvl2)


# Make the predictions
controlpredS <- data.frame(day=days,
              cells=lfun(mean(controlS$cells[controlS$day==4]),
                              controlS.rK$K,controlS.rK$r,days))

ribopredS <- data.frame(day=days,
              cells=lfun(mean(riboS$cells[riboS$day==4]),
                              riboS.rK$K,riboS.rK$r,days))

treatpredS <- data.frame(day=days,
              cells=lfun(mean(treatS$cells[treatS$day==4]),
                              treatS.rK$K,treatS.rK$r,days))

bothpredS <- data.frame(day=days,
              cells=lfun(mean(bothS$cells[bothS$day==4]),
                              bothS.rK$K,bothS.rK$r,days))

controlpredR <- data.frame(day=days,
              cells=lfun(mean(controlR$cells[controlR$day==4]),
                              controlR.rK$K,controlR.rK$r,days))

ribopredR <- data.frame(day=days,
              cells=lfun(mean(riboR$cells[riboR$day==4]),
                              riboR.rK$K,riboR.rK$r,days))

treatpredR <- data.frame(day=days,
              cells=lfun(mean(treatR$cells[treatR$day==4]),
                              treatR.rK$K,treatR.rK$r,days))

bothpredR <- data.frame(day=days,
              cells=lfun(mean(bothR$cells[bothR$day==4]),
                              bothR.rK$K,bothR.rK$r,days))

# Make the null models to compare
nullparms <- facilparms
nullparms["rS"] <- (treatS.rK$r)*(riboS.rK$r)/controlS.rK$r
nullparms["KS"] <- (treatS.rK$K)*(riboS.rK$K)/controlS.rK$K
nullparms["rR"] <- (treatR.rK$r)*(riboR.rK$r)/controlR.rK$r
nullparms["KR"] <- (treatR.rK$K)*(riboR.rK$K)/controlR.rK$K

nullpredS <- data.frame(day=days,
              cells=lfun(mean(bothS$cells[bothS$day==4]),
                              nullparms["KS"],nullparms["rS"],days))

nullpredR <- data.frame(day=days,
              cells=lfun(mean(bothR$cells[bothR$day==4]),
                              nullparms["KR"],nullparms["rR"],days))

# Plot up the data for S
ymax <- max(c(controlS$cells,riboS$cells,treatS$cells,bothS$cells))
plot(cells ~ day,controlS,pch="",col=camacols["S"],
     xlim=c(4,xmax),ylim=c(0,ymax),
     xlab="Day",ylab="Cell number",cex.lab=1.7,cex.axis=1.3)
points(cells ~ day,controlS,pch=1,col=camacols["S"])
points(cells ~ day,riboS,pch=19,col=camacols["S"])
points(cells ~ day,treatS,pch=1,col="purple")
points(cells ~ day,bothS,pch=19,col="purple")
legend("topleft",c("Monoculture: Untreated",
                   "Monoculture: Ribociclib",
                   "Monoculture: Estradiol",
                   "Monoculture: Ribociclib and estradiol",
                   "Expected combined effect"),lwd=3,seg.len=4,
       col=c(camacols["S"],camacols["S"],"purple","purple","black"),
       lty=c(2,1,2,1,1),cex=1.0)

# Add predicted to the graph
lines(cells ~ day,controlpredS,lwd=3,lty=2,col=camacols["S"])
lines(cells ~ day,ribopredS,lwd=3,col=camacols["S"])
lines(cells ~ day,treatpredS,lwd=3,lty=2,col="purple")
lines(cells ~ day,bothpredS,lwd=3,col="purple")
lines(cells ~ day,nullpredS,lwd=3,col="black")

# Add arrow
    lday <- 5
    fudge <- 0.02
    fp <- bothpredS$cells[lday]
    lp <- nullpredS$cells[lday]
    dy <- max(days)+dplus
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple4",code=1,lty=1,lwd=3,length=0.1)

# Do the stats
bothS$nullS <- nullpredS$cells
bothS$logratS <- with(bothS,log(cells/nullS))
summary(lm(logratS ~ day,bothS))
## (Intercept) -0.31465    0.05621    -5.6  8.7e-05 ***
## day          0.09491    0.00473    20.1  3.7e-11 ***

##################################################################
# Illustrate facilitation with other treatments (allfacil.R)
# Fulvestrant with lvl2=2: Point, no facilitation of S cells

# For some reason, ilvl1 and ilvl2 vanished.  Remake by hand
Fulvestrant3$ilvl1 <- 1
Fulvestrant3$ilvl1[Fulvestrant3$level1 > 0] <- 2
Fulvestrant3$ilvl2 <- 1
Fulvestrant3$ilvl2[Fulvestrant3$level2 > 0] <- 2
Fulvestrant3$ilvl2[Fulvestrant3$level2 > 1] <- 3
Fulvestrant3$ilvl2[Fulvestrant3$level2 > 3] <- 4

inm <- 2
  nm <- camafiles[inm]
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  lvl1 <- 1
  lvl2 <- lvl2s[inm]

# Pull out the S files
  controlS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                     ilvl1==1 & ilvl2==1)
  treatS <- subset(tfile,env=="S" & type=="S" & day %in% days &
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  competeS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                     ilvl1==1 & ilvl2==1)
  facilS <- subset(tfile,env=="M" & type=="S" & day %in% days &
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))

# Pull out the R files
  controlR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                     ilvl1==1 & ilvl2==1)
  treatR <- subset(tfile,env=="R" & type=="R" & day %in% days &
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  competeR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                     ilvl1==1 & ilvl2==1)
  facilR <- subset(tfile,env=="M" & type=="R" & day %in% days &
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))

# Pull out the fits
  controlS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                        ilvl1==1 & ilvl2==1)
  treatS.rK <- subset(tfile.rK,env=="S" & type=="S" & 
                      ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  competeS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                        ilvl1==1 & ilvl2==1)
  facilS.rK <- subset(tfile.rK,env=="M" & type=="S" & 
                      ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  
  controlR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                        ilvl1==1 & ilvl2==1)
  treatR.rK <- subset(tfile.rK,env=="R" & type=="R" & 
                      ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))
  competeR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                        ilvl1==1 & ilvl2==1)
  facilR.rK <- subset(tfile.rK,env=="M" & type=="R" & 
                    ilvl1==max(c(1,lvl1)) & ilvl2==max(c(1,lvl2)))

# Make the predictions
  controlpredS <- data.frame(day=days,
                  cells=lfun(mean(controlS$cells[controlS$day==4]),
                                  controlS.rK$K,controlS.rK$r,days))

  treatpredS <- data.frame(day=days,
                  cells=lfun(mean(treatS$cells[treatS$day==4]),
                                  treatS.rK$K,treatS.rK$r,days))

  controlpredR <- data.frame(day=days,
                  cells=lfun(mean(controlR$cells[controlR$day==4]),
                                  controlR.rK$K,controlR.rK$r,days))

  treatpredR <- data.frame(day=days,
                  cells=lfun(mean(treatR$cells[treatR$day==4]),
                                  treatR.rK$K,treatR.rK$r,days))

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
  nullparms["rS"] <- (treatS.rK$r)*(competeS.rK$r)/controlS.rK$r
  nullparms["KS"] <- (treatS.rK$K)*(competeS.rK$K)/controlS.rK$K
  nullparms["rR"] <- (treatR.rK$r)*(competeR.rK$r)/controlR.rK$r
  nullparms["KR"] <- (treatR.rK$K)*(competeR.rK$K)/controlR.rK$K
  LVout <- as.data.frame(ode(y=init, parms=nullparms,times=days, func=LVmodel))
  nullpredS <- LVout$S
  nullpredR <- LVout$R

# Plot up the data for S
  ymax <- max(c(controlS$cells,treatS$cells,competeS$cells,facilS$cells))
  plot(cells ~ day,controlS,pch="",col=camacols["S"],
       xlim=c(4,xmax),ylim=c(0,ymax),
       xlab="Day",ylab="Cell number",cex.lab=1.7,cex.axis=1.3)
  points(cells ~ day,controlS,pch=1,col=camacols["S"])
  points(cells ~ day,treatS,pch=19,col=camacols["S"])
  points(cells ~ day,competeS,pch=1,col="purple")
  points(cells ~ day,facilS,pch=19,col="purple")
  legend("topleft",c("Monoculture: Untreated",
                     "Monoculture: Fulvestrant",
                     "Coculture: Untreated",
                     "Coculture: Fulvestrant",
                     "Expected combined effect"),lwd=3,seg.len=4,
         col=c(camacols["S"],camacols["S"],"purple","purple","black"),
         lty=c(2,1,2,1,1),cex=1.0)

# Add predicted to the graph
  lines(cells ~ day,controlpredS,lwd=3,lty=2,col=camacols["S"])
  lines(cells ~ day,treatpredS,lwd=3,col=camacols["S"])
  lines(days,competepredS,lwd=3,lty=2,col="purple")
  lines(days,facilpredS,lwd=3,col="purple")
  lines(days,nullpredS,lwd=3,col="black")

# Add arrow
    lday <- 5
    fudge <- 0.02
    fp <- facilpredS[lday]
    lp <- nullpredS[lday]
    dy <- max(days)+dplus
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple4",code=1,lty=1,lwd=3,length=0.1)

# Do the stats
facilS$nullS <- nullpredS
facilS$logratS <- with(facilS,log(cells/nullS))
summary(lm(logratS ~ day,facilS))
## (Intercept)   0.0320     0.2179    0.15  0.88544    
## day           0.0849     0.0183    4.63  0.00047 ***

########################################################################
# Illustrate facilitation modification (allfacilmod.R)
# Estradiol1 with lvl2=2: Point, facilitation of S cells cancelled

inm <- 1
  nm <- camafiles[inm]
  tfile <- get(nm)
  tfile.rK <- get(paste0(nm,".rK"))
  lvl1 <- 2
  lvl2  <- lvl2s[inm]

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

# Plot up the data for S
  ymax <- max(c(controlS$cells,riboS$cells,competeS$cells,facilS$cells))
  plot(cells ~ day,controlS,pch="",col=camacols["S"],
       xlim=c(4,xmax),ylim=c(0,ymax),
       xlab="Day",ylab="Cell number",cex.lab=1.7,cex.axis=1.3)
  points(cells ~ day,controlS,pch=1,col=camacols["S"])
  points(cells ~ day,riboS,pch=19,col=camacols["S"])
  points(cells ~ day,competeS,pch=1,col="purple")
  points(cells ~ day,facilS,pch=19,col="purple")
  legend("topleft",c("Monoculture: Estradiol",
                     "Monoculture: Ribociclib and estradiol",
                     "Coculture: Estradiol",
                     "Coculture: Ribociclib and estradiol",
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
    dy <- max(days)+dplus
    arrows(dy,fp+fudge*(lp-fp),dy,fp+(1-fudge)*(lp-fp),
          col="purple4",code=1,lty=1,lwd=3,length=0.1)

# Do the stats
facilS$nullS <- nullpredS
facilS$logratS <- with(facilS,log(cells/nullS))
summary(lm(logratS ~ day,facilS))
## (Intercept) -0.06062    0.03812   -1.59     0.14    
## day          0.02780    0.00321    8.67  9.2e-07 ***
# Oops, but we really want to show a reduction in facilitation.
# This has to be compared with the facilitation in the absence of
# Estradiol.
estrafacilS <- facilS
# Just remake facilS without Estradiol
source("noestra.R")
noestrafacilS <- facilS
estrafacilS$flag <- "E"
noestrafacilS$flag <- "N"
tmp <- rbind(estrafacilS,noestrafacilS)
summary(lm(logratS ~ day*flag,tmp))
## (Intercept) -0.06062    0.10796   -0.56   0.5793    
## day          0.02780    0.00909    3.06   0.0051 ** 
## flagN       -0.34713    0.15268   -2.27   0.0315 *  
## day:flagN    0.12012    0.01285    9.35  8.5e-10 ***


# Make a boxplot for ribo and for the modifiers
# Set some sizes for the graphs
faxissize <- 1.2
flabsize <- 1.3
fnamesize <- 1.7
fmainsize <- 1.3
yfudge <- 0.3

Syrange <- c(-2,6)
Ryrange <- c(-2,6)
par(mar=c(9,5,4,2))

boxplot(log(cells/null) ~ doselabel,camafacilf,subset=type=="R",
        col=camacols["R"],
        ylim=Ryrange,ylab="Strength of facilitation",
        xaxt="n",xlab="",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
title(main="Facilitation of resistant cells under treatment",cex.main=fmainsize)
abline(h=0,lwd=2,lty=2,col="gray")
axis(1, labels = FALSE)
labels <- sort(unique(camafacilf$doselabel))
text(x =  seq_along(labels), y = par("usr")[3]-yfudge, srt = 45, adj = 1,
      labels = labels, xpd = TRUE,cex=flabsize)

boxplot(log(cells/null) ~ doselabel,camafacilmodf,subset=type=="R",
        col=camacols["R"],
        ylim=Ryrange,ylab="Strength of facilitation",
        xaxt="n",xlab="",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
title(main="Facilitation of resistant cells under \n combination treatment with ribociclib",
      cex.main=fmainsize)
abline(h=0,lwd=2,lty=2,col="gray")
axis(1, labels = FALSE)
tmplabels <- sort(unique(camafacilmodf$doselabel))
tmplabels <- as.character(tmplabels)
tmplabels[1] <- "No modifier"
labels <- paste0(tmplabels," + \n","Ribociclib 400nM")
# labels <- paste0(tmplabels," + ",\n,"Ribociclib 400nM")
text(x =  seq_along(labels), y = par("usr")[3]-yfudge, srt = 45, adj = 1,
      labels = labels, xpd = TRUE,cex=flabsize)

par(mar=c(5,5,4,2))
dev.off()

# Stats on facilitation: This compares everything with Ribociclib 200
# All but Fulvestrant show less facilitation of S cells
summary(lm(log(cells/null) ~ doselabel,camafacilf,subset=type=="S"))
## (Intercept)                  0.783      0.193    4.05  0.00014 ***
## doselabelRibociclib 400nM    1.953      0.306    6.38  1.9e-08 ***
## doselabelEstradiol 0.1nM    -0.883      0.387   -2.28  0.02559 *  
## doselabelFulvestrant 1nM     0.524      0.387    1.35  0.18044    
## doselabelFulvestrant 3nM     1.660      0.387    4.29  5.9e-05 ***

summary(lm(log(cells/null) ~ doselabel,camafacilf,subset=type=="R"))
## (Intercept)                 -0.214      0.076   -2.81   0.0065 ** 
## doselabelRibociclib 400nM   -0.147      0.120   -1.23   0.2247    
## doselabelEstradiol 0.1nM     0.823      0.152    5.42  8.9e-07 ***
## doselabelFulvestrant 1nM    -0.521      0.152   -3.43   0.0011 ** 
## doselabelFulvestrant 3nM    -0.991      0.152   -6.52  1.1e-08 ***

# Now compare with 0
facilsummf <- data.frame(doselabel=as.character(sort(unique(camafacilf$doselabel))),
                        estimateS=0,stderrS=0,pS=0,estimateR=0,stderrR=0,pR=0)
for (i in 1:nrow(facilsummf)) {
  nm <- facilsummf$doselabel[i]
  tmp.t <- with(subset(camafacilf,type=="S" & doselabel==nm),t.test(log(cells/null)))
  facilsummf$estimateS[i] <- tmp.t$estimate
  facilsummf$stderrS[i] <- tmp.t$stderr
  facilsummf$pS[i] <- tmp.t$p.value
  tmp.t <- with(subset(camafacilf,type=="R" & doselabel==nm),t.test(log(cells/null)))
  facilsummf$estimateR[i] <- tmp.t$estimate
  facilsummf$stderrR[i] <- tmp.t$stderr
  facilsummf$pR[i] <- tmp.t$p.value
}
facilsummf
##        doselabel estimateS  stderrS         pS estimateR  stderrR         pR
## Ribociclib 200nM  0.783154 0.064993 3.8231e-12  -0.21363 0.053845 0.00050888
## Ribociclib 400nM  2.736022 0.449230 1.2005e-05  -0.36090 0.080194 0.00031548
##  Estradiol 0.1nM -0.100246 0.034038 1.8563e-02   0.60960 0.164194 0.00593341
##  Fulvestrant 1nM  1.306806 0.069703 6.7659e-08  -0.73437 0.162385 0.00194375
##  Fulvestrant 3nM  2.442962 0.187844 1.1587e-06  -1.20434 0.190286 0.00022561
# See a nice dose response for Ribociclib and Fulvestrant in S cells

# Stats on facilitation modification: This compares everything with Ribociclib
summary(lm(log(cells/null) ~ doselabel,camafacilmodf,subset=type=="S"))
# Facilitation decreased by Estradiol
## (Intercept)                 2.736      0.321    8.52  1.3e-10 ***
## doselabelEstradiol 0.1nM   -2.478      0.556   -4.46  6.3e-05 ***
## doselabelFulvestrant 1nM   -0.459      0.556   -0.82   0.4144    
## doselabelFulvestrant 3nM    1.720      0.556    3.09   0.0036 ** 

summary(lm(log(cells/null) ~ doselabel,camafacilmodf,subset=type=="R"))
## (Intercept)               -0.3609     0.0566   -6.38  1.3e-07 ***
## doselabelEstradiol 0.1nM   0.2590     0.0980    2.64   0.0116 *  
## doselabelFulvestrant 1nM  -0.2874     0.0980   -2.93   0.0055 ** 
## doselabelFulvestrant 3nM  -0.0867     0.0980   -0.88   0.3814    
# Don't really care about the comparison with ribociclib.

# Compare with 0
facilmodsummf <- data.frame(doselabel=as.character(sort(unique(camafacilmodf$doselabel))),
                        estimateS=0,stderrS=0,pS=0,estimateR=0,stderrR=0,pR=0)
for (i in 1:nrow(facilmodsummf)) {
  nm <- facilmodsummf$doselabel[i]
  tmp.t <- with(subset(camafacilmodf,type=="S" & doselabel==nm),t.test(log(cells/null)))
  facilmodsummf$estimateS[i] <- tmp.t$estimate
  facilmodsummf$stderrS[i] <- tmp.t$stderr
  facilmodsummf$pS[i] <- tmp.t$p.value
  tmp.t <- with(subset(camafacilmodf,type=="R" & doselabel==nm),t.test(log(cells/null)))
  facilmodsummf$estimateR[i] <- tmp.t$estimate
  facilmodsummf$stderrR[i] <- tmp.t$stderr
  facilmodsummf$pR[i] <- tmp.t$p.value
}
facilmodsummf
##         doselabel estimateS  stderrS         pS estimateR  stderrR         pR
## Ribociclib 400nM   2.73602 0.449230 1.2005e-05  -0.36090 0.080194 3.1548e-04
##  Estradiol 0.1nM   0.25761 0.018308 6.3187e-07  -0.10187 0.033549 1.6146e-02
##  Fulvestrant 1nM   2.27742 0.200827 3.2948e-06  -0.64833 0.048452 9.3110e-07
##  Fulvestrant 3nM   4.45609 0.398099 3.6366e-06  -0.44764 0.045055 8.9073e-06
# Facilitation of S cells just about gone with Estradiol, 
# Stronger with Fulvestrant

# The message:
# S cells are hammered by Ribociclib and Fulvestrant.  They do better
# than expected with R cells around because R cells make enough
# Estradiol to overcome the effects of the drugs.

# R cells are harmed by Estradiol itself, and do a little better with S
# cells around.  Not sure why, and why there is no dose response.
# R cells do worse than expected in coculture when treated 
# with Ribociclib or Fulvestrant because of competition with the
# facilitated S cells.

