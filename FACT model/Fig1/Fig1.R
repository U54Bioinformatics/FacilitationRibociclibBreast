pdf("Fig1.pdf")
par(mar=c(5,5,4,2),font=2, font.axis=2, font.lab=2,
#    font.main=2,family="Arial")
    font.main=2)
faxissize <- 1.2
flabsize <- 1.7
fnamesize <- 1.7
fmainsize <- 1.7
#########################################
########## Costs of resistance ##########
#########################################

# CAMA-1
tmpS <- subset(cama.rK,ilvl1==1 & ilvl2==1 & env == "S")
tmpR <- subset(cama.rK,ilvl1==1 & ilvl2==1 & env == "R")
tnms <- c("file",rKnms)
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[names(tmpS) %in% rKnms] <- rKnmsS
names(tmpR)[names(tmpR) %in% rKnms] <- rKnmsR
camacosts <- cbind(tmpS,tmpR[,rKnmsR])

# Test for differences
t.test(camacosts$rS,camacosts$rR,paired=TRUE) #p=0.421
t.test(camacosts$KS,camacosts$KR,paired=TRUE) #p=0.0048
# Resistant cells mainly have lower K

#########################################
########## Benefit of resistance ########
#########################################

tmpS <- subset(cama.rK,ilvl1==2 & ilvl2==1 & env == "S")
tmpR <- subset(cama.rK,ilvl1==2 & ilvl2==1 & env == "R")
tnms <- c("file","level1",rKnms)
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[names(tmpS) %in% rKnms] <- rKnmsS
names(tmpR)[names(tmpR) %in% rKnms] <- rKnmsR
camabenefits <- cbind(tmpS,tmpR[,rKnmsR])
# Stats are annoying because Fulvestrant3 cells just quite growing.
# We'd like to show that S cells take a larger reduction in growth than
# R cells, but might need to use growth until some target day, like day
# 11, to do this

# Regressions show that there is an interaction, visible by day 7 and
# strongest at day 14
summary(lm(I(cells4/cells0) ~ type*(level1 > 0),camacompdat,subset=level2==0))
## (Intercept)             1.429      0.107   13.31   <2e-16 ***
## typeS                   0.221      0.152    1.45     0.15    
## level1 > 0TRUE         -0.144      0.152   -0.95     0.35    
## typeS:level1 > 0TRUE   -0.354      0.215   -1.65     0.11    

summary(lm(I(cells7/cells0) ~ type*(level1 > 0),camacompdat,subset=level2==0))
## (Intercept)             4.111      0.515    7.98  4.3e-10 ***
## typeS                   0.557      0.729    0.76    0.448    
## level1 > 0TRUE         -0.727      0.729   -1.00    0.324    
## typeS:level1 > 0TRUE   -2.237      1.031   -2.17    0.035 *  

summary(lm(I(cells11/cells0) ~ type*(level1 > 0),camacompdat,subset=level2==0))
## (Intercept)             11.96       1.19   10.01  6.5e-13 ***
## typeS                    5.82       1.69    3.45   0.0013 ** 
## level1 > 0TRUE          -3.19       1.69   -1.89   0.0656 .  
## typeS:level1 > 0TRUE   -12.23       2.39   -5.12  6.6e-06 ***

summary(lm(I(cells14/cells0) ~ type*(level1 > 0),camacompdat,subset=level2==0))
## (Intercept)             28.77       2.29   12.55  3.9e-16 ***
## typeS                   12.47       3.24    3.85  0.00038 ***
## level1 > 0TRUE         -10.04       3.24   -3.10  0.00339 ** 
## typeS:level1 > 0TRUE   -27.22       4.58   -5.94  4.1e-07 ***

summary(lm(I(cells18/cells0) ~ type*(level1 > 0),camacompdat,subset=level2==0))
## (Intercept)             50.02       4.04   12.37  6.4e-16 ***
## typeS                   20.34       5.72    3.56  0.00092 ***
## level1 > 0TRUE         -17.29       5.72   -3.02  0.00417 ** 
## typeS:level1 > 0TRUE   -43.72       8.09   -5.41  2.5e-06 ***

# This looks great -- at all times we see a significant interaction with
# S cells growing more slowly with treatment and faster without

# Check the effects of the experiment, keeping level1 as a number
summary(lm(I(cells14/cells0) ~ type*level1+as.factor(file),camacompdat,subset=level2==0))
## (Intercept)                  27.9459     3.2819    8.52  1.3e-10 ***
## typeS                         9.2110     3.2496    2.83   0.0071 ** 
## level1                       -0.0369     0.0105   -3.51   0.0011 ** 
## as.factor(file)FOHTamoxifen   5.8314     3.4923    1.67   0.1026    
## as.factor(file)Fulvestrant3   2.4340     3.4082    0.71   0.4792    
## as.factor(file)Raloxifene3   -2.8930     3.4923   -0.83   0.4122    
## typeS:level1                 -0.0690     0.0145   -4.75  2.5e-05 ***
# No significant effect of file on the ratio

# Ill-advised check of three-way interaction
summary(lm(I(cells14/cells0) ~ type*level1*as.factor(file),camacompdat,subset=level2==0))
## (Intercept)                               30.39126    1.86591   16.29  < 2e-16 ***
## typeS                                     -2.12058    2.63879   -0.80  0.42755    
## level1                                    -0.04337    0.00660   -6.57  2.1e-07 ***
## as.factor(file)FOHTamoxifen               12.53026    2.63879    4.75  4.1e-05 ***
## as.factor(file)Fulvestrant3               -6.42288    2.63879   -2.43  0.02069 *  
## as.factor(file)Raloxifene3               -12.60493    2.63879   -4.78  3.8e-05 ***
## typeS:level1                              -0.02399    0.00933   -2.57  0.01500 *  
## typeS:as.factor(file)FOHTamoxifen          7.07683    3.73181    1.90  0.06697 .  
## typeS:as.factor(file)Fulvestrant3         21.75138    3.73181    5.83  1.8e-06 ***
## typeS:as.factor(file)Raloxifene3          29.51902    3.73181    7.91  5.0e-09 ***
## level1:as.factor(file)FOHTamoxifen        -0.01966    0.01475   -1.33  0.19200    
## level1:as.factor(file)Fulvestrant3         0.01754    0.00933    1.88  0.06917 .  
## level1:as.factor(file)Raloxifene3          0.04401    0.01475    2.98  0.00542 ** 
## typeS:level1:as.factor(file)FOHTamoxifen  -0.13322    0.02086   -6.39  3.6e-07 ***
## typeS:level1:as.factor(file)Fulvestrant3  -0.05528    0.01319   -4.19  0.00021 ***
## typeS:level1:as.factor(file)Raloxifene3   -0.15676    0.02086   -7.51  1.5e-08 ***
# That was ill-advised.

# Graph growth curves with Estradiol
nm <- "Estradiol1"
tfile <- get(nm)
tfile <- subset(tfile,day < 21)
tfile.rK <- get(paste0(nm,".rK"))
plot(cells ~ day,tfile,col="white",
     xlab="Day",ylab="Cell Number",cex.lab=flabsize,cex.axis=faxissize)
tmpfile <- subset(tfile,env != "M" & ilvl1==2 & level2==0)
tmpfile.rK <- subset(tfile.rK,env != "M" & ilvl1==2 & level2==0)
for (irow in 1:nrow(tmpfile.rK)) {
  tmp0 <- subset(tmpfile,type ==tmpfile.rK$type[irow])
  tmp <- subset(tmp0,day > 0)
  days <- unique(tmp$day)
  Cpred <- lfun(mean(tmp$cells[tmp$day==4]),tmpfile.rK$K[irow],tmpfile.rK$r[irow],days)
  points(cells ~ day,tmp0,pch=19,col=camacols[tmp0$type],cex=1.3)
  lines(days,Cpred,lwd=3,col=camacols[tmp0$type],lty=1)
}

tmpfile <- subset(tfile,env != "M" & level1==0 & level2==0)
tmpfile.rK <- subset(tfile.rK,env != "M" & level1==0 & level2==0)
for (irow in 1:nrow(tmpfile.rK)) {
  tmp0 <- subset(tmpfile,type ==tmpfile.rK$type[irow])
  tmp <- subset(tmp0,day > 0)
  days <- unique(tmp$day)
  Cpred <- lfun(mean(tmp$cells[tmp$day==4]),tmpfile.rK$K[irow],tmpfile.rK$r[irow],days)
  points(cells ~ day,tmp0,pch=19,col=camacols[tmp0$type],cex=0.6)
  lines(days,Cpred,lwd=2,lty=2,col=camacols[tmp0$type])
}
legend(-1.0,210000,c("Sensitive","Resistant"),lwd=4,bty="n",
     col=camacols,cex=1.3,seg.len=4)
legend(-0.5,175000,c("Treated","Untreated"),bty="n",
     col=1,lty=c(1,2),lwd=c(3,2),pch=19,pt.cex=c(1.5,1.0),cex=1.3,seg.len=4)

# Make a graphable version
camacomprat <- subset(camacompdat,level2==0)
camacomprat$type[camacomprat$type=="S"] <- "B"  #Stupid way to get S before R
par(mar=c(9,5,4,2))
boxplot(log(cells14/cells0) ~ type+(level1>0),camacomprat,
        col=c(camacols,camacols),
        ylab="Growth by Day 14",xaxt="n",xlab="",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
axis(1, labels = FALSE)
# labels <- c("S untreated","R untreated","S treated","R treated")
# text(x =  seq_along(labels), y = par("usr")[3]-0.2, srt = 45, adj = 1,
#      labels = labels, xpd = TRUE,cex=flabsize)
labels <- c("Untreated","Treated")
text(x =  2*seq_along(labels)-0.5, y = par("usr")[3]-0.3, srt = 0, adj = 0.5,
      labels = labels, xpd = TRUE,cex=flabsize+0.2)
legend("left",c("Sensitive","Resistant"),fill=camacols,bty="n",cex=1.3)
par(mar=c(5,5,4,2))

summary(lm(I(cells14/cells0) ~ type*(level1 > 0),camacompdat,subset=level2==0))
## (Intercept)             28.77       2.29   12.55  3.9e-16 ***
## typeS                   12.47       3.24    3.85  0.00038 ***
## level1 > 0TRUE         -10.04       3.24   -3.10  0.00339 ** 
## typeS:level1 > 0TRUE   -27.22       4.58   -5.94  4.1e-07 ***

####################################################
##### Figure 1e, the pedagogical figure ############
####################################################
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
           pch=c(1,19,1,19,NA),lty=c(2,1,2,1,1),bty="n",cex=1.0)

    legend("topright",inset=c(0,-0.3),
                      c("a. Treatment cost",
                       "b. Competition cost",
                       "c. Observed combined cost",
                       "d. Expected combined cost",
                       "e. Facilitation"),
                        lwd=3,pch=rep(NA,5),
           col=c(camacols["S"],"purple","blue","black","purple4"),
           lty=rep(NA,5),bty="n",seg.len=2)
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

#########################################
########## Asymmetry of competition #####
#########################################

tmpS <- subset(cama.rK,ilvl1==1 & ilvl2==1 & type == "S" & env=="M")
tmpR <- subset(cama.rK,ilvl1==1 & ilvl2==1 & type == "R" & env=="M")
tnms <- c("file","alpha")
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[2] <- "alphaSR"
names(tmpR)[2] <- "alphaRS"
camaalphas <- cbind(tmpS,alphaRS=tmpR[,2])
ks.test(camaalphas$alphaRS,camaalphas$alphaSR)   #p-value = 0.029
wilcox.test(camaalphas$alphaRS,camaalphas$alphaSR)   #p-value = 0.029
t.test(camaalphas$alphaRS,camaalphas$alphaSR)   #p-value = 0.0023
t.test(camaalphas$alphaRS,camaalphas$alphaSR,paired=TRUE)   #p-value = 0.013

boxplot(camaalphas[,c("alphaRS","alphaSR")],
        ylim=c(0,2.1),xaxt="n",xlab="",
#        names=c(expression(alpha[RS],alpha[SR])),
        col=camacols,
        ylab="Competitive Effect",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
abline(h=1,lty=3,col="gray")
title(main="Competitive Effects in CAMA-1 Cells",cex.main=fmainsize)
labels=c("Effect of sensitive \n on resistant", "Effect of resistant \n on sensitive")
text(x = 0.0+seq_along(labels), y = -0.25, srt = 0, adj = 0.5,
   labels = labels, xpd = TRUE,cex=flabsize)
t.test(camaalphas[,2],camaalphas[,3])
# t = -6.53, df = 4.25, p-value = 0.0023

#########################################
############# Facilitation ##############
#########################################

controlS.rK <- subset(cama.rK,env=="S" & type=="S" & level1==0 & level2==0)
riboS.rK <- subset(cama.rK,env=="S" & type=="S" & level1>0 & level2==0)
competeS.rK <- subset(cama.rK,env=="M" & type=="S" & level1==0 & level2==0)
facilS.rK <- subset(cama.rK,env=="M" & type=="S" & level1>0 & level2==0)
facilS <- subset(camadat,env=="M" & type=="S" & level1>0 & level2==0)
tmpS <- aggregate(cells ~ day+file,subset(facilS,day > 0 & day < 21),mean)
names(tmpS)[names(tmpS)=="cells"] <- "S"

controlR.rK <- subset(cama.rK,env=="R" & type=="R" & level1==0 & level2==0)
riboR.rK <- subset(cama.rK,env=="R" & type=="R" & level1>0 & level2==0)
competeR.rK <- subset(cama.rK,env=="M" & type=="R" & level1==0 & level2==0)
facilR.rK <- subset(cama.rK,env=="M" & type=="R" & level1>0 & level2==0)
facilR <- subset(camadat,env=="M" & type=="R" & level1>0 & level2==0)
tmpR <- aggregate(cells ~ day+file,subset(facilR,day > 0 & day < 21),mean)
names(tmpR)[names(tmpR)=="cells"] <- "R"
camafacilitation <- cbind(tmpS,R=tmpR$R)

# Pull out the parameters 
facilparms <- controlS.rK[,c("file","level1")]
facilparms$rS <- facilS.rK$r
facilparms$KS <- facilS.rK$K
facilparms$rR <- facilR.rK$r
facilparms$KR <- facilR.rK$K
facilparms$alphaSR <- facilS.rK$alpha
facilparms$alphaRS <- facilR.rK$alpha

nullparms <- facilparms
nullparms["rS"] <- (riboS.rK$r)*(competeS.rK$r)/controlS.rK$r
nullparms["KS"] <- (riboS.rK$K)*(competeS.rK$K)/controlS.rK$K
nullparms["rR"] <- (riboR.rK$r)*(competeR.rK$r)/controlR.rK$r
nullparms["KR"] <- (riboR.rK$K)*(competeR.rK$K)/controlR.rK$K

# Make the null models to compare

nullpredS <- NULL
nullpredR <- NULL
for (nm in unique(camafacilitation$file)) {
  init <- with(camafacilitation,c(S[file==nm & day==4],
                           R[file==nm & day==4]))
  names(init) <- c("S","R")
  pnms <- c("rS","KS","rR","KR","alphaSR","alphaRS")
  parms <- nullparms[nullparms$file==nm,pnms]
  days <- camafacilitation$day[camafacilitation$file==nm]
  LVout <- as.data.frame(ode(y=init, parms=parms,times=days, func=LVmodel))
  nullpredS <- c(nullpredS,LVout$S)
  nullpredR <- c(nullpredR,LVout$R)
}
camafacilitation$nullS <- nullpredS
camafacilitation$nullR <- nullpredR

# Find the log of the ratios
camafacilitation <- within(camafacilitation,logSrat <- log(S/nullS))
camafacilitation <- within(camafacilitation,logRrat <- log(R/nullR))
labs <- substr(camafacilitation$file,1,1)
labs[camafacilitation$file=="FOHTamoxifen"] <- "T"
# These look really nice
plot(logSrat ~ day,camafacilitation,pch=19,cex=1,col=camacols["S"],ylim=c(-1,7),
     xlab="Day",ylab="Strength of Facilitation",
     cex.lab=flabsize,cex.axis=faxissize)
# with(camafacilitation,text(day,logSrat,labs,col=camacols["S"]))
lines(aggregate(logSrat ~ day,camafacilitation,mean),col=camacols["S"],lwd=3)
points(logRrat ~ day,camafacilitation,pch=19,cex=1,col=camacols["R"])
# with(camafacilitation,text(day,logRrat,labs,col=camacols["R"]))
lines(aggregate(logRrat ~ day,camafacilitation,mean),col=camacols["R"],lwd=3)
abline(h=0,lty=3,col="gray")
title(main="Facilitation in CAMA-1 Cells",cex.main=fmainsize)
legend("topleft",c("Sensitive","Resistant"),bty="n",lwd=3,
       col=c(camacols["S"],camacols["R"]),lty=1,cex=1.3)
Scama <-    camafacilitation[,c("day","file","S","nullS","logSrat")]
Rcama <-    camafacilitation[,c("day","file","R","nullR","logRrat")]
names(Scama) <- c("day","file","cells","null","lograt")
names(Rcama) <- c("day","file","cells","null","lograt")
Scama$type <- "S"
Rcama$type <- "R"
SRcama <- rbind(Scama,Rcama)
summary(lm(lograt ~ day*type,SRcama))
## (Intercept)   0.2032     0.4753    0.43  0.67153    
## day          -0.0309     0.0400   -0.77  0.44505    
## typeS        -0.9784     0.6722   -1.46  0.15417    
## day:typeS     0.2125     0.0566    3.76  0.00061 ***

dev.off()
