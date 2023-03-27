Mnms <- paste0("M",1:2)
################################################################
############## Get some files ready ############################
################################################################
rKnms <- c("r","K")
rKnmsS <- paste0(rKnms,"S") 
rKnmsR <- paste0(rKnms,"R") 

LY2dat <- LY2

gday <- 0
gdays <- c(gday,days)
LY2compdat <- subset(LY2dat,substr(env,1,1) != "M" & day==gday)
LY2compdat <- LY2compdat[,c("file","level1","level2","type","rep","cells")]
names(LY2compdat)[names(LY2compdat)=="cells"] <- paste0("cells",gday)

for (iday in 2:length(gdays)) {
  gday <- gdays[iday]
  tmp <- subset(LY2dat,substr(env,1,1) != "M" & day==gday)
  LY2compdat[,paste0("cells",gday)] <- tmp$cells
}

pdf("LY2figs.pdf")
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

# LY2  
tmpS <- subset(LY2.rK,ilvl1==1 & ilvl2==1 & env == "S0")
tmpR <- subset(LY2.rK,ilvl1==1 & ilvl2==1 & env == "R0")
tnms <- c("file",rKnms)
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[names(tmpS) %in% rKnms] <- rKnmsS
names(tmpR)[names(tmpR) %in% rKnms] <- rKnmsR
LY2costs <- cbind(tmpS,tmpR[,rKnmsR])
# This has only one row, so can't really test anything

#########################################
########## Benefit of resistance ########
#########################################

# LY2
tmpS <- subset(LY2.rK,ilvl1>1 & ilvl2==1 & env == "S0")
tmpR <- subset(LY2.rK,ilvl1>1 & ilvl2==1 & env == "R0")
tnms <- c("file","level1",rKnms)
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[names(tmpS) %in% rKnms] <- rKnmsS
names(tmpR)[names(tmpR) %in% rKnms] <- rKnmsR
LY2benefits <- cbind(tmpS,tmpR[,rKnmsR])

# Not sure of the right stats
summary(lm(I(cells14/cells4) ~ type*level1,LY2compdat))
summary(lm(I(cells14/cells4) ~ type+level1,LY2compdat))
summary(lm(I(cells18/cells4) ~ type*level1,LY2compdat))
summary(lm(I(cells18/cells4) ~ type+level1,LY2compdat))
# For LY2, plot the growth curves for Figure 1
fudge <- 0.6
cfun <- function(x) 1.0 + fudge*(x==3) + 2*fudge*(x==5)
ymax <- 200000
plot(cells ~ day,LY2,subset=day < 21,
     ylim=c(0,ymax),pch="",col="white",
     xlab="Day",ylab="Cell Number",cex.lab=flabsize,cex.axis=faxissize)
tmpfile <- subset(LY2,substr(env,1,1) != "M" & level1>0 & level2==0)
tmpfile.rK <- subset(LY2.rK,substr(env,1,1) != "M" & level1>0 & level2==0)
for (irow in 1:nrow(tmpfile.rK)) {
  tmp0 <- subset(tmpfile,level1==tmpfile.rK$level1[irow] & 
          type==tmpfile.rK$type[irow] & day < 21)
  tmp <- subset(tmp0,day > 0)
  days <- unique(tmp$day)
  Cpred <- lfun(mean(tmp$cells[tmp$day==4]),tmpfile.rK$K[irow],
          tmpfile.rK$r[irow],days)
  points(cells ~ day,tmp0,pch=19,col=LY2cols[tmp0$type],
         cex=cfun(tmpfile.rK$level1[irow]))
  lines(days,Cpred,col=LY2cols[tmp0$type],lwd=cfun(tmpfile.rK$level1[irow]))
}

tmpfile <- subset(LY2,substr(env,1,1) != "M" & level1==0 & level2==0)
tmpfile.rK <- subset(LY2.rK,substr(env,1,1) != "M" & level1==0 & level2==0)
for (irow in 1:nrow(tmpfile.rK)) {
  tmp0 <- subset(tmpfile,level1==tmpfile.rK$level1[irow] & type==tmpfile.rK$type[irow] & day < 21)
  tmp <- subset(tmp0,day > 0)
  days <- unique(tmp$day)
  Cpred <- lfun(mean(tmp$cells[tmp$day==4]),tmpfile.rK$K[irow],tmpfile.rK$r[irow],days)
  points(cells ~ day,tmp0,pch=19,col=LY2cols[tmp0$type],cex=0.6)
  lines(days,Cpred,col=LY2cols[tmp0$type],lwd=2,lty=2)
}
legend("topleft",c("Sensitive","Resistant"),lwd=4,bty="n",
       col=LY2cols,cex=1.3)
legend(-0.2,0.9*ymax,c("Treated 5.0 nM","Treated 3.0 nM","Untreated"),
       col=1,bty="n",
       lty=c(1,1,2),lwd=c(3,2,1.5),pch=19,pt.cex=c(2,1,0.8),cex=1.3)

# Do the stats.  Remove all mixed, day 0 and day 21
tmpfile <- subset(LY2,substr(env,1,1) != "M" & day > 0 & day < 21)
summary(lm(cells ~ type*(day+as.factor(level1)),tmpfile))
## (Intercept)                -32523       8251   -3.94  0.00017 ***
## typeS                       72245      11669    6.19  2.3e-08 ***
## day                         10384        598   17.36  < 2e-16 ***
## as.factor(level1)3         -14278       7260   -1.97  0.05263 .  
## as.factor(level1)5         -23032       7260   -3.17  0.00213 ** 
## typeS:day                   -6326        846   -7.48  7.4e-11 ***
## typeS:as.factor(level1)3   -49437      10268   -4.81  6.6e-06 ***
## typeS:as.factor(level1)5   -42184      10268   -4.11  9.4e-05 ***

#########################################
########## Asymmetry of competition #####
#########################################

tmpS <- subset(LY2.rK,ilvl1==1 & ilvl2==1 & type == "S" & env %in% Mnms)
tmpR <- subset(LY2.rK,ilvl1==1 & ilvl2==1 & type == "R" & env %in% Mnms)
tnms <- c("file","alpha")
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[2] <- "alphaSR"
names(tmpR)[2] <- "alphaRS"
LY2alphas <- cbind(tmpS,alphaRS=tmpR[,2])
t.test(LY2alphas$alphaRS,LY2alphas$alphaSR,paired=TRUE)   #p-value = 0.038
# There are just two values here...

boxplot(LY2alphas[,c("alphaRS","alphaSR")],
        ylim=c(0,2),xaxt="n",xlab="",
        col=LY2cols,
        ylab="Competitive Effect",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
abline(h=1,lty=3,col="gray")
abline(h=1,lty=3,col="gray")
title(main="Competitive Effects in LY2 Cells",cex.main=fmainsize)
labels=c("Effect of sensitive \n on resistant", "Effect of resistant \n on sensitive")
text(x = 0.0+seq_along(labels), y = -0.45, srt = 0, adj = 0.5,
   labels = labels, xpd = TRUE,cex=flabsize)

#########################################
############# Facilitation ##############
#########################################

## LY2
LY2facilitation <- NULL
for (ilvl in 2:3) {
  controlS.rK <- subset(LY2.rK,env=="S0" & type=="S" & level1==0 & level2==0)
  tmp <- subset(LY2.rK,env=="S0" & type=="S" & ilvl1 %in% ilvl & level2==0)
  riboS.rK  <- tmp[1,] 
  riboS.rK[,c("r","K")] <- apply(tmp[,c("r","K")],2,mean)
  competeS.rK <- subset(LY2.rK,env=="M1" & type=="S" & level1==0 & level2==0)
  tmp <- subset(LY2.rK,env=="M1" & type=="S" & ilvl1 %in% ilvl & level2==0)
  facilS.rK  <- tmp[1,] 
  facilS.rK[,c("r","K")] <- apply(tmp[,c("r","K")],2,mean)
  facilS <- subset(LY2dat,env=="M1" & type=="S" & ilvl1 %in% ilvl & level2==0)
  tmpS <- aggregate(cells ~ day+file,subset(facilS,day > 0 & day < 21),mean)
  names(tmpS)[names(tmpS)=="cells"] <- "S"

  controlR.rK <- subset(LY2.rK,env=="R0" & type=="R" & level1==0 & level2==0)
  tmp <- subset(LY2.rK,env=="R0" & type=="R" & ilvl1 %in% ilvl & level2==0)
  riboR.rK  <- tmp[1,] 
  riboR.rK[,c("r","K")] <- apply(tmp[,c("r","K")],2,mean)
  competeR.rK <- subset(LY2.rK,env=="M1" & type=="R" & level1==0 & level2==0)
  tmp <- subset(LY2.rK,env=="M1" & type=="R" & ilvl1 %in% ilvl & level2==0)
  facilR.rK  <- tmp[1,] 
  facilR.rK[,c("r","K")] <- apply(tmp[,c("r","K")],2,mean)
  facilR <- subset(LY2dat,env=="M1" & type=="R" & ilvl1 %in% ilvl & level2==0)
  tmpR <- aggregate(cells ~ day+file,subset(facilR,day > 0 & day < 21),mean)
  names(tmpR)[names(tmpR)=="cells"] <- "R"

  tmpm <- cbind(ilvl1=ilvl,tmpS,R=tmpR$R)

# Pull out the parameters 
  facilparms <- controlS.rK[,c("file","level1")]
  facilparms$rS <- facilS.rK$r
  facilparms$KS <- facilS.rK$K
  facilparms$rR <- facilR.rK$r
  facilparms$KR <- facilR.rK$K
  facilparms$alphaSR <- facilS.rK$alphaval
  facilparms$alphaRS <- facilR.rK$alphaval

  nullparms <- facilparms
  nullparms["rS"] <- (riboS.rK$r)*(competeS.rK$r)/controlS.rK$r
  nullparms["KS"] <- (riboS.rK$K)*(competeS.rK$K)/controlS.rK$K
  nullparms["rR"] <- (riboR.rK$r)*(competeR.rK$r)/controlR.rK$r
  nullparms["KR"] <- (riboR.rK$K)*(competeR.rK$K)/controlR.rK$K

# Make the null models to compare

  nullpredS <- NULL
  nullpredR <- NULL
  nm <- "LY2"
    init <- with(tmpm,c(S[file==nm & day==4],
                             R[file==nm & day==4]))
    names(init) <- c("S","R")
    pnms <- c("rS","KS","rR","KR","alphaSR","alphaRS")
    parms <- nullparms[nullparms$file==nm,pnms]
    days <- tmpm$day[tmpm$file==nm]
    LVout <- as.data.frame(ode(y=init, parms=parms,times=days, func=LVmodel))
    nullpredS <- c(nullpredS,LVout$S)
    nullpredR <- c(nullpredR,LVout$R)

  tmpm$nullS <- nullpredS
  tmpm$nullR <- nullpredR

# Find the log of the ratios
  tmpm <- within(tmpm,logSrat <- log(S/nullS))
  tmpm <- within(tmpm,logRrat <- log(R/nullR))
  LY2facilitation <- rbind(LY2facilitation,tmpm)
}

# These look really nice
plot(logSrat ~ day,LY2facilitation,pch=19,cex=0.5*ilvl1,
     col=LY2cols["S"],ylim=c(-1.2,2.0),
     xlab="Day",ylab="Strength of Facilitation",
     cex.lab=flabsize,cex.axis=faxissize)
# lines(aggregate(logSrat ~ day,LY2facilitation,mean),col=LY2cols["S"],lwd=3)
lines(logSrat ~ day,LY2facilitation,col=LY2cols["S"],subset=ilvl1==2,lwd=1,lty=1)
lines(logSrat ~ day,LY2facilitation,col=LY2cols["S"],subset=ilvl1==3,lwd=2,lty=1)
lines(logSrat ~ day,LY2facilitation,col=LY2cols["S"],subset=ilvl1==4,lwd=3,lty=1)
points(logRrat ~ day,LY2facilitation,pch=19,cex=0.5*ilvl1,col=LY2cols["R"])
# lines(aggregate(logRrat ~ day,LY2facilitation,mean),col=LY2cols["R"],lwd=3)
lines(logRrat ~ day,LY2facilitation,col=LY2cols["R"],subset=ilvl1==2,lwd=1,lty=1)
lines(logRrat ~ day,LY2facilitation,col=LY2cols["R"],subset=ilvl1==3,lwd=2,lty=1)
lines(logRrat ~ day,LY2facilitation,col=LY2cols["R"],subset=ilvl1==4,lwd=3,lty=1)
abline(h=0,lty=3,col="gray")
title(main="Facilitation in LY2 Cells",cex.main=fmainsize)
legend("topleft",c("Sensitive","Resistant"),bty="n",lwd=3,seg.len=4,
       col=c(LY2cols["S"],LY2cols["R"]),lty=1,cex=1.3)
legend(4,1.5,c("Treated 5.0 nM","Treated 3.0 nM"),
       col=1,bty="n",seg.len=4,
       lty=1,lwd=c(3,2),pch=19,pt.cex=c(2,1.5),cex=1.3)

# Analysis
SLY2 <-    LY2facilitation[LY2facilitation$ilvl1==2,
                             c("day","file","S","nullS","logSrat")]
RLY2 <-    LY2facilitation[LY2facilitation$ilvl1==2,
                             c("day","file","R","nullR","logRrat")]
names(SLY2) <- c("day","file","cells","null","lograt")
names(RLY2) <- c("day","file","cells","null","lograt")
SLY2$type <- "S"
RLY2$type <- "R"
SRLY2 <- rbind(SLY2,RLY2)
summary(lm(lograt ~ day*type,SRLY2))  
## (Intercept)  -0.1908     0.2433   -0.78    0.463  
## day          -0.0254     0.0205   -1.24    0.260  
## typeS        -0.0970     0.3441   -0.28    0.787  
## day:typeS     0.0837     0.0290    2.89    0.028 *
# Yes. Interaction.  Faciliatation
summary(lm(lograt ~ day,SRLY2,subset=type=="S"))  
## (Intercept) -0.28786    0.07387   -3.90   0.0300 * 
## day          0.05830    0.00622    9.38   0.0026 **
t.test(SRLY2$lograt[SRLY2$type=="R" & SRLY2$day>0])
# p = 0.027
t.test(SRLY2$lograt[SRLY2$type=="S" & SRLY2$day>0])
# p = 0.081

SLY2 <-    LY2facilitation[LY2facilitation$ilvl1==3,
                             c("day","file","S","nullS","logSrat")]
RLY2 <-    LY2facilitation[LY2facilitation$ilvl1==3,
                             c("day","file","R","nullR","logRrat")]
names(SLY2) <- c("day","file","cells","null","lograt")
names(RLY2) <- c("day","file","cells","null","lograt")
SLY2$type <- "S"
RLY2$type <- "R"
SRLY2 <- rbind(SLY2,RLY2)
summary(lm(lograt ~ day*type,SRLY2))  
## (Intercept)  -0.2565     0.4381   -0.59    0.580  
## day          -0.0334     0.0369   -0.91    0.400  
## typeS        -0.6654     0.6195   -1.07    0.324  
## day:typeS     0.1832     0.0521    3.51    0.013 *
summary(lm(lograt ~ day,SRLY2,subset=type=="S"))  
## (Intercept)  -0.9219     0.4009   -2.30    0.105  
## day           0.1498     0.0337    4.44    0.021 *
t.test(SRLY2$lograt[SRLY2$type=="R" & SRLY2$day>0])
# p = 0.031
t.test(SRLY2$lograt[SRLY2$type=="S" & SRLY2$day>0])
# p = 0.16

dev.off()
