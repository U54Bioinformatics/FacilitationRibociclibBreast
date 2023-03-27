Mnms <- paste0("M",1:4)
################################################################
############## Get some files ready ############################
################################################################
rKnms <- c("r","K")
rKnmsS <- paste0(rKnms,"S") 
rKnmsR <- paste0(rKnms,"R") 

mcf7.rK <- newMCF7.rK
mcf7dat <- newMCF7

gday <- 0
gdays <- c(gday,days)
mcf7compdat <- subset(mcf7dat,substr(env,1,1) != "M" & day==gday)
mcf7compdat <- mcf7compdat[,c("file","level1","level2","type","rep","cells")]
names(mcf7compdat)[names(mcf7compdat)=="cells"] <- paste0("cells",gday)

for (iday in 2:length(gdays)) {
  gday <- gdays[iday]
  tmp <- subset(mcf7dat,substr(env,1,1) != "M" & day==gday)
  mcf7compdat[,paste0("cells",gday)] <- tmp$cells
}

pdf("MCF7figs.pdf")
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

# newMCF7  
tmpS <- subset(mcf7.rK,ilvl1==1 & ilvl2==1 & env == "S0")
tmpR <- subset(mcf7.rK,ilvl1==1 & ilvl2==1 & env == "R0")
tnms <- c("file",rKnms)
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[names(tmpS) %in% rKnms] <- rKnmsS
names(tmpR)[names(tmpR) %in% rKnms] <- rKnmsR
mcf7costs <- cbind(tmpS,tmpR[,rKnmsR])
# This has only one row, so can't really test anything

#########################################
########## Benefit of resistance ########
#########################################

# newMCF7
tmpS <- subset(mcf7.rK,ilvl1>1 & ilvl2==1 & env == "S0")
tmpR <- subset(mcf7.rK,ilvl1>1 & ilvl2==1 & env == "R0")
tnms <- c("file","level1",rKnms)
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[names(tmpS) %in% rKnms] <- rKnmsS
names(tmpR)[names(tmpR) %in% rKnms] <- rKnmsR
mcf7benefits <- cbind(tmpS,tmpR[,rKnmsR])

# Not sure of the right stats
summary(lm(I(cells14/cells4) ~ type*level1,mcf7compdat))
summary(lm(I(cells14/cells4) ~ type+level1,mcf7compdat))
summary(lm(I(cells18/cells4) ~ type*level1,mcf7compdat))
summary(lm(I(cells18/cells4) ~ type+level1,mcf7compdat))
# For newMCF7, plot the growth curves for Figure 1
fudge <- 0.6
cfun <- function(x) 1.0 + fudge*(x==24) + 2*fudge*(x==50)
ymax <- 120000
plot(cells ~ day,newMCF7,subset=day < 21,
     ylim=c(0,ymax),pch="",col="white",
     xlab="Day",ylab="Cell Number",cex.lab=flabsize,cex.axis=faxissize)
tmpfile <- subset(newMCF7,substr(env,1,1) != "M" & level1>0 & level2==0)
tmpfile.rK <- subset(mcf7.rK,substr(env,1,1) != "M" & level1>0 & level2==0)
for (irow in 1:nrow(tmpfile.rK)) {
  tmp0 <- subset(tmpfile,level1==tmpfile.rK$level1[irow] & 
          type==tmpfile.rK$type[irow] & day < 21)
  tmp <- subset(tmp0,day > 0)
  days <- unique(tmp$day)
  Cpred <- lfun(mean(tmp$cells[tmp$day==4]),tmpfile.rK$K[irow],
          tmpfile.rK$r[irow],days)
  points(cells ~ day,tmp0,pch=19,col=mcf7cols[tmp0$type],
         cex=cfun(tmpfile.rK$level1[irow]))
  lines(days,Cpred,col=mcf7cols[tmp0$type],lwd=cfun(tmpfile.rK$level1[irow]))
}

tmpfile <- subset(newMCF7,substr(env,1,1) != "M" & level1==0 & level2==0)
tmpfile.rK <- subset(mcf7.rK,substr(env,1,1) != "M" & level1==0 & level2==0)
for (irow in 1:nrow(tmpfile.rK)) {
  tmp0 <- subset(tmpfile,level1==tmpfile.rK$level1[irow] & type==tmpfile.rK$type[irow] & day < 21)
  tmp <- subset(tmp0,day > 0)
  days <- unique(tmp$day)
  Cpred <- lfun(mean(tmp$cells[tmp$day==4]),tmpfile.rK$K[irow],tmpfile.rK$r[irow],days)
  points(cells ~ day,tmp0,pch=19,col=mcf7cols[tmp0$type],cex=0.6)
  lines(days,Cpred,col=mcf7cols[tmp0$type],lwd=2,lty=2)
}
legend("topleft",c("Sensitive","Resistant"),lwd=4,bty="n",
       col=mcf7cols,cex=1.3)
legend(-0.2,0.9*ymax,c("Treated 5.0 nM","Treated 2.4 nM","Treated 1.5 nM","Untreated"),
       col=1,bty="n",
       lty=c(1,1,1,2),lwd=c(3,2,1.5,1.5),pch=19,pt.cex=c(2,1.5,1,0.8),cex=1.3)

# Do the stats.  Remove all mixed, day 0 and day 21
tmpfile <- subset(newMCF7,substr(env,1,1) != "M" & day > 0 & day < 21)
summary(lm(cells ~ day*type*as.factor(level1),tmpfile))
summary(lm(cells ~ type*(day+as.factor(level1)),tmpfile))

just4 <- subset(tmpfile,day==4)[,c("level1","type","rep","cells")]
names(just4)[4] <- "cells4"
tmp <- merge(tmpfile,just4,by=c("level1","type","rep"))
tmp <- within(tmp,lograt <- log(cells/cells4))
summary(lm(exp(lograt) ~ type*as.factor(level1),tmp,subset=day==18))
## (Intercept)                 0.7968     0.0862    9.25  8.1e-08 ***
## typeS                       1.4089     0.1219   11.56  3.5e-09 ***
## as.factor(level1)15         0.0829     0.1219    0.68    0.506    
## as.factor(level1)24        -0.0975     0.1219   -0.80    0.435    
## as.factor(level1)50        -0.3065     0.1219   -2.51    0.023 *  
## typeS:as.factor(level1)15  -0.4998     0.1724   -2.90    0.010 *  
## typeS:as.factor(level1)24  -0.4891     0.1724   -2.84    0.012 *  
## typeS:as.factor(level1)50  -0.3661     0.1724   -2.12    0.050 *  

#########################################
########## Asymmetry of competition #####
#########################################

tmpS <- subset(mcf7.rK,ilvl1==1 & ilvl2==1 & type == "S" & env %in% Mnms)
tmpR <- subset(mcf7.rK,ilvl1==1 & ilvl2==1 & type == "R" & env %in% Mnms)
tnms <- c("file","alpha")
tmpS <- tmpS[,tnms]
tmpR <- tmpR[,tnms]
names(tmpS)[2] <- "alphaSR"
names(tmpR)[2] <- "alphaRS"
mcf7alphas <- cbind(tmpS,alphaRS=tmpR[,2])
ks.test(mcf7alphas$alphaRS,mcf7alphas$alphaSR)   #p-value = 1
wilcox.test(mcf7alphas$alphaRS,mcf7alphas$alphaSR)   #p-value = 1
t.test(mcf7alphas$alphaRS,mcf7alphas$alphaSR)   #p-value = 0.57
t.test(mcf7alphas$alphaRS,mcf7alphas$alphaSR,paired=TRUE)   #p-value = 0.52

boxplot(mcf7alphas[,c("alphaRS","alphaSR")],
        ylim=c(0,3.3),xaxt="n",xlab="",
        col=mcf7cols,
        ylab="Competitive Effect",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
abline(h=1,lty=3,col="gray")
abline(h=1,lty=3,col="gray")
title(main="Competitive Effects in MCF7 Cells",cex.main=fmainsize)
labels=c("Effect of sensitive \n on resistant", "Effect of resistant \n on sensitive")
text(x = 0.0+seq_along(labels), y = -0.45, srt = 0, adj = 0.5,
   labels = labels, xpd = TRUE,cex=flabsize)

#########################################
############# Facilitation ##############
#########################################

## newMCF7
mcf7facilitation <- NULL
for (ilvl in 2:4) {
  controlS.rK <- subset(mcf7.rK,env=="S0" & type=="S" & level1==0 & level2==0)
  tmp <- subset(mcf7.rK,env=="S0" & type=="S" & ilvl1 %in% ilvl & level2==0)
  riboS.rK  <- tmp[1,] 
  riboS.rK[,c("r","K")] <- apply(tmp[,c("r","K")],2,mean)
  competeS.rK <- subset(mcf7.rK,env=="M1" & type=="S" & level1==0 & level2==0)
  tmp <- subset(mcf7.rK,env=="M1" & type=="S" & ilvl1 %in% ilvl & level2==0)
  facilS.rK  <- tmp[1,] 
  facilS.rK[,c("r","K")] <- apply(tmp[,c("r","K")],2,mean)
  facilS <- subset(mcf7dat,env=="M1" & type=="S" & ilvl1 %in% ilvl & level2==0)
  tmpS <- aggregate(cells ~ day+file,subset(facilS,day > 0 & day < 21),mean)
  names(tmpS)[names(tmpS)=="cells"] <- "S"

  controlR.rK <- subset(mcf7.rK,env=="R0" & type=="R" & level1==0 & level2==0)
  tmp <- subset(mcf7.rK,env=="R0" & type=="R" & ilvl1 %in% ilvl & level2==0)
  riboR.rK  <- tmp[1,] 
  riboR.rK[,c("r","K")] <- apply(tmp[,c("r","K")],2,mean)
  competeR.rK <- subset(mcf7.rK,env=="M1" & type=="R" & level1==0 & level2==0)
  tmp <- subset(mcf7.rK,env=="M1" & type=="R" & ilvl1 %in% ilvl & level2==0)
  facilR.rK  <- tmp[1,] 
  facilR.rK[,c("r","K")] <- apply(tmp[,c("r","K")],2,mean)
  facilR <- subset(mcf7dat,env=="M1" & type=="R" & ilvl1 %in% ilvl & level2==0)
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
  nm <- "newMCF7"
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
  mcf7facilitation <- rbind(mcf7facilitation,tmpm)
}

# These look really nice
plot(logSrat ~ day,mcf7facilitation,pch=19,cex=0.5*ilvl1,
     col=mcf7cols["S"],ylim=c(-0.8,0.8),
     xlab="Day",ylab="Strength of Facilitation",
     cex.lab=flabsize,cex.axis=faxissize)
# lines(aggregate(logSrat ~ day,mcf7facilitation,mean),col=mcf7cols["S"],lwd=3)
lines(logSrat ~ day,mcf7facilitation,col=mcf7cols["S"],subset=ilvl1==2,lwd=1,lty=1)
lines(logSrat ~ day,mcf7facilitation,col=mcf7cols["S"],subset=ilvl1==3,lwd=2,lty=1)
lines(logSrat ~ day,mcf7facilitation,col=mcf7cols["S"],subset=ilvl1==4,lwd=3,lty=1)
points(logRrat ~ day,mcf7facilitation,pch=19,cex=0.5*ilvl1,col=mcf7cols["R"])
# lines(aggregate(logRrat ~ day,mcf7facilitation,mean),col=mcf7cols["R"],lwd=3)
lines(logRrat ~ day,mcf7facilitation,col=mcf7cols["R"],subset=ilvl1==2,lwd=1,lty=1)
lines(logRrat ~ day,mcf7facilitation,col=mcf7cols["R"],subset=ilvl1==3,lwd=2,lty=1)
lines(logRrat ~ day,mcf7facilitation,col=mcf7cols["R"],subset=ilvl1==4,lwd=3,lty=1)
abline(h=0,lty=3,col="gray")
title(main="Facilitation in MCF7 Cells",cex.main=fmainsize)
legend("topleft",c("Sensitive","Resistant"),bty="n",lwd=3,seg.len=4,
       col=c(mcf7cols["S"],mcf7cols["R"]),lty=1,cex=1.3)
legend("bottomleft",c("Treated 5.0 nM","Treated 2.4 nM","Treated 1.5 nM"),
       col=1,bty="n",seg.len=4,
       lty=1,lwd=c(3,2,1),pch=19,pt.cex=c(2,1.5,1,0.8),cex=1.3)

# Analysis
Smcf7 <-    mcf7facilitation[mcf7facilitation$ilvl1>=2,
                             c("day","file","S","nullS","logSrat")]
Rmcf7 <-    mcf7facilitation[mcf7facilitation$ilvl1>=2,
                             c("day","file","R","nullR","logRrat")]
names(Smcf7) <- c("day","file","cells","null","lograt")
names(Rmcf7) <- c("day","file","cells","null","lograt")
Smcf7$type <- "S"
Rmcf7$type <- "R"
SRmcf7 <- rbind(Smcf7,Rmcf7)
summary(lm(lograt ~ day*type,SRmcf7))  #No interaction at this dose
summary(lm(lograt ~ day+type,SRmcf7))  #No type effect at this dose
summary(lm(lograt ~ day,SRmcf7))  #No type effect at this dose
## (Intercept) -0.00704    0.09600   -0.07    0.943  
## day          0.02425    0.00808    3.00    0.017 *
# Both types benefit equally

Smcf7 <-    mcf7facilitation[mcf7facilitation$ilvl1==3,
                             c("day","file","S","nullS","logSrat")]
Rmcf7 <-    mcf7facilitation[mcf7facilitation$ilvl1==3,
                             c("day","file","R","nullR","logRrat")]
names(Smcf7) <- c("day","file","cells","null","lograt")
names(Rmcf7) <- c("day","file","cells","null","lograt")
Smcf7$type <- "S"
Rmcf7$type <- "R"
SRmcf7 <- rbind(Smcf7,Rmcf7)
summary(lm(lograt ~ day*type,SRmcf7))  #No interaction at this dose
summary(lm(lograt ~ day+type,SRmcf7))  
## (Intercept)  0.06749    0.09257    0.73   0.4896   
## day          0.02726    0.00719    3.79   0.0068 **
## typeS       -0.17737    0.07126   -2.49   0.0417 * 
# This is spurious because of the curvature

Smcf7 <-    mcf7facilitation[mcf7facilitation$ilvl1==4,
                             c("day","file","S","nullS","logSrat")]
Rmcf7 <-    mcf7facilitation[mcf7facilitation$ilvl1==4,
                             c("day","file","R","nullR","logRrat")]
names(Smcf7) <- c("day","file","cells","null","lograt")
names(Rmcf7) <- c("day","file","cells","null","lograt")
Smcf7$type <- "S"
Rmcf7$type <- "R"
SRmcf7 <- rbind(Smcf7,Rmcf7)
summary(lm(lograt ~ day*type,SRmcf7))
## (Intercept)   0.3775     0.1706    2.21   0.0689 . 
## day          -0.0447     0.0144   -3.12   0.0207 * 
## typeS        -0.6317     0.2412   -2.62   0.0397 * 
## day:typeS     0.0971     0.0203    4.78   0.0031 **

# Combined analysis with quadratic terms
Smcf7 <-    mcf7facilitation[,c("day","file","S","ilvl1","nullS","logSrat")]
Rmcf7 <-    mcf7facilitation[,c("day","file","R","ilvl1","nullR","logRrat")]
names(Smcf7) <- c("day","file","cells","ilvl1","null","lograt")
names(Rmcf7) <- c("day","file","cells","ilvl1","null","lograt")
Smcf7$type <- "S"
Rmcf7$type <- "R"
SRmcf7 <- rbind(Smcf7,Rmcf7)
summary(lm(lograt ~ day+I(day^2)*type,SRmcf7))
## (Intercept)    -0.21306    0.21423   -0.99    0.329  
## day             0.09093    0.04249    2.14    0.042 *
## I(day^2)       -0.00430    0.00194   -2.22    0.036 *
## typeS          -0.18255    0.12719   -1.44    0.164  
## I(day^2):typeS  0.00197    0.00071    2.78    0.010 *
# No dose effect, but we do have the quadratic effect with R curving
# down, but starting out steeper.

dev.off()
