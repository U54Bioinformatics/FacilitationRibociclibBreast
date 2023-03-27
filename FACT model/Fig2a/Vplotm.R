# Plot up the results

cols <- c("chartreuse2","black","cyan3")
tfun <- function(x) I(x)
pdf("Vplot.pdf")
par(mfrow=c(1,2),pty="s",mar=c(5,5,4,2),font=2, font.axis=2, 
    font.lab=2,font.main=2)
faxissize <- 1.2
flabsize <- 1.7
fnamesize <- 1.7
fmainsize <- 1.7
tmp <- subset(Vdat,treat=="U")
plot(tfun(cells) ~ day,tmp,pch="",col="white",
     xlab="Day",ylab="Cell Number",cex.lab=flabsize,cex.axis=faxissize)
points(tfun(cells) ~ day,tmp,subset=env=="S",pch=19,col=cols[1])
points(tfun(cells) ~ day,tmp,subset=env=="M",pch=19,col=cols[2])
points(tfun(cells) ~ day,tmp,subset=env=="R",pch=19,col=cols[3])
lines(tfun(SUmean) ~ day,VU,lwd=2,col=cols[1])
lines(tfun(MUmean) ~ day,VU,lwd=2,col=cols[2])
lines(tfun(RUmean) ~ day,VU,lwd=2,col=cols[3])
legend("topleft",c("Sensitive","Mixed","Resistant"),lty=1,bty="n",
       lwd=2,col=cols[1:3])
title(main="Untreated \n conditioned media",cex.main=fmainsize)

tmp <- subset(Vdat,treat=="T")
plot(tfun(cells) ~ day,tmp,pch="",col="white",
     xlab="Day",ylab="",cex.lab=flabsize,cex.axis=faxissize)
points(tfun(cells) ~ day,tmp,subset=env=="S",pch=19,col=cols[1])
points(tfun(cells) ~ day,tmp,subset=env=="M",pch=19,col=cols[2])
points(tfun(cells) ~ day,tmp,subset=env=="R",pch=19,col=cols[3])
lines(tfun(STmean) ~ day,VT,lwd=2,col=cols[1])
lines(tfun(MTmean) ~ day,VT,lwd=2,col=cols[2])
lines(tfun(RTmean) ~ day,VT,lwd=2,col=cols[3])
title(main="Treated \n conditioned media",cex.main=fmainsize)
par(mfrow=c(1,1),pty="m")
dev.off()

# Ok, now can do some analysis

tfun <- function(x) log(x)
summary(lm(tfun(cells) ~ day*env,Vdat,subset=treat=="U"))
## (Intercept)  7.65056    0.07938   96.38   <2e-16 ***
## day          0.15635    0.00732   21.37   <2e-16 ***
## envR         0.13323    0.11226    1.19    0.241    
## envS         0.03510    0.11226    0.31    0.756    
## day:envR    -0.00756    0.01035   -0.73    0.469    
## day:envS    -0.02070    0.01035   -2.00    0.051 .  
# Feeble, cells from S do almost worse than those from R and M 

summary(lm(tfun(cells) ~ day*env,Vdat,subset=treat=="T"))
## (Intercept)  7.82012    0.06827  114.55  < 2e-16 ***
## day          0.09462    0.00629   15.03  < 2e-16 ***
## envR        -0.18451    0.09655   -1.91  0.06198 .  
## envS         0.21691    0.09655    2.25  0.02930 *  
## day:envR     0.02738    0.00890    3.08  0.00345 ** 
## day:envS    -0.03256    0.00890   -3.66  0.00063 ***
# Yay, cells from the S env do worse than M, while R do better

