Eagg <- aggregate(Estradiol ~ level1+env,Egood,mean)
Eagg$Epred2 <- aggregate(Epred2 ~ level1+env,Egood,mean)$Epred2
Eagg$Efudge <- Eagg$Estradiol
Eagg$Efudge[3] <- Eagg$Efudge[3] + 2

pdf("Fig2c.pdf")
par(mar=c(5,5,4,2))
plot(Estradiol ~ Epred2,Egood,pch=1,col=as.factor(level1),cex=0.6,
     subset=env=="R",
     xlab="Predicted Estradiol",
     ylab="Observed Estradiol",
     cex.lab=1.4,cex.axis=1.4,
     xlim=range(c(Estradiol,Epred2)),ylim=range(c(Estradiol,Epred2)))
points(Epred2 ~ Estradiol,Egood,pch=8,col=as.factor(level1),cex=0.6,
       subset=env=="S")
points(Epred2 ~ Estradiol,Egood,pch=19,col=as.factor(level1),cex=0.6,
       subset=env=="M")
points(Epred2 ~ Efudge,Eagg,pch=1,col=as.factor(level1),cex=2.0,
     subset=env=="R")
points(Epred2 ~ Efudge,Eagg,pch=8,col=as.factor(level1),cex=2.0,
       subset=env=="S")
points(Epred2 ~ Efudge,Eagg,pch=19,col=as.factor(level1),cex=2.0,
       subset=env=="M")
abline(0,1,lwd=2,col="gray",lty=2)
legend("topleft",c("Untreated","Treated"),pch=19, col=1:2,cex=1.3)
legend(25,117.5,c("S cells only","R cells only","Mixed"),pch=c(8,1,19),col=1,cex=1.3)
dev.off()

print(summary(lm(Estradiol ~ Epred2,Egood)))
## (Intercept) 0.002579   3.748648   0.001    0.999    
## Epred2      0.999955   0.057300  17.451   <2e-16 ***
## Adjusted R-squared:  0.8537 

print(summary(lm(Estradiol ~ Epred2-1,Egood)))

