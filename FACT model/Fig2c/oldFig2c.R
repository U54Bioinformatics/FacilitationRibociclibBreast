Eagg <- aggregate(Estradiol ~ level1+env,Egood,mean)
Eagg$Epred2 <- aggregate(Epred2 ~ level1+env,Egood,mean)$Epred2
Eagg$Efudge <- Eagg$Estradiol
Eagg$Efudge[3] <- Eagg$Efudge[3] + 2

pdf("Fig2c.pdf")
par(mar=c(5,5,4,2))
plot(Epred2 ~ Estradiol,Egood,pch=1,col=as.factor(level1),cex=0.6,
     subset=env=="R",
     xlab="Observed Estradiol",
     ylab="Predicted Estradiol",
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

# The original version of Figure 2c was incomprehensible

# pdf("oldFig2c.pdf")
# Rmax <- max(Egood[,c("S","R")])
# par(mar=c(5,5,4,2))
# plot(Estradiol ~ R,Egood,pch=1,col=as.factor(level1),
#      xlab="Number of Resistant Cells",
#      cex.lab=1.4,cex.axis=1.4,
#      subset= R > 0 & S==0,
#      xlim=c(0,Rmax),ylim=range(c(Estradiol,Epred2)))
# points(Estradiol ~ S,Egood,pch=8,col=as.factor(level1),
#      subset= S > 0 & R==0)
# points(Estradiol ~ R,Egood,pch=19,col=as.factor(level1),
#      subset= R > 0 & S>0)
# Rfun <- function(R) Efun2(0,R)
# Sfun <- function(S) Efun2(S,0)
# curve(Rfun,add=TRUE)
# curve(Sfun,add=TRUE)
# points(Epred2 ~ R,Egood,pch=19,col=3)
# legend("topleft",c("Untreated","Treated","Model"),pch=19, col=1:3)
# legend(-10000,117.1,c("S cells only","R cells only","Mixed"),pch=c(8,1,19),col=1)
# legend("topright",c("R cells only","Mixed"),pch=c(1,19),col=1)
# dev.off()

