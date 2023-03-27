# The growth curves for MCF7 cells are tricky

mcf7cols <- c("chartreuse2","red")
names(mcf7cols) <- c("S","R")
gmn <- read.csv("../DATA/growth_curves_feb25_2021.csv",header=TRUE)

tmp <- t(gmn[,3:ncol(gmn)])
tmp <- as.data.frame(tmp)
tmp <- cbind(Cells=as.numeric(substr(rownames(tmp),2,nchar(rownames(tmp)))),tmp)
rownames(tmp) <- NULL
pcts <- c(100,90,80,70,60,0)
colnms <- paste0("S",rep(pcts,each=3),"_",rep(1:3,length(pcts)))
names(tmp) <- c("Cells",colnms)
gmn <- tmp

# The sprawling R cells are glued back together by the S cells, so we
# can use an overall correction based solely on the S cells.

mean2 <- function(x) mean(na.omit(x))
S100nms <- c("S100_1","S100_2","S100_3")
S90nms <- c("S90_1","S90_2","S90_3")
S80nms <- c("S80_1","S80_2","S80_3")
S70nms <- c("S70_1","S70_2","S70_3")
S60nms <- c("S60_1","S60_2","S60_3")
S0nms <- c("S0_1","S0_2","S0_3")

gmn$S100ave <- apply(gmn[,S100nms],1,mean2)
gmn$S90ave <- apply(gmn[,S90nms],1,mean2)
gmn$S80ave <- apply(gmn[,S80nms],1,mean2)
gmn$S70ave <- apply(gmn[,S70nms],1,mean2)
gmn$S60ave <- apply(gmn[,S60nms],1,mean2)
gmn$S0ave <- apply(gmn[,S0nms],1,mean2)
gmn$Save <- apply(gmn[,c(S100nms,S90nms,S80nms,S70nms,S60nms)],1,mean2)

newgm <- gmn[,c("Cells","S100ave","S90ave","S80ave","S70ave","S60ave","S0ave","Save")]
  
# Fit each of the sets with a Michaelis-Menten curve
gmells <- gmn$Cells
gfit <- function(x) {
 K <- x[1]
 A <- x[2]
 sum((A*tmp$Cells/(K+tmp$Cells)-tmp$Area)^2)
}

fitdat <- data.frame(pct=pcts,K=0,A=0)
for (ipct in 1:length(pcts)) {
  pct <- pcts[ipct]
  gtnms <- paste0("S",pct,"_",1:3)
  tmp <- matrix(unlist(gmn[,gtnms]),3*length(gmells),1)
  tmp <- data.frame(Cells=rep(gmells,3),Area=tmp)
  tmp <- subset(tmp,is.na(tmp$Area)==FALSE)
  Kstart <- mean(gmells)
  Astart <- max(tmp$Area)
  tmp.fit <- optim(c(Kstart,Astart),gfit)
  fitdat[ipct,c("K","A")] <-tmp.fit$par
}

# Fit to all the points with S cells
gtnms <- paste0("S",rep(pcts,each=3),"_",1:3)
tmp <- matrix(unlist(gmn[,gtnms]),3*length(pcts)*length(gmells),1)
tmp <- data.frame(Cells=rep(gmells,3*length(pcts)),Area=tmp)
tmp <- subset(tmp,is.na(tmp$Area)==FALSE)
Kstart <- mean(gmells)
Astart <- max(tmp$Area)
all.fit <- optim(c(Kstart,Astart),gfit)$par

cols <- c(mcf7cols["S"],"gray30","gray40","gray50","gray60",mcf7cols["R"])

# Superimpose the MM fits
mmfun <- function(x) A*x/(K+x)
for (ipct in 1:length(pcts)) {
  K <- fitdat$K[ipct]
  A <- fitdat$A[ipct]
}
K <- all.fit[1]
A <- all.fit[2]

# Now define it as a function
AS <- function(C) all.fit[2]*C/(all.fit[1]+C)
AR <- function(C) fitdat$A[fitdat$pct==0]*C/(fitdat$K[fitdat$pct==0]+C)

# Might as well save the inverse
ASinv <- function(A) all.fit[1]*A/(all.fit[2]-A)
ARinv <- function(A) (fitdat$K[fitdat$pct==0]*A)/(fitdat$A[fitdat$pct==0]-A)

# Huzzah, now we can use these in FACT!


############################
########  Areas ############
############################

# This has a different structure, and we'll just make the ones with
# different proportions into different replicates

days <- c(0, 4, 7, 11, 14, 18)
Areasm <- read.csv("../DATA/AreanewMCF7.csv",header=FALSE)
Fractionsm <- read.csv("../DATA/FractionnewMCF7.csv",header=FALSE)
Sfluorsm <- read.csv("../DATA/MCF7_Sfluor.csv",header=FALSE)
Rfluorsm <- read.csv("../DATA/MCF7_Rfluor.csv",header=FALSE)
Areasm[,1] <- days
Fractionsm[,1] <- days
Sfluorsm[,1] <- days
Rfluorsm[,1] <- days

# Fix the names: This has 4 ribo doses, and 7 initial conditions
tmp <- c("day",rep(c(rep("R000",3),
               rep("R015",3),
               rep("R024",3),
               rep("R050",3)),7))
tmp2 <- c("",c(rep(".S0",12),rep(".R0",12),rep(paste0(".M",1:5),each=12)))
tmp3 <- c("",rep(1:3,28))
colnms <- paste0(tmp,tmp2,tmp3)
names(Areasm) <- colnms
names(Fractionsm) <- colnms
names(Sfluorsm) <- colnms
names(Rfluorsm) <- colnms

Scolnms <- colnms[substr(colnms,6,6)=="S"]
Rcolnms <- colnms[substr(colnms,6,6)=="R"]
Mcolnms <- colnms[substr(colnms,6,6)=="M"]

# Fluorescence analysis from the LY2 file
# Replace the NA's with 0's in the fluorescence data
Sfluorsm[is.na(Sfluorsm)] <- 0
Rfluorsm[is.na(Rfluorsm)] <- 0

# Make fractions and see what happens
rawFractionsm <- Areasm
rawFractionsm <- Sfluorsm/(Sfluorsm+Rfluorsm)
rawFractionsm$day <- days

Scolnms <- colnms[substr(colnms,6,6)=="S"]
Rcolnms <- colnms[substr(colnms,6,6)=="R"]
Mcolnms <- colnms[substr(colnms,6,6)=="M"]
M1colnms <- Mcolnms[substr(Mcolnms,7,7)=="1"]
M2colnms <- Mcolnms[substr(Mcolnms,7,7)=="2"]
M3colnms <- Mcolnms[substr(Mcolnms,7,7)=="3"]
M4colnms <- Mcolnms[substr(Mcolnms,7,7)=="4"]
M5colnms <- Mcolnms[substr(Mcolnms,7,7)=="5"]

# Invert areas to total numbers
Totalsm <- Areasm
Totalsm[,Scolnms] <- ASinv(Totalsm[,Scolnms])
Totalsm[,Rcolnms] <- ARinv(Totalsm[,Rcolnms])
Totalsm[,Mcolnms] <- ASinv(Totalsm[,Mcolnms])

# Make the figure 
pdf("MCF7fits.pdf")

  par(mar=c(5,5,4,2),font=2, font.axis=2, font.lab=2,
#      font.main=2,family="Arial")
    font.main=2)
  faxissize <- 1.2
  flabsize <- 1.7
  fnamesize <- 1.7
  cols <- c(mcf7cols["S"],"gray30","gray40","gray50","gray60",mcf7cols["R"])
  plot(S0ave ~ Cells,newgm,col=cols[6],pch=19,type="p",
       xlab="Cells",ylab="Area",cex.lab=1.3,cex.axis=1.3)
  points(S100ave ~ Cells,newgm,col=cols[1],pch=19,type="p")
  points(S90ave ~ Cells,newgm,col=cols[2],pch=19,type="p")
  points(S80ave ~ Cells,newgm,col=cols[3],pch=19,type="p")
  points(S70ave ~ Cells,newgm,col=cols[4],pch=19,type="p")
  points(S60ave ~ Cells,newgm,col=cols[5],pch=19,type="p")
  legend("topleft",c("100% S","90% S","80% S","70% S","60% S","0% S"),
         col=cols,lwd=2)
# Superimpose the MM fits
  mmfun <- function(x) A*x/(K+x)
  for (ipct in 1:length(pcts)) {
    K <- fitdat$K[ipct]
    A <- fitdat$A[ipct]
    lines(mmfun(Cells) ~ Cells,newgm,col=cols[ipct],lwd=2)
  }
  K <- all.fit[1]
  A <- all.fit[2]

  plot(unlist(Totalsm[,Scolnms]),unlist(Sfluorsm[,Scolnms]),
       pch=19,col=mcf7cols["S"],
       xlab="Cell Number",ylab="Fluorescence",cex.lab=1.3,cex.axis=1.3)
  points(unlist(Totalsm[,Rcolnms]),unlist(Rfluorsm[,Rcolnms]),
       pch=19,col=mcf7cols["R"])
  Stmpm.lm <- lm(unlist(Sfluorsm[,Scolnms])~unlist(Totalsm[,Scolnms])-1)
  summary(Stmpm.lm)
  abline(Stmpm.lm,col=mcf7cols["S"])
  Rtmpm.lm <- lm(unlist(Rfluorsm[,Rcolnms])~unlist(Totalsm[,Rcolnms])-1)
  summary(Rtmpm.lm)
  abline(Rtmpm.lm,col=mcf7cols["R"])
  legend("topleft",c("Sensitive","Resistant"),pch=19,col=mcf7cols,cex=1.2)

dev.off()
