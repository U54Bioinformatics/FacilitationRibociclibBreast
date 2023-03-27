# I edited the csv output of growth_curves_LY2_2.0.xlsx by hand and will
# finish cleaning here
LY2cols <- c("chartreuse2","darkblue")
names(LY2cols) <- c("S","R")

LYraw <- read.csv("../DATA/coculture.csv",header=TRUE)
LYraw$X <- NULL
tmp <- t(LYraw)
tmp <- tmp[2:nrow(tmp),]
tmp <- as.data.frame(tmp)
tmp$Cells <- tmp[,"V4"]
pcts <- seq(from=100,to=0,by=-20)
tmp[,paste0("V",seq(from=4,to=20,by=4))] <- NULL
names(tmp) <- c(
        paste0(rep(paste0("S",pcts),each=3),
               rep(paste0("_",1:3),6)),"Cells")
gcl <- tmp

# The sprawling R cells are glued back together by the S cells, so we
# can use an overall correction based solely on the S cells.

mean2 <- function(x) mean(na.omit(x))
S100nms <- c("S100_1","S100_2","S100_3")
S80nms <- c("S80_1","S80_2","S80_3")
S60nms <- c("S60_1","S60_2","S60_3")
S40nms <- c("S40_1","S40_2","S40_3")
S20nms <- c("S20_1","S20_2","S20_3")
S0nms <- c("S0_1","S0_2","S0_3")

for (i in 1:ncol(gcl)) gcl[,i] <- as.numeric(gcl[,i])

gcl$S100ave <- apply(gcl[,S100nms],1,mean2)
gcl$S80ave <- apply(gcl[,S80nms],1,mean2)
gcl$S60ave <- apply(gcl[,S60nms],1,mean2)
gcl$S40ave <- apply(gcl[,S40nms],1,mean2)
gcl$S20ave <- apply(gcl[,S20nms],1,mean2)
gcl$S0ave <- apply(gcl[,S0nms],1,mean2)

newgcl <- gcl[,c("Cells","S100ave","S80ave","S60ave",
                         "S40ave","S20ave","S0ave")]
  
# Fit each of the sets with a Michaelis-Menten curve
gcells <- gcl$Cells
gfit <- function(x) {
 K <- x[1]
 A <- x[2]
 sum((A*tmp$Cells/(K+tmp$Cells)-tmp$Area)^2)
}

LYfitdat <- data.frame(pct=pcts,K=0,A=0)
for (ipct in 1:length(pcts)) {
  pct <- pcts[ipct]
  gtnms <- paste0("S",pct,"_",1:3)
  tmp <- matrix(unlist(gcl[,gtnms]),3*length(gcells),1)
  tmp <- data.frame(Cells=rep(gcells,3),Area=tmp)
  tmp <- subset(tmp,is.na(tmp$Area)==FALSE)
  Kstart <- mean(gcells)
  Astart <- max(tmp$Area)
  tmp.fit <- optim(c(Kstart,Astart),gfit)
  LYfitdat[ipct,c("K","A")] <-tmp.fit$par
}

# Fit to all the points with S cells
gtnms <- paste0("S",rep(pcts,each=3),"_",1:3)
tmp <- matrix(unlist(gcl[,gtnms]),3*length(pcts)*length(gcells),1)
tmp <- data.frame(Cells=rep(gcells,3*length(pcts)),Area=tmp)
tmp <- subset(tmp,is.na(tmp$Area)==FALSE)
Kstart <- mean(gcells)
Astart <- max(tmp$Area)
all.fit <- optim(c(Kstart,Astart),gfit)$par

cols <- c(LY2cols["S"],"gray30","gray40","gray50","gray60",LY2cols["R"])

# Superimpose the MM fits
mmfun <- function(x) A*x/(K+x)
for (ipct in 1:length(pcts)) {
  K <- LYfitdat$K[ipct]
  A <- LYfitdat$A[ipct]
}
KM <- mean(LYfitdat$K[2:5])
AM <- mean(LYfitdat$A[2:5])
K <- KM
A <- AM

# Now define it as a function
LYAS <- function(C) LYfitdat$A[LYfitdat$pct==100]*C/(LYfitdat$K[LYfitdat$pct==100]+C)
LYAM <- function(C) AM*C/(KM+C)
LYAR <- function(C) LYfitdat$A[LYfitdat$pct==0]*C/(LYfitdat$K[LYfitdat$pct==0]+C)

# Might as well save the inverse
LYASinv <- function(A) (LYfitdat$K[LYfitdat$pct==100]*A)/(LYfitdat$A[LYfitdat$pct==100]-A)
LYAMinv <- function(A) (KM*A)/(AM-A)
LYARinv <- function(A) (LYfitdat$K[LYfitdat$pct==0]*A)/(LYfitdat$A[LYfitdat$pct==0]-A)

# Huzzah, now we can use these in FACT!

# This is for the fluorescence

days <- c(0, 4, 7, 11, 14, 18)
Areasc <- read.csv("../DATA/LYarea.csv",header=FALSE)
Sfluorsc <- read.csv("../DATA/LYSfluor.csv",header=FALSE)
Rfluorsc <- read.csv("../DATA/LYRfluor.csv",header=FALSE)
Areasc[,1] <- days
Sfluorsc[,1] <- days
Rfluorsc[,1] <- days

# Fix the names: This has 2 ribo doses, and 4 initial conditions
tmp <- c("day",c(rep("R000",12),
               rep("R003",12),
               rep("R005",12)))
tmp2 <- c("",rep(c(rep(".S0",3),rep(".R0",3),rep(paste0(".M",1:2),each=3)),3))
tmp3 <- c("",rep(1:3,12))
colnms <- paste0(tmp,tmp2,tmp3)
names(Areasc) <- colnms
names(Sfluorsc) <- colnms
names(Rfluorsc) <- colnms

# Replace the NA's with 0's in the fluorescence data
Sfluorsc[is.na(Sfluorsc)] <- 0
Rfluorsc[is.na(Rfluorsc)] <- 0

# Make fractions and see what happens
rawFractionsc <- Areasc
rawFractionsc <- Sfluorsc/(Sfluorsc+Rfluorsc)
rawFractionsc$day <- days

Scolnms <- colnms[substr(colnms,6,6)=="S"]
Rcolnms <- colnms[substr(colnms,6,6)=="R"]
Mcolnms <- colnms[substr(colnms,6,6)=="M"]
M1colnms <- Mcolnms[substr(Mcolnms,7,7)=="1"]
M2colnms <- Mcolnms[substr(Mcolnms,7,7)=="2"]

# Invert areas to total numbers
Totalsc <- Areasc
Totalsc[,Scolnms] <- LYASinv(Totalsc[,Scolnms])
Totalsc[,Rcolnms] <- LYARinv(Totalsc[,Rcolnms])
Totalsc[,Mcolnms] <- LYAMinv(Totalsc[,Mcolnms])

pdf("LY2fits.pdf")
# The curve fits
  par(mar=c(5,5,4,2),font=2, font.axis=2, font.lab=2,
#      font.main=2,family="Arial")
      font.main=2)
  cols <- c(LY2cols["S"],"gray30","gray40","gray50","gray60",LY2cols["R"])
  plot(S0ave ~ Cells,newgcl,col=cols[6],pch=19,type="p",
       xlab="Cells",ylab="Area",cex.lab=1.3,cex.axis=1.3)
  points(S100ave ~ Cells,newgcl,col=cols[1],pch=19,type="p")
  points(S80ave ~ Cells,newgcl,col=cols[2],pch=19,type="p")
  points(S60ave ~ Cells,newgcl,col=cols[3],pch=19,type="p")
  points(S40ave ~ Cells,newgcl,col=cols[4],pch=19,type="p")
  points(S20ave ~ Cells,newgcl,col=cols[5],pch=19,type="p")
  legend("topleft",paste0(pcts,"% S"),col=cols,lwd=2)

# Superimpose the MM fits
  mmfun <- function(x) A*x/(K+x)
  for (ipct in 1:length(pcts)) {
    K <- LYfitdat$K[ipct]
    A <- LYfitdat$A[ipct]
    lines(mmfun(Cells) ~ Cells,newgcl,col=cols[ipct],lwd=2)
  }
  KM <- mean(LYfitdat$K[2:5])
  AM <- mean(LYfitdat$A[2:5])
  K <- KM
  A <- AM
  lines(mmfun(Cells) ~ Cells,newgcl,col="red",lwd=3,lty=2)

# Show the fluorescence to cell relatioship
  par(mar=c(5,5,4,2),font=2, font.axis=2, font.lab=2,
#      font.main=2,family="Arial")
      font.main=2)
  plot(unlist(Totalsc[,Scolnms]),unlist(Sfluorsc[,Scolnms]),
       pch=19,col=LY2cols["S"],
       xlab="Cell Number",ylab="Fluorescence",cex.lab=1.3,cex.axis=1.3)
  points(unlist(Totalsc[,Rcolnms]),unlist(Rfluorsc[,Rcolnms]),
       pch=19,col=LY2cols["R"])
  Stmp.lm <- lm(unlist(Sfluorsc[,Scolnms])~unlist(Totalsc[,Scolnms])-1)
  abline(Stmp.lm,col=LY2cols["S"],lwd=3)
  Rtmp.lm <- lm(unlist(Rfluorsc[,Rcolnms])~unlist(Totalsc[,Rcolnms])-1)
  abline(Rtmp.lm,col=LY2cols["R"],lwd=3)
  legend("topleft",c("Sensitive","Resistant"),pch=19,col=LY2cols,cex=1.2)

dev.off()

