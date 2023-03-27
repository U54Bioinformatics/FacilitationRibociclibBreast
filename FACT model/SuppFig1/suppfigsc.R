filenms <- c("Estradiol1","Fulvestrant3","Raloxifene3","4OH_Tamoxifen")
allareas <- NULL
allyfpnorm2 <- NULL
allcfpnorm2 <- NULL
pdf("CAMAfits.pdf")
  par(mar=c(5,5,4,2),font=2, font.axis=2, font.lab=2,
#      font.main=2,family="Arial")
      font.main=2)
# The curve fits
  a <- -0.5222
  b <- 11041
  k <- -240505
  f2 <- function(x) a*x+b*sqrt(x)+k
  finv <- function(A) ((-b+sqrt(b^2-4*a*(k-A)))/(2*a))^2

  cdat <- read.csv("../DATA/fitting_equations.csv")
  cdat <- cdat[,1:3]
  names(cdat) <- c("cells","area1","area2")
  plot(area1 ~ cells,cdat,pch=19,col=1,
       xlab="Cells",ylab="Area",cex.lab=1.3,cex.axis=1.3)
  points(area2 ~ cells,cdat,pch=8,col=1)
  lines(f2(cells) ~ cells,cdat,lwd=1,lty=2,col=1)

# Show the fluorescence to cell relationship
  for (filenm in filenms) {
    print(filenm)
    areas <- read.csv(paste0("../DATA/",filenm,"_area.csv"))
    fluor <- read.csv(paste0("../DATA/",filenm,"_fluor.csv"))
    areas <- areas[,1:11]
    fluor <- fluor[,1:11]
    Anms <- paste0("A",2:10)
    Snms <- paste0("A",2:4)
    Mnms <- paste0("A",5:7)
    Rnms <- paste0("A",8:10)
    names(areas) <- c("day","time",Anms)
    names(fluor) <- c("day","time",Anms)
    days <- c(0,4,7,10,14,18)
    if (nrow(areas)==7) days <- c(days,21)
    areas$day <- days
    areas$time <- NULL
    fluor$time <- NULL
    
    # Find cell numbers by inverting the function f2 with finv
    Acells <- cbind(day=days,finv(areas[,Anms]))
    
    # Pull out the untreated data
    if (nrow(areas)==6) {
      yfp <- fluor[3:8,]
      cfp <- fluor[12:17,]
      yfpnorm <- fluor[22:27,]
      cfpnorm <- fluor[32:37,]
      yfpnorm2 <- fluor[42:47,]
      cfpnorm2 <- fluor[52:57,]
      yfpprop <- fluor[62:67,]
    } else {
      yfp <- fluor[3:9,]
      cfp <- fluor[13:19,]
      yfpnorm <- fluor[24:30,]
      cfpnorm <- fluor[35:41,]
      yfpnorm2 <- fluor[46:52,]
      cfpnorm2 <- fluor[57:63,]
      yfpprop <- fluor[68:74,]
    }
    
    for (nm in Anms) {
      yfp[,nm] <- as.numeric(yfp[,nm])
      cfp[,nm] <- as.numeric(cfp[,nm])
      yfpnorm[,nm] <- as.numeric(yfpnorm[,nm])
      cfpnorm[,nm] <- as.numeric(cfpnorm[,nm])
      yfpnorm2[,nm] <- as.numeric(yfpnorm2[,nm])
      cfpnorm2[,nm] <- as.numeric(cfpnorm2[,nm])
      yfpprop[,nm] <- as.numeric(yfpprop[,nm])
    }
    yfp[,1] <- days
    cfp[,1] <- days
    yfpnorm[,1] <- days
    cfpnorm[,1] <- days
    yfpnorm2[,1] <- days
    cfpnorm2[,1] <- days
    yfpprop[,1] <- days
  
  # Save in a big file
    allareas <- rbind(allareas,cbind(file=filenm,areas))
    allyfpnorm2 <- rbind(allyfpnorm2,cbind(file=filenm,yfpnorm2))
    allcfpnorm2 <- rbind(allcfpnorm2,cbind(file=filenm,cfpnorm2))
  }
    
  plot(unlist(allyfpnorm2[,Snms]) ~ unlist(allareas[,Snms]),pch=19,col="green",
           xlim=c(0,4e6),xlab="Area",ylab="Normalized fluorescence",
           cex.axis=1.3,cex.lab=1.3)
  points(unlist(allcfpnorm2[,Rnms]) ~ unlist(allareas[,Rnms]),pch=19,col="red")
  abline(lm(unlist(allyfpnorm2[,Snms]) ~ unlist(allareas[,Snms])),col="green")
  abline(lm(unlist(allcfpnorm2[,Rnms]) ~ unlist(allareas[,Rnms])),col="red")
  legend("topleft",c("S","R"),pch=19,col=c("green","red"),cex=1.3)
dev.off()

