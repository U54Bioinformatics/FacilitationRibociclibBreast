source("readcharcoal.R")  #Read in the data
source("functions.R")
camacols <- c("chartreuse2","cyan3")
names(camacols) <- c("S","R")
usefiles <- c("Charcoal")

# Create files to save the results
anms <- c("level1","level2","type","env")
for (nm in usefiles) {
  tmp <- get(nm)
  tmp$ilvl1 <- (rank(tmp$level1,ties.method="min")-1)/(nrow(tmp)/length(unique(tmp$level1)))+1
  tmp$ilvl2 <- (rank(tmp$level2,ties.method="min")-1)/(nrow(tmp)/length(unique(tmp$level2)))+1
# Reorder columns
  tmp <- tmp[,c("file","day","level1","level2","type","env","ilvl1","ilvl2","rep","cells")]
  assign(nm,tmp)
# Set up to store parameters: Add columns for r, K, alpha
  tmp <- tmp[duplicated(tmp[,anms])==FALSE,]
  tmp <- tmp[,c("file",anms,"ilvl1","ilvl2")]
# Now add on the data
  tmp[,c("r","K","alpha","alphaval")] <- 0
# Rename
  assign(paste0(nm,".rK"),tmp)
}
rm(anms)

# Find the parameter fits
plotm <- 0
source("step125charcoal.R")

# Now make the file
pdf("suppfig7.pdf")
par(mar=c(5,5,4,2),font=2, font.axis=2, font.lab=2,
#    font.main=2,family="Arial")
    font.main=2)
faxissize <- 1.2
flabsize <- 1.7
fnamesize <- 1.7
fmainsize <- 1.7

# Graph growth curves with charcoal
nm <- "Charcoal"
tfile <- get(nm)
tfile.rK <- get(paste0(nm,".rK"))
plot(cells ~ day,tfile,col="white",
     xlab="Day",ylab="Cell Number",cex.lab=flabsize,cex.axis=faxissize)
tmpfile <- subset(tfile,env != "M" & level1==0 & level2==1)
tmpfile.rK <- subset(tfile.rK,env != "M" & level1==0 & level2==1)
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
legend(-1.0,11500,c("Sensitive","Resistant"),lwd=4,bty="n",
     col=camacols,cex=1.3,seg.len=4)
legend(-0.5,9500,c("Added estradiol","No added estradiol"),bty="n",
     col=1,lty=c(1,2),lwd=c(3,2),pch=19,pt.cex=c(1.5,1.0),cex=1.3,seg.len=4)
dev.off()

