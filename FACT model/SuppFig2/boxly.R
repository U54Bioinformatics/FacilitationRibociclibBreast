# Need to merge on the null
# Ran ly2figs.R to create facilS and facilR. And LY2facilitation
# tmp <- subset(LY2facilitation,ilvl1==3)
# tnms <- names(facilS)[1:10]
# facilS <- merge(facilS,tmp[,c("day","nullS")],by="day")
# facilS$logSrat <- log(facilS$cells/facilS$nullS)
# names(facilS) <- c(tnms,"null","lograt")
# facilR <- merge(facilR,tmp[,c("day","nullR")],by="day")
# facilR$logRrat <- log(facilR$cells/facilR$nullR)
# names(facilR) <- c(tnms,"null","lograt")
# facilly <- rbind(facilS,facilR)
# facilly$type <- factor(facilly$type,levels=c("S","R"))

pdf("boxLY2.pdf")
vyrange <- c(-1.5,1.5)
faxissize <- 1.2
flabsize <- 1.3
fnamesize <- 1.7
fmainsize <- 1.3
yfudge <- 0.3
par(mar=c(9,5,4,2))
dy <- 18
    boxplot(lograt ~ as.factor(type),subset(facilly,day==dy),
        col=LY2cols,
        ylim=vyrange,
        ylab="Strength of facilitation",
        xaxt="n",xlab="",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
    title(main="Facilitation under treatment: LY2 cells",
         cex.main=fmainsize)
    abline(h=0,lwd=2,lty=2,col="gray")
    axis(1, labels = FALSE)
    labels <- names(LY2cols)
    text(x =  seq_along(labels), y = par("usr")[3]-yfudge, srt = 45, 
         adj = 1,labels = labels, xpd = TRUE,cex=flabsize)
par(mar=c(5,5,4,2))
dev.off()

