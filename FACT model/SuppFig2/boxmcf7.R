# Need to merge on the null
# Ran mcf7figs.R to create facilS and facilR. And mcf7facilitation
# tmp <- subset(mcf7facilitation,ilvl1==4)
# tnms <- names(facilS)[1:10]
# facilS <- merge(facilS,tmp[,c("day","nullS")],by="day")
# facilS$logSrat <- log(facilS$cells/facilS$nullS)
# names(facilS) <- c(tnms,"null","lograt")
# facilR <- merge(facilR,tmp[,c("day","nullR")],by="day")
# facilR$logRrat <- log(facilR$cells/facilR$nullR)
# names(facilR) <- c(tnms,"null","lograt")
# facilmcf7 <- rbind(facilS,facilR)
# facilmcf7$type <- factor(facilmcf7$type,levels=c("S","R"))
# names(facilmcf7)[1:2] <- names(facilmcf7)[2:1]

pdf("boxmcf7.pdf")
vyrange <- c(-1.5,1.5)
faxissize <- 1.2
flabsize <- 1.3
fnamesize <- 1.7
fmainsize <- 1.3
yfudge <- 0.3
par(mar=c(9,5,4,2))
dy <- 18
    boxplot(lograt ~ as.factor(type),subset(facilmcf7,day==dy),
        col=mcf7cols,
        ylim=vyrange,
        ylab="Strength of facilitation",
        xaxt="n",xlab="",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
    title(main="Facilitation under treatment: MCF7 cells",
         cex.main=fmainsize)
    abline(h=0,lwd=2,lty=2,col="gray")
    axis(1, labels = FALSE)
    labels <- names(mcf7cols)
    text(x =  seq_along(labels), y = par("usr")[3]-yfudge, srt = 45, 
         adj = 1,labels = labels, xpd = TRUE,cex=flabsize)
par(mar=c(5,5,4,2))
dev.off()

