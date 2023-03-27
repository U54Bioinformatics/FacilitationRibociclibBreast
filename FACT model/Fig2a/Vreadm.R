# Read in the data
VT <- read.csv("../DATA/Vince_Treated.csv")
VU <- read.csv("../DATA/Vince_Untreated.csv")
STnms <- paste0("ST",1:3)
MTnms <- paste0("MT",1:3)
RTnms <- paste0("RT",1:3)
SUnms <- paste0("SU",1:3)
MUnms <- paste0("MU",1:3)
RUnms <- paste0("RU",1:3)
names(VT) <- c("day",STnms,MTnms,RTnms)
names(VU) <- c("day",SUnms,MUnms,RUnms)

# Create the averages
VT$STmean <- apply(VT[,STnms],1,mean)
VT$MTmean <- apply(VT[,MTnms],1,mean)
VT$RTmean <- apply(VT[,RTnms],1,mean)

VU$SUmean <- apply(VU[,SUnms],1,mean)
VU$MUmean <- apply(VU[,MUnms],1,mean)
VU$RUmean <- apply(VU[,RUnms],1,mean)

# Looks pretty good.  Now figure out how to plot these sothey are
# visible.

# Actually, this was dumb.  Put into a long format
VT <- subset(VT,day < 21)
VU <- subset(VU,day < 21)

VTlong <- data.frame(day=rep(VT$day,9),env=rep(c("S","M","R"),each=18),
                     treat="T",rep=rep(1:3,each=6),
                     cells=unlist(VT[,c(STnms,MTnms,RTnms)]))
VUlong <- data.frame(day=rep(VU$day,9),env=rep(c("S","M","R"),each=18),
                     treat="U",rep=rep(1:3,each=6),
                     cells=unlist(VU[,c(SUnms,MUnms,RUnms)]))
Vdat <- rbind(VTlong,VUlong)
