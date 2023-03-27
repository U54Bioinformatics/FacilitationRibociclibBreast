EresultsM <- read.csv("../DATA/EresultsM.csv")
# This one has all the info, the others will get added on
# I did a fair amount of hand editing
EresultsF <- read.csv("../DATA/EresultsF.csv")
EresultsFGF <- read.csv("../DATA/EresultsFGF.csv")

tmp <- EresultsM
tmp$Estrone <- EresultsF$V1
tmp$Estradiol <- EresultsF$V2
tmp$Sample <- paste0("F",1:18)
EresultsF <- tmp

tmp <- EresultsM
tmp$Estrone <- EresultsFGF$V1
tmp$Estradiol <- EresultsFGF$V2
tmp$Sample <- paste0("FGF",1:18)
EresultsFGF <- tmp

EresultsF$rep <- rep(1:3)
EresultsFGF$rep <- rep(1:3)
EresultsM$rep <- rep(1:3)

# Read in the count data

source("Ereadm.R")
dy <- 21   
EcountsF <- subset(EcountsF,day==dy)
EcountsFGF <- subset(EcountsFGF,day==dy)
EcountsM <- subset(EcountsM,day==dy)

# Now put the number of cells of each type on day dy
# Need to add columns for S and R cells
EresultsF$S <- EresultsFGF$S <- EresultsM$S <- 0
EresultsF$R <- EresultsFGF$R <- EresultsM$R <-  0

# Now we need to merge on the number of cells
EresultsF$S[EresultsF$env %in% c("S","M")] <- 
  EcountsF$cells[EcountsF$type =="S"]
EresultsF$R[EresultsF$env %in% c("R","M")] <- 
  EcountsF$cells[EcountsF$type =="R"]

EresultsFGF$S[EresultsFGF$env %in% c("S","M")] <- 
  EcountsFGF$cells[EcountsFGF$type =="S"]
EresultsFGF$R[EresultsFGF$env %in% c("R","M")] <- 
  EcountsFGF$cells[EcountsFGF$type =="R"]

EresultsM$S[EresultsM$env %in% c("S","M")] <- 
  EcountsM$cells[EcountsM$type =="S"]
EresultsM$R[EresultsM$env %in% c("R","M")] <- 
  EcountsM$cells[EcountsM$type =="R"]

# Now bind it all together
EresultsF <- cbind(file="F",EresultsF)
EresultsFGF <- cbind(file="FGF",EresultsFGF)
EresultsM <- cbind(file="M",EresultsM)

Eresults <- rbind(EresultsF,EresultsFGF,EresultsM)

# Fix the 0's
oldLLOD <- 0.0
LLOD <- 62.5
Eresults$Estradiol[Eresults$Estradiol==oldLLOD/2] <- LLOD/2

