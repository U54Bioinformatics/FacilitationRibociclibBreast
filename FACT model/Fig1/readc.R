# unpack the colnms in Sc
Sunpack <- function(colnm) {
             l <- nchar(colnm)
             x <- list(type="",env="",level1=0,level2=0,rep=1)
             x["type"] <- "S"
             x["env"] <- ifelse(substr(colnm,l-1,l-1)=="S","S","M")
             x["level1"] <- as.numeric(substr(colnm,2,4))
             x["level2"] <- as.numeric(substr(colnm,6,8))
             x["rep"] <- as.numeric(substr(colnm,l,l))
             x <- as.data.frame(x)
             x <- cbind(day=days,x)
             return(x)
           }

# unpack the colnms in Rc
Runpack <- function(colnm) {
             l <- nchar(colnm)
             x <- list(type="",env="",level1=0,level2=0,rep=1)
             x["type"] <- "R"
             x["env"] <- ifelse(substr(colnm,l-1,l-1)=="R","R","M")
             x["level1"] <- as.numeric(substr(colnm,2,4))
             x["level2"] <- as.numeric(substr(colnm,6,8))
             x["rep"] <- as.numeric(substr(colnm,l,l))
             x <- as.data.frame(x)
             x <- cbind(day=days,x)
             return(x)
           }

#############################################
############  Estradiol1.csv ################
#############################################

#  001,010,100 tenths of nM

Sc <- read.csv("../DATA/SEstradiol1.csv",header=FALSE)
Rc <- read.csv("../DATA/REstradiol1.csv",header=FALSE)

# Change the first row:
tmp <- c("day",rep("R000E000",9),
               rep("R400E000",9),
               rep("R000E001",9),
               rep("R000E010",9),
               rep("R000E100",9),
               rep("R400E001",9),
               rep("R400E010",9),
               rep("R400E100",9))
tmp2 <- c("",rep(c(rep(".S",3),rep(".M",3),rep(".R",3)),8))
tmp3 <- c("",rep(1:3,24))
colnms <- paste0(tmp,tmp2,tmp3)
names(Sc) <- colnms
names(Rc) <- colnms

days <- c(0, 4, 7, 11, 14, 18, 21)
Sc$day <- days
Rc$day <- days

# Get rid of empty columns in Sc and Rc
means <- apply(Sc,2,mean)
badmeans <- names(means[is.na(means)])
Sc[,badmeans] <- NULL

means <- apply(Rc,2,mean)
badmeans <- names(means[is.na(means)])
Rc[,badmeans] <- NULL

# Put all the Sc data together
allSc <- NULL
for (colnm in names(Sc[,2:ncol(Sc)])) {
  x <- Sunpack(colnm)
  x <- cbind(x,cells=Sc[,colnm])
  allSc <- rbind(allSc,x)
}

# Put all the Rc data together
allRc <- NULL
for (colnm in names(Rc[,2:ncol(Rc)])) {
  x <- Runpack(colnm)
  x <- cbind(x,cells=Rc[,colnm])
  allRc <- rbind(allRc,x)
}

# Put it together
Estradiol1 <- rbind(allSc,allRc)
# Fix the units to nM
Estradiol1$level2 <- Estradiol1$level2/10

# Add the name
Estradiol1 <- cbind(file="Estradiol1",Estradiol1)

#######################################
###########  Raloxifene3.csv ##########
#######################################

#  Raloxifene3  25 nm
days <- c(0, 4, 7, 11, 14, 18, 21)
Sc <- read.csv("../DATA/SRaloxifene3.csv",header=FALSE)
Rc <- read.csv("../DATA/RRaloxifene3.csv",header=FALSE)
Sc[,1] <- days
Rc[,1] <- days

# Fix the names
tmp <- c("day",rep("R000E000",9),
               rep("R200E000",9),
               rep("R000E010",9),
               rep("R000E025",9),
               rep("R000E075",9),
               rep("R200E010",9),
               rep("R200E025",9),
               rep("R200E075",9))
tmp2 <- c("",rep(c(rep(".S",3),rep(".M",3),rep(".R",3)),8))
tmp3 <- c("",rep(1:3,24))
colnms <- paste0(tmp,tmp2,tmp3)
names(Sc) <- colnms
names(Rc) <- colnms

# Get rid of empty columns in Sc and Rc
means <- apply(Sc,2,mean)
badmeans <- names(means[is.na(means)])
Sc[,badmeans] <- NULL

means <- apply(Rc,2,mean)
badmeans <- names(means[is.na(means)])
Rc[,badmeans] <- NULL

# Put all the Sc data together
allSc <- NULL
for (colnm in names(Sc[,2:ncol(Sc)])) {
  x <- Sunpack(colnm)
  x <- cbind(x,cells=Sc[,colnm])
  allSc <- rbind(allSc,x)
}

# Put all the Rc data together
allRc <- NULL
for (colnm in names(Rc[,2:ncol(Rc)])) {
  x <- Runpack(colnm)
  x <- cbind(x,cells=Rc[,colnm])
  allRc <- rbind(allRc,x)
}

# Put it together
Raloxifene3 <- rbind(allSc,allRc)
Raloxifene3 <- cbind(file="Raloxifene3",Raloxifene3)

#######################################
##########  4OH_Tamoxifen.csv #########
#######################################

#  4OH_Tamoxifen  1 and 5 mu
days <- c(0, 4, 7, 11, 14, 18)
Sc <- read.csv("../DATA/S4OH_Tamoxifen.csv",header=FALSE)
Rc <- read.csv("../DATA/R4OH_Tamoxifen.csv",header=FALSE)
Sc[,1] <- days
Rc[,1] <- days

# Fix the names
tmp <- c("day",rep("R000E000",9),
               rep("R200E000",9),
               rep("R000E001",9),
               rep("R000E005",9),
               rep("R200E001",9),
               rep("R200E005",9))
tmp2 <- c("",rep(c(rep(".S",3),rep(".M",3),rep(".R",3)),6))
tmp3 <- c("",rep(1:3,18))
colnms <- paste0(tmp,tmp2,tmp3)
names(Sc) <- colnms
names(Rc) <- colnms

# Get rid of empty columns in Sc and Rc
means <- apply(Sc,2,mean)
badmeans <- names(means[is.na(means)])
Sc[,badmeans] <- NULL

means <- apply(Rc,2,mean)
badmeans <- names(means[is.na(means)])
Rc[,badmeans] <- NULL

# Put all the Sc data together
allSc <- NULL
for (colnm in names(Sc[,2:ncol(Sc)])) {
  x <- Sunpack(colnm)
  x <- cbind(x,cells=Sc[,colnm])
  allSc <- rbind(allSc,x)
}

# Put all the Rc data together
allRc <- NULL
for (colnm in names(Rc[,2:ncol(Rc)])) {
  x <- Runpack(colnm)
  x <- cbind(x,cells=Rc[,colnm])
  allRc <- rbind(allRc,x)
}

# Put it together
FOHTamoxifen <- rbind(allSc,allRc)
FOHTamoxifen <- cbind(file="FOHTamoxifen",FOHTamoxifen)

###############################################
############  Fulvestrant3.csv ################
###############################################
#  001,003,010 uM

# I think the layout is slightly different
# unpack the colnms in Sc
Sunpack <- function(colnm) {
             l <- nchar(colnm)
             x <- list(type="",env="",level1=0,level2=0,rep=1)
             x["type"] <- "S"
             x["env"] <- ifelse(substr(colnm,l-1,l-1)=="S","S","M")
             x["level1"] <- as.numeric(substr(colnm,2,4))
             x["level2"] <- as.numeric(substr(colnm,6,8))
             x["rep"] <- as.numeric(substr(colnm,l,l))
             x <- as.data.frame(x)
             x <- cbind(day=days,x)
             return(x)
           }

# unpack the colnms in Rc
Runpack <- function(colnm) {
             l <- nchar(colnm)
             x <- list(type="",env="",level1=0,level2=0,rep=1)
             x["type"] <- "R"
             x["env"] <- ifelse(substr(colnm,l-1,l-1)=="R","R","M")
             x["level1"] <- as.numeric(substr(colnm,2,4))
             x["level2"] <- as.numeric(substr(colnm,6,8))
             x["rep"] <- as.numeric(substr(colnm,l,l))
             x <- as.data.frame(x)
             x <- cbind(day=days,x)
             return(x)
           }

Sc <- read.csv("../DATA/SFulvestrant3.csv",header=FALSE)
Rc <- read.csv("../DATA/RFulvestrant3.csv",header=FALSE)

# Make the names
tmp <- c("day",rep("R000R000",9),
               rep("R400R000",9),
               rep("R000R001",9),
               rep("R000R003",9),
               rep("R000R010",9),
               rep("R400R001",9),
               rep("R400R003",9),
               rep("R400R010",9))
tmp2 <- c("",rep(c(rep(".S",3),rep(".M",3),rep(".R",3)),8))
tmp3 <- c("",rep(1:3,24))
colnms <- paste0(tmp,tmp2,tmp3)
names(Sc) <- colnms
names(Rc) <- colnms

days <- c(0, 4, 7, 11, 14, 18)
Sc$day <- days
Rc$day <- days

# Get rid of empty columns in Sc and Rc
means <- apply(Sc,2,mean)
badmeans <- names(means[is.na(means)])
Sc[,badmeans] <- NULL

means <- apply(Rc,2,mean)
badmeans <- names(means[is.na(means)])
Rc[,badmeans] <- NULL

# Put all the Sc data together
allSc <- NULL
for (colnm in names(Sc[,2:ncol(Sc)])) {
  x <- Sunpack(colnm)
  x <- cbind(x,cells=Sc[,colnm])
  allSc <- rbind(allSc,x)
}

# Put all the Rc data together
allRc <- NULL
for (colnm in names(Rc[,2:ncol(Rc)])) {
  x <- Runpack(colnm)
  x <- cbind(x,cells=Rc[,colnm])
  allRc <- rbind(allRc,x)
}

# Put it together
Fulvestrant3 <- rbind(allSc,allRc)
Fulvestrant3 <- cbind(file="Fulvestrant3",Fulvestrant3)

# There's a mistake in here
Fulvestrant3[457,] #14467179
# Fix it by dividing by 10000
Fulvestrant3[457,"cells"] <- Fulvestrant3[457,"cells"]/10000

#############################################
############  Letrozole_4.0.csv ################
#############################################

#  5,1,0.2 tenths of uM

Sc <- read.csv("../DATA/SLetrozole_4.0.csv",header=FALSE)
Rc <- read.csv("../DATA/RLetrozole_4.0.csv",header=FALSE)

# Change the first row:
tmp <- c("day",rep("R000L000",9),
               rep("R200L000",9),
               rep("R000L050",9),
               rep("R000L010",9),
               rep("R000L002",9),
               rep("R200L050",9),
               rep("R200L010",9),
               rep("R200L002",9))
tmp2 <- c("",rep(c(rep(".S",3),rep(".M",3),rep(".R",3)),8))
tmp3 <- c("",rep(1:3,24))
colnms <- paste0(tmp,tmp2,tmp3)
names(Sc) <- colnms
names(Rc) <- colnms

days <- c(0, 4, 7, 11, 14, 18)
Sc$day <- days
Rc$day <- days

# Get rid of empty columns in Sc and Rc
means <- apply(Sc,2,mean)
badmeans <- names(means[is.na(means)])
Sc[,badmeans] <- NULL

means <- apply(Rc,2,mean)
badmeans <- names(means[is.na(means)])
Rc[,badmeans] <- NULL

# Put all the Sc data together
allSc <- NULL
for (colnm in names(Sc[,2:ncol(Sc)])) {
  x <- Sunpack(colnm)
  x <- cbind(x,cells=Sc[,colnm])
  allSc <- rbind(allSc,x)
}

# Put all the Rc data together
allRc <- NULL
for (colnm in names(Rc[,2:ncol(Rc)])) {
  x <- Runpack(colnm)
  x <- cbind(x,cells=Rc[,colnm])
  allRc <- rbind(allRc,x)
}

# Put it together
Letrozole4 <- rbind(allSc,allRc)

# Fix the units to uM
Letrozole4$level2 <- Letrozole4$level2/10

# Add the name
Letrozole4 <- cbind(file="Letrozole4",Letrozole4)


################################
########## Clean up ############
################################

rm(Rc)
rm(allRc)
rm(allSc)
rm(badmeans)
rm(colnm)
rm(colnms)
rm(days)
rm(means)
rm(Sc)
rm(tmp)
rm(tmp2)
rm(tmp3)
rm(x)
rm(Sunpack)
rm(Runpack)

