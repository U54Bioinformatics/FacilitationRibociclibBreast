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

#######################################
###########  Charcoal.csv ##########
#######################################

days <- c(0, 4, 7, 11, 14, 18)
Sc <- read.csv("../DATA/SCharcoal.csv",header=FALSE)
Rc <- read.csv("../DATA/RCharcoal.csv",header=FALSE)
Sc[,1] <- days
Rc[,1] <- days

# Fix the names
tmp <- c("day",rep("R000E001",9),
               rep("R000E000",9),
               rep("R000F000",9),
               rep("R200E000",9),
               rep("R200F000",9),
               rep("R200G000",9),
               rep("R400E000",9),
               rep("R400F000",9))
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
Charcoal <- rbind(allSc,allRc)
Charcoal <- cbind(file="Charcoal",Charcoal)

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

