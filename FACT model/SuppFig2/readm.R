##################################
############  MCF7.csv ###########
##################################

# This has a different structure, and we'll just make the ones with
# different proportions into different replicates
# unpack the colnms in Sc
Sunpack <- function(colnm) {
             l <- nchar(colnm)
             x <- list(type="",env="",level1=0,level2=0,rep=1)
             x["type"] <- "S"
             x["env"] <- substr(colnm,l-2,l-1)
             x["level1"] <- as.numeric(substr(colnm,2,4))
             x["level2"] <- 0
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
             x["env"] <- substr(colnm,l-2,l-1)
             x["level1"] <- as.numeric(substr(colnm,2,4))
             x["level2"] <- 0
             x["rep"] <- as.numeric(substr(colnm,l,l))
             x <- as.data.frame(x)
             x <- cbind(day=days,x)
             return(x)
           }

#  MCF7  0.2 um
days <- c(0, 4, 7, 11, 14, 18, 21)
Sc <- read.csv("../DATA/SMCF7.csv",header=FALSE)
Rc <- read.csv("../DATA/RMCF7.csv",header=FALSE)
Sc[,1] <- days
Rc[,1] <- days

# Fix the names: This has 4 ribo doses, and 7 initial conditions
tmp <- c("day",rep(c(rep("R000",3),
               rep("R015",3),
               rep("R024",3),
               rep("R050",3)),7))
tmp2 <- c("",c(rep(".S0",12),rep(".R0",12),rep(paste0(".M",1:5),each=12)))
tmp3 <- c("",rep(1:3,28))
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
MCF7 <- rbind(allSc,allRc)
MCF7 <- cbind(file="MCF7",MCF7)

# Divide the huge ones by 1e5
MCF7$cells[MCF7$env == "M5"] <- MCF7$cells[MCF7$env == "M5"]/1e5

##################################
############  newMCF7.csv ########
##################################

#  newMCF7  
days <- c(0, 4, 7, 11, 14, 18)
Sc <- read.csv("../DATA/SnewMCF7.csv",header=FALSE)
Rc <- read.csv("../DATA/RnewMCF7.csv",header=FALSE)
Areasc <- read.csv("../DATA/AreanewMCF7.csv",header=FALSE)
Fractionsc <- read.csv("../DATA/FractionnewMCF7.csv",header=FALSE)
Sc[,1] <- days
Rc[,1] <- days
Areasc[,1] <- days
Fractionsc[,1] <- days

# Fix the names: This has 4 ribo doses, and 7 initial conditions
tmp <- c("day",rep(c(rep("R000",3),
               rep("R015",3),
               rep("R024",3),
               rep("R050",3)),7))
tmp2 <- c("",c(rep(".S0",12),rep(".R0",12),rep(paste0(".M",1:5),each=12)))
tmp3 <- c("",rep(1:3,28))
colnms <- paste0(tmp,tmp2,tmp3)
names(Sc) <- colnms
names(Rc) <- colnms
names(Areasc) <- colnms
names(Fractionsc) <- colnms

Scolnms <- colnms[substr(colnms,6,6)=="S"]
Rcolnms <- colnms[substr(colnms,6,6)=="R"]
Mcolnms <- colnms[substr(colnms,6,6)=="M"]
# Replace the fractions with 1.0 for sensitive and 0.0 with resistant
Fractionsc[,colnms[substr(colnms,6,6)=="S"]] <- 1.0
Fractionsc[,colnms[substr(colnms,6,6)=="R"]] <- 0.0

# Use this and ARinv and AR made in growthread.R to convert Areas to
# total cell numbers
Totalsc <- Areasc
Totalsc[,Scolnms] <- ASinv(Totalsc[,Scolnms])
Totalsc[,Rcolnms] <- ARinv(Totalsc[,Rcolnms])
Totalsc[,Mcolnms] <- ASinv(Totalsc[,Mcolnms])

Scnew <- Totalsc*Fractionsc
Rcnew <- Totalsc*(1-Fractionsc)
Scnew$day <- days
Rcnew$day <- days

# Get rid of empty columns in Sc and Rc
Sc[,Rcolnms] <- NULL
Rc[,Scolnms] <- NULL

# Get rid of zeros in Scnew and Rcnew
Scnew[,Rcolnms] <- NULL
Rcnew[,Scolnms] <- NULL

# Put all the Sc and Scnew data together
allSc <- NULL
for (colnm in names(Sc[,2:ncol(Sc)])) {
  x <- Sunpack(colnm)
  x <- cbind(x,oldcells=Sc[,colnm])
  x <- cbind(x,cells=Scnew[,colnm])
  allSc <- rbind(allSc,x)
}

# Put all the Rc data together
allRc <- NULL
for (colnm in names(Rc[,2:ncol(Rc)])) {
  x <- Runpack(colnm)
  x <- cbind(x,oldcells=Rc[,colnm])
  x <- cbind(x,cells=Rcnew[,colnm])
  allRc <- rbind(allRc,x)
}

# oldMCF7 <- newMCF7
# Put it together
newMCF7 <- rbind(allSc,allRc)
newMCF7 <- cbind(file="newMCF7",newMCF7)

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

# Create file to save the results
anms <- c("level1","level2","type","env")
nm <- "newMCF7"
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

rm(anms)

