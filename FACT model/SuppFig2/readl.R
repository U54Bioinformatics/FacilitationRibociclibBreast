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

##################################
############  LY.csv ########
##################################

#  LY  
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

# plot(unlist(rawFractionsc[1,Mcolnms]),col=1+is.na(match(Mcolnms,M1colnms)))
t.test(unlist(rawFractionsc[1,M1colnms]),unlist(rawFractionsc[1,M2colnms]))
# mean of x mean of y 
#   0.52197   0.48666 
# p = 0.0013
# The M1's are supposed to be 80-20, and the M2's are 70-30.
# S/(S+kR)
# 70/(70+30k) = 0.48666 => 1+3k/7=1/0.4866 => k=2.4613
# 80/(80+20k) = 0.52197 => 1+2k/8=1/0.52197 => k=3.6633
# This looks pretty bad, but might not matter if the initial cell
# numbers are not that well controlled.  However, we still need to know
# the relative fluorescence of the two cell types.
# Idea: We convert the areas of the pure cultures to numbers and then
# see how the fluorescences compare

# Invert areas to total numbers
Totalsc <- Areasc
Totalsc[,Scolnms] <- LYASinv(Totalsc[,Scolnms])
Totalsc[,Rcolnms] <- LYARinv(Totalsc[,Rcolnms])
Totalsc[,Mcolnms] <- LYAMinv(Totalsc[,Mcolnms])

plot(unlist(Totalsc[Totalsc$day==0,Scolnms]),
     unlist(Sfluorsc[Totalsc$day==0,Scolnms]),col=LY2cols["S"],
     xlim=c(0,12000),ylim=c(3e8,4e8))
points(unlist(Totalsc[Totalsc$day==0,Rcolnms]),
     unlist(Rfluorsc[Totalsc$day==0,Rcolnms]),col=LY2cols["R"])
allrats <- (unlist(Rfluorsc[Totalsc$day==0,Rcolnms])/
             unlist(Totalsc[Totalsc$day==0,Rcolnms]))/
           (unlist(Sfluorsc[Totalsc$day==0,Scolnms])/
             unlist(Totalsc[Totalsc$day==0,Scolnms]))

# Ignoring the M's for the moment, let's look at the fluorescence
# plot(unlist(Totalsc[,Scolnms]),unlist(Sfluorsc[,Scolnms]),col=LY2cols["S"])
points(unlist(Totalsc[,Rcolnms]),unlist(Rfluorsc[,Rcolnms]),col=LY2cols["R"])
Stmp.lm <- lm(unlist(Sfluorsc[,Scolnms])~unlist(Totalsc[,Scolnms])-1)
summary(Stmp.lm)
abline(Stmp.lm,col=LY2cols["S"])
Rtmp.lm <- lm(unlist(Rfluorsc[,Rcolnms])~unlist(Totalsc[,Rcolnms])-1)
summary(Rtmp.lm)
abline(Rtmp.lm,col=LY2cols["R"])
Frat <- Rtmp.lm$coef/Stmp.lm$coef  # 1.6774
# Not completely out of line with previous, we assume resistant cells
# make this much more fluorescence. Then we compute the corrected
# fractions.  Then we invert to get total cell numbers and the estimated
# numbers.  The only real problem is that the curves are not in order
# for the mixed cells, but I think it would make sense to just find one
# curve for the mixed ones and use that for the inverse.

Frat <- Rtmp.lm$coef/Stmp.lm$coef  # 1.6774
Frat <- mean(allrats)
Frat <- 3.0
# Make new fractions 
Fractionsc <- Areasc
Fractionsc <- Sfluorsc/(Sfluorsc+Rfluorsc/Frat)
Fractionsc$day <- days

# Summarize the fractions issue
Ftest <- data.frame(M1target=rep(0.8,length(M1colnms)),
          M1raw=unlist(rawFractionsc[rawFractionsc$day==0,M1colnms]),
          M1corr=unlist(Fractionsc[Fractionsc$day==0,M1colnms]),
          M2target=rep(0.7,length(M2colnms)),
          M2raw=unlist(rawFractionsc[rawFractionsc$day==0,M2colnms]),
          M2corr=unlist(Fractionsc[Fractionsc$day==0,M2colnms]))

# plot(unlist(Fractionsc[1,Mcolnms]),col=1+is.na(match(Mcolnms,M1colnms)))
t.test(unlist(Fractionsc[1,M1colnms]),unlist(Fractionsc[1,M2colnms]))

# Separate numbers
Scnew <- Totalsc*Fractionsc
Rcnew <- Totalsc*(1-Fractionsc)
Scnew$day <- days
Rcnew$day <- days

# Get rid of zeros in Scnew and Rcnew
Scnew[,Rcolnms] <- NULL
Rcnew[,Scolnms] <- NULL

# Put all the Scnew data together
allSc <- NULL
for (colnm in names(Scnew[,2:ncol(Scnew)])) {
  x <- Sunpack(colnm)
  x <- cbind(x,cells=Scnew[,colnm])
  allSc <- rbind(allSc,x)
}

# Put all the Rc data together
allRc <- NULL
for (colnm in names(Rcnew[,2:ncol(Rcnew)])) {
  x <- Runpack(colnm)
  x <- cbind(x,cells=Rcnew[,colnm])
  allRc <- rbind(allRc,x)
}

# Put it together
LY2 <- rbind(allSc,allRc)
LY2 <- cbind(file="LY2",LY2)

################################
########## Clean up ############
################################

rm(colnm)
rm(colnms)
rm(days)
rm(Scnew)
rm(Rcnew)
rm(tmp)
rm(tmp2)
rm(tmp3)
rm(x)
rm(Sunpack)
rm(Runpack)

# Create file to save the results
anms <- c("level1","level2","type","env")
nm <- "LY2"
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

