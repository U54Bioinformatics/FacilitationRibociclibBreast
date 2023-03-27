################################################################
############## Get some files ready ############################
################################################################
cama.rK <- NULL
for (nm in sort(camafiles)) cama.rK <- rbind(cama.rK,get(paste0(nm,".rK")))
camadat <- NULL
for (nm in sort(camafiles)) camadat <- rbind(camadat,get(nm))

rKnms <- c("r","K")
rKnmsS <- paste0(rKnms,"S") 
rKnmsR <- paste0(rKnms,"R") 

# Make a file with all the monoculture data in wide form 
gdays <- c(0,4,7,11,14,18)
gday <- gdays[1]
camacompdat <- subset(camadat,env != "M" & day==gday)
camacompdat <- camacompdat[,c("file","level1","level2","type","rep","cells")]
names(camacompdat)[names(camacompdat)=="cells"] <- paste0("cells",gday)

for (iday in 2:length(gdays)) {
  gday <- gdays[iday]
  tmp <- subset(camadat,env != "M" & day==gday)
  camacompdat[,paste0("cells",gday)] <- tmp$cells
}

# All the data are in camadat.  We will make a full file that has the
# preds and the null model for all cases with two interesting things going on:
#  logistic fits : monoculture (env=="S"" or "R") with one treatment
#  LV fits : coculture (env=="M") with no treatment
#  facilitation: coculture (env=="M") and one of level1>0 or level2 > 0
#  synergy: monoculture (env!="M") and level1>0 and level2 > 0
#  facilitation modification: coculture and level1>0 and level2 > 0

# Restrict to days 4-18 to make fulldat from camadat
# Loop through cama.rK.  Create fits and null model and merge on
days <- c(4,7,11,14,18)
fulldat <- subset(camadat,day %in% days)
fulldat$pred <- 0
fulldat$null <- 0
fulldat <- cbind(label="a",fulldat)

for (i in 1:nrow(cama.rK)) {
  frows <- ((15*(i-1)+1):(15*i))
  tmp.rK <- cama.rK[i,]
  tmp <- subset(fulldat,fulldat$file==tmp.rK$file &
                        fulldat$level1==tmp.rK$level1 &
                        fulldat$level2==tmp.rK$level2 &
                        fulldat$type==tmp.rK$type &
                        fulldat$env==tmp.rK$env)
  if ((tmp.rK$level1==0 | tmp.rK$level2==0) & tmp.rK$env!="M") {   #use logistic fit
    tmp$label <- "control"
    pred <- lfun(mean(tmp$cells[tmp$day==4]),tmp.rK$K,tmp.rK$r,days)
    tmp$pred <- pred
    fulldat[frows,] <- tmp
  }
  if (tmp.rK$level1==0 & tmp.rK$level2==0 & tmp.rK$env=="M") {    #use LV fit, redundant
    tmp$label <- "compete"
    competeS <- subset(fulldat,file==tmp.rK$file & env=="M" & 
                       type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    competeR <- subset(fulldat,file==tmp.rK$file & env=="M" & 
                       type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    competeS.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                       type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    competeR.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                       type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
# Pull out the parameters 
    competeparms <- c(rS=competeS.rK$r,
                      KS=competeS.rK$K,
                      rR=competeR.rK$r,
                      KR=competeR.rK$K,
                      alphaRS=competeR.rK$alpha,
                      alphaSR=competeS.rK$alpha)
    init <- c(mean(competeS$cells[competeS$day==4]),
              mean(competeR$cells[competeR$day==4]))
    names(init) <- c("S","R")
    LVout <- as.data.frame(ode(y=init,parms=competeparms,times=days,func=LVmodel))
    if (tmp.rK$type=="S") tmp$pred <- LVout$S
    if (tmp.rK$type=="R") tmp$pred <- LVout$R
    fulldat[frows,] <- tmp
  }
  if (tmp.rK$level1>0 & tmp.rK$level2==0 & tmp.rK$env=="M") { # facilitation, ribo
    tmp$label <- "ribofacil"
    facilS <- subset(fulldat,file==tmp.rK$file & env=="M" & 
                       type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    facilR <- subset(fulldat,file==tmp.rK$file & env=="M" & 
                       type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    controlS.rK <- subset(cama.rK,file==tmp.rK$file & env=="S" & 
                          type=="S" & level1==0 & level2==tmp.rK$level2)
    treatS.rK <- subset(cama.rK,file==tmp.rK$file & env=="S" & 
                        type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    competeS.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                          type=="S" & level1==0 & level2==tmp.rK$level2)
    facilS.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                        type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    controlR.rK <- subset(cama.rK,file==tmp.rK$file & env=="R" & 
                          type=="R" & level1==0 & level2==tmp.rK$level2)
    treatR.rK <- subset(cama.rK,file==tmp.rK$file & env=="R" & 
                        type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    competeR.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                          type=="R" & level1==0 & level2==tmp.rK$level2)
    facilR.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                        type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    facilparms <- c(rS=facilS.rK$r,
                    KS=facilS.rK$K,
                    rR=facilR.rK$r,
                    KR=facilR.rK$K,
                    alphaRS=facilR.rK$alpha,
                    alphaSR=facilS.rK$alpha)
    init <- c(mean(facilS$cells[facilS$day==4]),
              mean(facilR$cells[facilR$day==4]))
    names(init) <- c("S","R")
    LVout <- as.data.frame(ode(y=init,parms=facilparms,times=days,func=LVmodel))
    if (tmp.rK$type=="S") tmp$pred <- LVout$S
    if (tmp.rK$type=="R") tmp$pred <- LVout$R
# Also need the null model for this one, which requires more parameters
    nullparms <- facilparms
    nullparms["rS"] <- (treatS.rK$r)*(competeS.rK$r)/controlS.rK$r
    nullparms["KS"] <- (treatS.rK$K)*(competeS.rK$K)/controlS.rK$K
    nullparms["rR"] <- (treatR.rK$r)*(competeR.rK$r)/controlR.rK$r
    nullparms["KR"] <- (treatR.rK$K)*(competeR.rK$K)/controlR.rK$K
    LVout <- as.data.frame(ode(y=init, parms=nullparms,times=days, func=LVmodel))
    if (tmp.rK$type=="S") tmp$null <- LVout$S
    if (tmp.rK$type=="R") tmp$null <- LVout$R
    fulldat[frows,] <- tmp
  }
  if (tmp.rK$level1==0 & tmp.rK$level2>0 & tmp.rK$env=="M") { # facilitation, other
    tmp$label <- "modfacil"
    facilS <- subset(fulldat,file==tmp.rK$file & env=="M" & 
                       type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    facilR <- subset(fulldat,file==tmp.rK$file & env=="M" & 
                       type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    controlS.rK <- subset(cama.rK,file==tmp.rK$file & env=="S" & 
                          type=="S" & level1==0 & level2==0)
    treatS.rK <- subset(cama.rK,file==tmp.rK$file & env=="S" & 
                        type=="S" & level1==0 & level2==tmp.rK$level2)
    competeS.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                          type=="S" & level1==0 & level2==0)
    facilS.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                        type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    controlR.rK <- subset(cama.rK,file==tmp.rK$file & env=="R" & 
                          type=="R" & level1==0 & level2==0)
    treatR.rK <- subset(cama.rK,file==tmp.rK$file & env=="R" & 
                        type=="R" & level1==0 & level2==tmp.rK$level2)
    competeR.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                          type=="R" & level1==0 & level2==0)
    facilR.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                        type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    facilparms <- c(rS=facilS.rK$r,
                    KS=facilS.rK$K,
                    rR=facilR.rK$r,
                    KR=facilR.rK$K,
                    alphaRS=facilR.rK$alpha,
                    alphaSR=facilS.rK$alpha)
    init <- c(mean(facilS$cells[facilS$day==4]),
              mean(facilR$cells[facilR$day==4]))
    names(init) <- c("S","R")
    LVout <- as.data.frame(ode(y=init,parms=facilparms,times=days,func=LVmodel))
    if (tmp.rK$type=="S") tmp$pred <- LVout$S
    if (tmp.rK$type=="R") tmp$pred <- LVout$R
# Also need the null model for this one, which requires more parameters
    nullparms <- facilparms
    nullparms["rS"] <- (treatS.rK$r)*(competeS.rK$r)/controlS.rK$r
    nullparms["KS"] <- (treatS.rK$K)*(competeS.rK$K)/controlS.rK$K
    nullparms["rR"] <- (treatR.rK$r)*(competeR.rK$r)/controlR.rK$r
    nullparms["KR"] <- (treatR.rK$K)*(competeR.rK$K)/controlR.rK$K
    LVout <- as.data.frame(ode(y=init, parms=nullparms,times=days, func=LVmodel))
    if (tmp.rK$type=="S") tmp$null <- LVout$S
    if (tmp.rK$type=="R") tmp$null <- LVout$R
    fulldat[frows,] <- tmp
  }

  if (tmp.rK$level1>0 & tmp.rK$level2>0 & tmp.rK$env!="M") {  # synergy
    tmp$label <- "synergy"
# Can use tmp for dynamics with both and tmp.rK for parameters
    control.rK <- subset(cama.rK,file==tmp.rK$file & env==tmp.rK$env & 
                                 type==tmp.rK$type & level1==0 & level2==0)
    ribo.rK <- subset(cama.rK,file==tmp.rK$file & env==tmp.rK$env & 
                              type==tmp.rK$type & level1==tmp.rK$level1 & level2==0)
    treat.rK <- subset(cama.rK,file==tmp.rK$file & env==tmp.rK$env & 
                               type==tmp.rK$type & level1==0 & level2==tmp.rK$level2)
# Can use tmp.rK for parameters

    pred <- lfun(mean(tmp$cells[tmp$day==4]),tmp.rK$K,tmp.rK$r,days)
    tmp$pred <- pred

# Make the null models to compare
    nullparms <- tmp.rK
    nullparms["r"] <- (treat.rK$r)*(ribo.rK$r)/control.rK$r
    nullparms["K"] <- (treat.rK$K)*(ribo.rK$K)/control.rK$K
    null <- lfun(mean(tmp$cells[tmp$day==4]),nullparms$K,nullparms$r,days)
    tmp$null <- null
    fulldat[frows,] <- tmp
  }
  if (tmp.rK$level1>0 & tmp.rK$level2>0 & tmp.rK$env=="M") { # facilitation modification
    tmp$label <- "facilmod"
# Same as ribo facilitation but with a background of the modifier
    facilS <- subset(fulldat,file==tmp.rK$file & env=="M" & 
                       type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    facilR <- subset(fulldat,file==tmp.rK$file & env=="M" & 
                       type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    controlS.rK <- subset(cama.rK,file==tmp.rK$file & env=="S" & 
                          type=="S" & level1==0 & level2==tmp.rK$level2)
    treatS.rK <- subset(cama.rK,file==tmp.rK$file & env=="S" & 
                        type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    competeS.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                          type=="S" & level1==0 & level2==tmp.rK$level2)
    facilS.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                        type=="S" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    controlR.rK <- subset(cama.rK,file==tmp.rK$file & env=="R" & 
                          type=="R" & level1==0 & level2==tmp.rK$level2)
    treatR.rK <- subset(cama.rK,file==tmp.rK$file & env=="R" & 
                        type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)
    competeR.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                          type=="R" & level1==0 & level2==tmp.rK$level2)
    facilR.rK <- subset(cama.rK,file==tmp.rK$file & env=="M" & 
                        type=="R" & level1==tmp.rK$level1 & level2==tmp.rK$level2)

    facilparms <- c(rS=facilS.rK$r,
                    KS=facilS.rK$K,
                    rR=facilR.rK$r,
                    KR=facilR.rK$K,
                    alphaRS=facilR.rK$alpha,
                    alphaSR=facilS.rK$alpha)
    init <- c(mean(facilS$cells[facilS$day==4]),
              mean(facilR$cells[facilR$day==4]))
    names(init) <- c("S","R")
    LVout <- as.data.frame(ode(y=init,parms=facilparms,times=days,func=LVmodel))
    if (tmp.rK$type=="S") tmp$pred <- LVout$S
    if (tmp.rK$type=="R") tmp$pred <- LVout$R
# Also need the null model for this one, which requires more parameters
    nullparms <- facilparms
    nullparms["rS"] <- (treatS.rK$r)*(competeS.rK$r)/controlS.rK$r
    nullparms["KS"] <- (treatS.rK$K)*(competeS.rK$K)/controlS.rK$K
    nullparms["rR"] <- (treatR.rK$r)*(competeR.rK$r)/controlR.rK$r
    nullparms["KR"] <- (treatR.rK$K)*(competeR.rK$K)/controlR.rK$K
    LVout <- as.data.frame(ode(y=init, parms=nullparms,times=days, func=LVmodel))
    if (tmp.rK$type=="S") tmp$null <- LVout$S
    if (tmp.rK$type=="R") tmp$null <- LVout$R
    fulldat[frows,] <- tmp
  }
}

# Make the summary files
gdays <- c(11,14,18)
camafacil <- subset(fulldat,label %in% c("ribofacil","modfacil") & day %in% gdays)
# Need good labels
camafacil$goodlabel <- camafacil$file
camafacil$goodlabel[camafacil$file=="Estradiol1"] <- "Estradiol"
camafacil$goodlabel[camafacil$file=="Fulvestrant3"] <- "Fulvestrant"
camafacil$goodlabel[camafacil$file=="Letrozole4"] <- "Letrozole"
camafacil$goodlabel[camafacil$file=="Raloxifene3"] <- "Raloxifene"
camafacil$goodlabel[camafacil$file=="FOHTamoxifen"] <- "Tamoxifen"
camafacil$goodlabel[camafacil$label=="ribofacil"] <- "Ribociclib"
camafacil$goodlabel <- factor(camafacil$goodlabel,
          levels=c("Ribociclib","Estradiol","Fulvestrant","Letrozole","Raloxifene","Tamoxifen"))
camafacil$type <- factor(camafacil$type,levels=c("S","R"))

# Now do the same for facilitation modification
camafacilmod <- subset(fulldat,label %in% c("ribofacil","facilmod") & day %in% gdays)
# Need good labels
camafacilmod$goodlabel <- camafacilmod$file
camafacilmod$goodlabel[camafacilmod$file=="Estradiol1"] <- "Estradiol"
camafacilmod$goodlabel[camafacilmod$file=="Fulvestrant3"] <- "Fulvestrant"
camafacilmod$goodlabel[camafacilmod$file=="Letrozole4"] <- "Letrozole"
camafacilmod$goodlabel[camafacilmod$file=="Raloxifene3"] <- "Raloxifene"
camafacilmod$goodlabel[camafacilmod$file=="FOHTamoxifen"] <- "Tamoxifen"
camafacilmod$goodlabel[camafacilmod$label=="ribofacil"] <- "Ribociclib"
camafacilmod$goodlabel <- factor(camafacilmod$goodlabel,
          levels=c("Ribociclib","Estradiol","Fulvestrant","Letrozole","Raloxifene","Tamoxifen"))
camafacilmod$type <- factor(camafacilmod$type,levels=c("S","R"))
# Get rid of the negative values in null
camafacilmod$null[camafacilmod$null < 1] <- 1
subset(camafacilmod,log(cells/null) < -4)
# This is all the ones with level2==10.  How about we get rid of these?
camafacilmod <- subset(camafacilmod,file != "Fulvestrant3" | level2 != 10)

