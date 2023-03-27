# Make a boxplot for ribo and for the modifiers
# Set some sizes for the graphs
faxissize <- 1.2
flabsize <- 1.3
fnamesize <- 1.7
fmainsize <- 1.3
yfudge <- 0.3

camafacilf <- subset(camafacil,goodlabel %in% c("Ribociclib","Estradiol","Fulvestrant"))
camafacilf$doselabel <- as.character(camafacilf$goodlabel)
# Get rid of highest dose of Fulvestrant
camafacilf$doselabel[camafacilf$goodlabel=="Ribociclib"] <- 
  paste0("Ribociclib ",camafacilf$level1[camafacilf$goodlabel=="Ribociclib"],"nM")
camafacilf$doselabel[camafacilf$goodlabel=="Estradiol"] <- 
  paste0("Estradiol ",camafacilf$level2[camafacilf$goodlabel=="Estradiol"],"nM")
camafacilf$doselabel[camafacilf$goodlabel=="Fulvestrant"] <- 
  paste0("Fulvestrant ",camafacilf$level2[camafacilf$goodlabel=="Fulvestrant"],"nM")
flevels <- c("Ribociclib 200nM","Ribociclib 400nM",
               "Estradiol 0.1nM","Estradiol 1nM","Estradiol 10nM",
               "Fulvestrant 1nM","Fulvestrant 3nM","Fulvestrant 10nM")
camafacilf$doselabel <- factor(camafacilf$doselabel,levels=flevels)

table(camafacilmod$level1,camafacilmod$level2,camafacilmod$goodlabel)
camafacilmodf <- subset(camafacilmod,
                   goodlabel %in% c("Ribociclib","Estradiol","Fulvestrant") & level1 == 400)
camafacilmodf$doselabel <- as.character(camafacilmodf$goodlabel)
# Get rid of highest dose of Fulvestrant
camafacilmodf <- subset(camafacilmodf,goodlabel != "Fulvestrant" | level2 != 10)
camafacilmodf$doselabel[camafacilmodf$goodlabel=="Ribociclib"] <- 
  paste0("Ribociclib ",camafacilmodf$level1[camafacilmodf$goodlabel=="Ribociclib"],"nM")
camafacilmodf$doselabel[camafacilmodf$goodlabel=="Estradiol"] <- 
  paste0("Estradiol ",camafacilmodf$level2[camafacilmodf$goodlabel=="Estradiol"],"nM")
camafacilmodf$doselabel[camafacilmodf$goodlabel=="Fulvestrant"] <- 
  paste0("Fulvestrant ",camafacilmodf$level2[camafacilmodf$goodlabel=="Fulvestrant"],"nM")
modflevels <- c("Ribociclib 400nM",
               "Estradiol 0.1nM","Estradiol 1nM","Estradiol 10nM",
               "Fulvestrant 1nM","Fulvestrant 3nM","Fulvestrant 10nM")
camafacilmodf$doselabel <- factor(camafacilmodf$doselabel,levels=modflevels)

fuselevels <- c("Ribociclib 200nM","Ribociclib 400nM",
                "Estradiol 0.1nM","Fulvestrant 1nM","Fulvestrant 3nM")

modfuselevels <- c("Ribociclib 400nM","Estradiol 0.1nM",
                  "Fulvestrant 1nM","Fulvestrant 3nM")

camafacilf <- subset(camafacilf,doselabel %in% fuselevels)
camafacilf$doselabel <- factor(camafacilf$doselabel,levels=fuselevels)

camafacilmodf <- subset(camafacilmodf,doselabel %in% modfuselevels)
camafacilmodf$doselabel <- factor(camafacilmodf$doselabel,levels=modfuselevels)

pdf("Fig3.pdf")
Syrange <- c(-2,6)
Ryrange <- c(-2,6)
par(mar=c(9,5,4,2))
boxplot(log(cells/null) ~ doselabel,camafacilf,subset=type=="S",
        col=camacols["S"],
        ylim=Syrange,ylab="Strength of facilitation",
        xaxt="n",xlab="",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
title(main="Facilitation of sensitive cells under treatment",cex.main=fmainsize)
abline(h=0,lwd=2,lty=2,col="gray")
axis(1, labels = FALSE)
labels <- sort(unique(camafacilf$doselabel))
text(x =  seq_along(labels), y = par("usr")[3]-yfudge, srt = 45, adj = 1,
      labels = labels, xpd = TRUE,cex=flabsize)

boxplot(log(cells/null) ~ doselabel,camafacilmodf,subset=type=="S",
        col=camacols["S"],
        ylim=Syrange,ylab="Strength of facilitation",
        xaxt="n",xlab="",
        cex.lab=flabsize,cex.axis=faxissize,cex.names=fnamesize)
title(main="Facilitation of sensitive cells under \n combination treatment with ribociclib",
      cex.main=fmainsize)
abline(h=0,lwd=2,lty=2,col="gray")
axis(1, labels = FALSE)
tmplabels <- sort(unique(camafacilmodf$doselabel))
tmplabels <- as.character(tmplabels)
tmplabels[1] <- "No modifier"
labels <- paste0(tmplabels," + \n","Ribociclib 400nM")
text(x =  seq_along(labels), y = par("usr")[3]-yfudge, srt = 45, adj = 1,
      labels = labels, xpd = TRUE,cex=flabsize)

par(mar=c(5,5,4,2))
dev.off()

# Stats on facilitation: This compares everything with Ribociclib 200
# All but Fulvestrant show less facilitation of S cells
summary(lm(log(cells/null) ~ doselabel,camafacilf,subset=type=="S"))
## (Intercept)                  0.783      0.193    4.05  0.00014 ***
## doselabelRibociclib 400nM    1.953      0.306    6.38  1.9e-08 ***
## doselabelEstradiol 0.1nM    -0.883      0.387   -2.28  0.02559 *  
## doselabelFulvestrant 1nM     0.524      0.387    1.35  0.18044    
## doselabelFulvestrant 3nM     1.660      0.387    4.29  5.9e-05 ***

summary(lm(log(cells/null) ~ doselabel,camafacilf,subset=type=="R"))
## (Intercept)                 -0.214      0.076   -2.81   0.0065 ** 
## doselabelRibociclib 400nM   -0.147      0.120   -1.23   0.2247    
## doselabelEstradiol 0.1nM     0.823      0.152    5.42  8.9e-07 ***
## doselabelFulvestrant 1nM    -0.521      0.152   -3.43   0.0011 ** 
## doselabelFulvestrant 3nM    -0.991      0.152   -6.52  1.1e-08 ***

# Now compare with 0
facilsummf <- data.frame(doselabel=as.character(sort(unique(camafacilf$doselabel))),
                        estimateS=0,stderrS=0,pS=0,estimateR=0,stderrR=0,pR=0)
for (i in 1:nrow(facilsummf)) {
  nm <- facilsummf$doselabel[i]
  tmp.t <- with(subset(camafacilf,type=="S" & doselabel==nm),t.test(log(cells/null)))
  facilsummf$estimateS[i] <- tmp.t$estimate
  facilsummf$stderrS[i] <- tmp.t$stderr
  facilsummf$pS[i] <- tmp.t$p.value
  tmp.t <- with(subset(camafacilf,type=="R" & doselabel==nm),t.test(log(cells/null)))
  facilsummf$estimateR[i] <- tmp.t$estimate
  facilsummf$stderrR[i] <- tmp.t$stderr
  facilsummf$pR[i] <- tmp.t$p.value
}
facilsummf
##        doselabel estimateS  stderrS         pS estimateR  stderrR         pR
## Ribociclib 200nM  0.783154 0.064993 3.8231e-12  -0.21363 0.053845 0.00050888
## Ribociclib 400nM  2.736022 0.449230 1.2005e-05  -0.36090 0.080194 0.00031548
##  Estradiol 0.1nM -0.100246 0.034038 1.8563e-02   0.60960 0.164194 0.00593341
##  Fulvestrant 1nM  1.306806 0.069703 6.7659e-08  -0.73437 0.162385 0.00194375
##  Fulvestrant 3nM  2.442962 0.187844 1.1587e-06  -1.20434 0.190286 0.00022561
# See a nice dose response for Ribociclib and Fulvestrant in S cells

# Stats on facilitation modification: This compares everything with Ribociclib
summary(lm(log(cells/null) ~ doselabel,camafacilmodf,subset=type=="S"))
# Facilitation decreased by Estradiol
## (Intercept)                 2.736      0.321    8.52  1.3e-10 ***
## doselabelEstradiol 0.1nM   -2.478      0.556   -4.46  6.3e-05 ***
## doselabelFulvestrant 1nM   -0.459      0.556   -0.82   0.4144    
## doselabelFulvestrant 3nM    1.720      0.556    3.09   0.0036 ** 

summary(lm(log(cells/null) ~ doselabel,camafacilmodf,subset=type=="R"))
## (Intercept)               -0.3609     0.0566   -6.38  1.3e-07 ***
## doselabelEstradiol 0.1nM   0.2590     0.0980    2.64   0.0116 *  
## doselabelFulvestrant 1nM  -0.2874     0.0980   -2.93   0.0055 ** 
## doselabelFulvestrant 3nM  -0.0867     0.0980   -0.88   0.3814    
# Don't really care about the comparison with ribociclib.

# Compare with 0
facilmodsummf <- data.frame(doselabel=as.character(sort(unique(camafacilmodf$doselabel))),
                        estimateS=0,stderrS=0,pS=0,estimateR=0,stderrR=0,pR=0)
for (i in 1:nrow(facilmodsummf)) {
  nm <- facilmodsummf$doselabel[i]
  tmp.t <- with(subset(camafacilmodf,type=="S" & doselabel==nm),t.test(log(cells/null)))
  facilmodsummf$estimateS[i] <- tmp.t$estimate
  facilmodsummf$stderrS[i] <- tmp.t$stderr
  facilmodsummf$pS[i] <- tmp.t$p.value
  tmp.t <- with(subset(camafacilmodf,type=="R" & doselabel==nm),t.test(log(cells/null)))
  facilmodsummf$estimateR[i] <- tmp.t$estimate
  facilmodsummf$stderrR[i] <- tmp.t$stderr
  facilmodsummf$pR[i] <- tmp.t$p.value
}
facilmodsummf
##         doselabel estimateS  stderrS         pS estimateR  stderrR         pR
## Ribociclib 400nM   2.73602 0.449230 1.2005e-05  -0.36090 0.080194 3.1548e-04
##  Estradiol 0.1nM   0.25761 0.018308 6.3187e-07  -0.10187 0.033549 1.6146e-02
##  Fulvestrant 1nM   2.27742 0.200827 3.2948e-06  -0.64833 0.048452 9.3110e-07
##  Fulvestrant 3nM   4.45609 0.398099 3.6366e-06  -0.44764 0.045055 8.9073e-06
# Facilitation of S cells just about gone with Estradiol, 
# Stronger with Fulvestrant

# The message:
# S cells are hammered by Ribociclib and Fulvestrant.  They do better
# than expected with R cells around because R cells make enough
# Estradiol to overcome the effects of the drugs.

# R cells are harmed by Estradiol itself, and do a little better with S
# cells around.  Not sure why, and why there is no dose response.
# R cells do worse than expected in coculture when treated 
# with Ribociclib or Fulvestrant because of competition with the
# facilitated S cells.

