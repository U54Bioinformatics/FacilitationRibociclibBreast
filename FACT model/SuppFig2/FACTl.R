###########################################################
###################### Read in the data ###################
###########################################################

require(deSolve)
LY2cols <- c("chartreuse2","darkblue")
names(LY2cols) <- c("S","R")

source("growthreadl.R")  #Read and fit LY2 growth data
source("readl.R")     #Read in all other data
source("functions.R")   #Key functions
plotm <- 0  #Control plotting

###########################################################
################### step1, 2, 5: Fit logistic #############
###########################################################

source("step125l.R")  #For LY2 cells

###########################################################
############## step3: Estimate LV coefficients ###########
###########################################################

source("step3l.R")  #For LY2 cells

###############################################################
############## step4, 6: Estimate rR and KR  ##################
###############################################################

source("step46l.R")  #For LY2 cells

