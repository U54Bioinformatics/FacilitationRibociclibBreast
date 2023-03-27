###########################################################
###################### Read in the data ###################
###########################################################

require(deSolve)
mcf7cols <- c("chartreuse2","red")
names(mcf7cols) <- c("S","R")

source("growthreadm.R")  #Read and fit MCF7 growth data
source("readm.R")     #Read in all other data
source("functions.R")   #Key functions
plotm <- 0  #Control plotting

###########################################################
################### step1, 2, 5: Fit logistic #############
###########################################################

source("step125m.R")  #For MCF7 cells

###########################################################
############## step3: Estimate LV coefficients ###########
###########################################################

source("step3m.R")  #For MCF7 cells

###############################################################
############## step4, 6: Estimate rR and KR  ##################
###############################################################

source("step46m.R")  #For MCF7 cells

