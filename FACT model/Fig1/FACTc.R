require(deSolve)
###########################################################
###################### Read in the data ###################
###########################################################

# Set some nice colors
camacols <- c("chartreuse2","cyan3")
names(camacols) <- c("S","R")

source("allread.R")     #Read in all other data
source("functions.R")   #Key functions
plotm <- 0  #Control plotting

###########################################################
################### step1, 2, 5: Fit logistic #############
###########################################################

source("step125c.R")  #For CAMA-1 cells

###########################################################
############## step3: Estimate LV coefficients ###########
###########################################################

source("step3c.R")  #For CAMA-1 cells

###############################################################
############## step4, 6: Estimate rR and KR  ##################
###############################################################

source("step46c.R")  #For CAMA-1 cells

################################################
############## Organize files ##################
################################################
source("makedat.R")
