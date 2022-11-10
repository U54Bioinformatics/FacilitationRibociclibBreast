rm(list=ls())
require(deSolve) ;require(ggplot2) ; require(tidyr) ; require(dplyr) ; require(data.table) 
require(Matrix); require(fda); library("readr"); require(boot);
require(stringr)

# Summarised results of the fitted full model of estradiol mediated facilitation.
# To be loaded after running the fitting code in :Full Estradiol mediated facilitation model Fit Bayesian Facilitation asymmetric 3D to coculutre spheroid growth trajectories under Ribociclib.R
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/Lab Facilitation/parameterised LH full pars and core data for simulation studies.RData")     

#Load data 
Experiment_Data <- data.table(read.csv(file="~/Dropbox/Vince data/processed data/riboadavo/Sorted_190711_start_ribociclib_exp_analyzed_8-13-19.csv"))
Experiment_Data[,Day:=Day - min(Day)]
Experiment_Data[,Compostition:="polyculture"]
Experiment_Data[Resistance%in%c("100% sensitive","100% resistant"),Compostition:="monoculture"]

# Uniquely label each experimental unit (one microcosm of cells)
experim_unit <- unique(Experiment_Data%>%dplyr::select(DoseNum, Replicate ,    Resistance,Compostition))
experim_unit[,experim_unit:=1:nrow(experim_unit)]
Experiment_Data <- merge(Experiment_Data,experim_unit,by=c("DoseNum", "Replicate" ,    "Resistance", "Compostition"))
Experiment_Data <- Experiment_Data[order(experim_unit)]  #plot(Experiment_Data$X,Experiment(desiredexp)$X)
plot(Experiment_Data$DoseNum,ylab="dose")
plot(as.numeric(Experiment_Data$Resistance),ylab="response")
plot(Experiment_Data$Replicate,ylab="replicate")

# Extract initial conditions and metadata
lu_table<-data.table(Experiment_Data[Day==0]%>%group_by(Resistance, DoseNum)%>%mutate(N=mean(SCellNum),A=mean(RCellNum),"Z_N_"=0,"Z_A_"=0,"S_N_"=0,"S_A_"=0,E=0))
ndose <- length(unique(Experiment_Data$DoseNum))
doses_seq <- unique(Experiment_Data$DoseNum)#rep(0:1,each=1)         

# 1008 obs of: 2 states , 3 compositions , 8 doses 3 reps 7 times    = 2*3*8*3*7
# put initial data in long format ... useful for plotting .. contains metadata for only inits

# ode initial conditions and indexing
inits <- y0 <- c("N_"=lu_table$N,  "A_"=lu_table$A,
                 "Z_N_"=rep(0,length(lu_table$N)) ,
                 "Z_A_"=rep(0,length(lu_table$A)),
                 "S_N_"=rep(0,length(lu_table$N)),
                 "S_A_"=rep(0,length(lu_table$A)),
                 "E_"=lu_table$E )
plot(y0,Info$Initial_Condition);         if(!all(y0==Info$Initial_Condition)){print("DANGER:: Wrong data entry. \n initial conditions not being parsed correctly to stan")}
plot(lu_table$DoseNum,Info[State=="N"]$DoseNum);         if(!all(lu_table$DoseNum==Info[State=="N"]$DoseNum)){print("DANGER:: Wrong data entry. \n dose data not being parsed correctly to stan. \n mismatch identified")}


# Simulation evaluation times (hours) for 21 days
end.day <- max(Experiment_Data$Day)
times   <- seq(0,end.day,by=1)   

### Reformat data for analysis
fit.dd <- data.table(Experiment_Data %>%
                       dplyr::select(X, Day, SCellNum, RCellNum,
                                     Replicate, Compostition, DoseNum, experim_unit) %>%
                       gather(StateOld, value, SCellNum, RCellNum))

# Note that A in the data = sum A+ Z_A
fit.dd[, State:= "A"]
fit.dd[StateOld== "SCellNum", State:= "N"]
fit.dd <- data.table(merge(fit.dd, Info, by= c("DoseNum", "Replicate", "experim_unit", "Compostition", "State")) )
fit.dd <- fit.dd[order(ode_ID, Day)]
fit.dd[, State2:= "TotA"]
fit.dd[State== "N", State2:= "TotN"]

Info[, isNpresent:= TRUE]; Info[, isApresent:= TRUE]
Info[experim_unit%in% which(y0[index_N]==0), isNpresent:= FALSE]
Info[experim_unit%in% which(y0[index_A]==0), isApresent:= FALSE]
Info$State <- gsub("N_", "N",    gsub("A_", "A", Info$State))

pars.outlong <- data.table(gather(pars.out, par, val, r_RbyN:k_RbyN))
pars.outlong$par <- factor(pars.outlong$par, levels = rev(c("gamma_RbyN", "k_RbyN", "r_RbyN", "lambda_RbyN", "B_RbyN", "delta_SRbySN", "delta_RbyN")))

p1Blank<-ggplot(pars.outlong[!par%in%c("delta_SRbySN","delta_RbyN")], aes(y=log(val),x= par ,col=par,fill=par) ) + geom_violin(scale="width",bw=.015)+theme_classic(base_size = 19)+
  theme(aspect.ratio=1,legend.position = "none")+coord_flip()+
  labs(x="Process",y="Resistant cell performance \n (relative to sensitive cells)")+
  geom_hline(yintercept=0,linetype=2,size=2)+
  scale_x_discrete(labels=rev(c("Facilitation \n contribution","Competition \n effect","Division \n (Baseline)","Quiescence \n (Baseline)","Drug sensitivity")))+
  scale_y_continuous(breaks=log(c(0.125,0.25,0.5,1,2,4,8)),labels=c(0.125,0.25,0.5,1,2,4,8))+
  theme(axis.title = element_blank(),axis.text =  element_blank())
ggsave(p1Blank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Resist vs Sensitive BLANK.png")

p1<-ggplot(pars.outlong[!par%in%c("delta_SRbySN","delta_RbyN")], aes(y=log(val),x= par ,col=par,fill=par) ) + geom_violin(scale="width",bw=.015)+theme_classic(base_size = 19)+
  theme(aspect.ratio=1,legend.position = "none")+coord_flip()+
  labs(x="Process",y="Resistant cell performance \n (relative to sensitive cells)")+
  geom_hline(yintercept=0,linetype=2,size=2)+
  scale_x_discrete(labels=rev(c("Facilitation \n contribution","Competition \n effect","Division \n (Baseline)","Quiescence \n (Baseline)","Drug sensitivity")))+
  scale_y_continuous(breaks=log(c(0.125,0.25,0.5,1,2,4,8)),labels=c(0.125,0.25,0.5,1,2,4,8))
ggsave(p1,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Resist vs Sensitive.png")

# Extract max likelihood param set
par.vec_mcmc_max <- as.vector(unlist(pars.out[lp__==max(lp__)]))
names(par.vec_mcmc_max) <-colnames(pars.out)

# Adjust param naming
pars1 <- par.vec_mcmc_max #par.vec_mcmc_i
names(pars1)[names(pars1)=="k_N"] <- "K_N"
names(pars1)[names(pars1)=="k_A"] <- "K_A"
names(pars1)[names(pars1)=="B_N"] <- "K_q_N"
names(pars1)[names(pars1)=="B_A"] <- "K_q_A"
parstt <- pars1
parstt["omega"] <- 1
dY_fun2 <- function(t, y, pars, doses) {
  with( as.list( c(pars, y) ), {
    N <- y[index_N]
    A <- y[index_A]
    Z_N<- y[index_Z_N]
    Z_A<- y[index_Z_A]
    S_N<- y[index_S_N]
    S_A<- y[index_S_A]
    E <- y[index_Eg]
    RC <- (1 - (N + Z_N + S_N )/K_N - (A + Z_A + S_A )/K_A) # Resource limitation
    
    q_x_N <- (doses)/(K_q_N + (doses)) # Quiesce
    q_x_A <- doses/(K_q_A + doses) # Quiesce
    Gphase_N <- r_N*(1 + omega*E/(1+c_N*E))*RC # Enter G1/S checkpoint phase
    Gphase_A <- r_A*(1 + E/(1+c_A*E))*RC # Enter G1/S checkpoint phase
    
    d_N <- ( Gphase_N * (1 - q_x_N) - Gphase_N*q_x_N - lambda_N )*N     #/(1+B_N*doses)
    d_A <- ( Gphase_A * (1 - q_x_A) - Gphase_A*q_x_A - lambda_A )*A     #/(1+B_A*doses)
    
    d_Z_N <- (lambda_N + Gphase_N*q_x_N)*N - delta_N*Z_N
    d_Z_A <- (lambda_A + Gphase_A*q_x_A)*A - delta_A*Z_A 
    
    d_S_N <- delta_N*Z_N - delta_SN*S_N
    d_S_A <- delta_A*Z_A - delta_SA*S_A 
    
    d_E <- gamma_N*N +  gamma_A*A  - delta*E    
    list(c(d_N, d_A,d_Z_N, d_Z_A,d_S_N, d_S_A,d_E))
  } )
}


#rho <- 1
#alpha <- 2
luinter <- expand.grid(rho=1,alpha= 1  ,phi=1, omega=1) 


### Run simulations fixing competition and facilitation
res <- rbindlist(lapply(1:nrow(luinter),function(i){
  cat(i/nrow(luinter)); cat("   ")
  rho<-luinter[i,]$rho
  alpha<-luinter[i,]$alpha
  phi<-luinter[i,]$phi
  omega <- luinter[i,]$omega
  parstt[names(parstt)=="gamma_N"] <- rho*pars1[names(parstt)=="gamma_N"]
  parstt[names(parstt)=="gamma_A"] <- rho*pars1[names(parstt)=="gamma_A"]
  
  parstt[names(parstt)=="K_N"] <- (alpha)*pars1[names(parstt)=="K_N"] 
  parstt[names(parstt)=="K_A"] <- pars1[names(parstt)=="K_A"]
  
  parstt[names(parstt)=="r_N"] <- pars1[names(parstt)=="r_N"] 
  parstt[names(parstt)=="r_A"] <- phi*pars1[names(parstt)=="r_A"]
  parstt["omega"]  <- omega
  CompA_relB <- parstt[names(parstt)=="K_A"]/parstt[names(parstt)=="K_N"]
  # simuate ode add relevant ode state data to merge withExperimental metadata 
  out_N_N <- data.table(ode(y=y0, parms=parstt, times=times, func=dY_fun2,doses=doses))
  out_long0 <- data.table(gather(out_N_N,Variable,y0,-c(time)))
  out_long0[,State:="E"]
  out_long0[Variable%in%c(paste0("N_",1:length(index_N) )),State:="N"];out_long0[Variable%in%c(paste0("Z_N_",1:length(index_N) )),State:="Z_N"];out_long0[Variable%in%c(paste0("S_N_",1:length(index_N) )),State:="S_N"]
  out_long0[Variable%in%c(paste0("A_",1:length(index_N) )),State:="A"];out_long0[Variable%in%c(paste0("Z_A_",1:length(index_A) )),State:="Z_A"];out_long0[Variable%in%c(paste0("S_A_",1:length(index_A) )),State:="S_A"]
  out_long0[, experim_unit := as.numeric(str_extract(Variable, "([0-9]+)"))]
  out_long0$State <- factor(out_long0$State, levels = c("N","A","Z_N", "Z_A", "S_N", "S_A","E"))
  out_long0[,experim_unit:=as.numeric(as.character(experim_unit))]
  out_long0[,Variable:=NULL]
  out_long <- merge(out_long0, Info, by= c("State","experim_unit"))[order(State, experim_unit, time, DoseNum)][Replicate==1]
  out_long$rho <- rho
  out_long$alpha <- alpha
  out_long$phi <- phi
  out_long$omega <- omega
  out_long$CompA_relB <- CompA_relB
  out_longsum <- data.table(gather(out_long %>% 
                                     dplyr::select(-c(Initial_Condition, ode_ID)) %>% 
                                     spread( State,y0)%>%group_by(experim_unit, time) %>%
                                     dplyr::mutate(TotN= N + S_N + Z_N, TotA= A + S_A + Z_A ) , State, y0, A:TotA))
    return(out_longsum)
}))

finabund_res2 <- unique(res[time==max(time)][State!="E"]%>%dplyr::select(-c(experim_unit,Replicate)))
finabund_res <- unique(res[Replicate==1][time==max(time)][State!="E"])
finabund_res[,subclonespresent:="both"]
finabund_res[isNpresent==FALSE,subclonespresent:="A"]
finabund_res[isApresent==FALSE,subclonespresent:="N"]
finabund_res[,is_Totstate:=FALSE]
finabund_res[State%in%c("TotA","TotN"),is_Totstate:=TRUE]
finabund_res <- data.table( finabund_res %>% group_by(Replicate, experim_unit, subclonespresent, rho, alpha, phi, omega, CompA_relB, DoseNum) %>% mutate(Total= sum(y0*TRUE)) )
fit.dd_tt <- fit.dd
fit.dd_tt$State <- fit.dd_tt$State2
qss <- res[Compostition=="polyculture"][rho%in%c(1)][alpha==1][omega==1]%>%spread(State,y0)
qss$qss_E <- (par.vec_mcmc_max["gamma_N"]*(qss$N+qss$S_N)+par.vec_mcmc_max["gamma_A"]*(qss$A+qss$S_A))/par.vec_mcmc_max["delta"]

medE<-  median(qss$E)
max_mu<- 10
nu_N <- max_mu*(medE-par.vec_mcmc_max["gamma_N"] )/(medE*(max_mu-1))
nu_A <- max_mu*(medE-par.vec_mcmc_max["gamma_A"] )/(medE*(max_mu-1))

dY_fun3 <- function(t, y, pars,doses) {
  with( as.list( c(pars, y) ), {
    N <- y[index_N]
    A <- y[index_A]
    Z_N<- y[index_Z_N]
    Z_A<- y[index_Z_A]
    S_N<- y[index_S_N]
    S_A<- y[index_S_A]
    E <-     y[index_Eg]
    EA <-     mu*(gamma_A+nu_A*E) /(mu+nu_A)
    EN <-     mu*(gamma_N+nu_N*E) /(mu+nu_N)
    RC <- (1 - (N + Z_N + S_N )/K_N - (A + Z_A + S_A )/K_A) # Resource limitation
    
    q_x_N <- (doses)/(K_q_N + (doses)) # Quiesce
    q_x_A <- doses/(K_q_A + doses) # Quiesce
    Gphase_N <- r_N*(1 + EN/(1+c_N*EN))*RC # Enter G1/S checkpoint phase
    Gphase_A <- r_A*(1 + EA/(1+c_A*EA))*RC # Enter G1/S checkpoint phase
    
    d_N <- ( Gphase_N * (1 - q_x_N) - Gphase_N*q_x_N - lambda_N )*N     #/(1+B_N*doses)
    d_A <- ( Gphase_A * (1 - q_x_A) - Gphase_A*q_x_A - lambda_A )*A     #/(1+B_A*doses)
    
    d_Z_N <- (lambda_N + Gphase_N*q_x_N)*N - delta_N*Z_N
    d_Z_A <- (lambda_A + Gphase_A*q_x_A)*A - delta_A*Z_A 
    
    d_S_N <- delta_N*Z_N - delta_SN*S_N
    d_S_A <- delta_A*Z_A - delta_SA*S_A 
    
    d_E <- gamma_N*N +  gamma_A*A  - delta*E   #+ Z_N + Z_A
    list(c(d_N, d_A,d_Z_N, d_Z_A,d_S_N, d_S_A,d_E))
  } )
}

luinter<-expand.grid(rho=1,alpha= 1,mu=c(1:max_mu),nu_A=nu_A,nu_N=nu_N
                     ,phi=1)

### Run simulations varying competition and facilitation
res <- rbindlist(lapply(1:nrow(luinter),function(i){
  cat(i/nrow(luinter)); cat("   ")
  rho<-luinter[i,]$rho
  alpha<-luinter[i,]$alpha
  phi<-luinter[i,]$phi
  nu_A <- luinter[i,]$nu_A
  nu_N <- luinter[i,]$nu_N
  mu <- luinter[i,]$mu
  parstt[names(parstt)=="gamma_N"] <- rho*pars1[names(parstt)=="gamma_N"]
  parstt[names(parstt)=="gamma_A"] <- rho*pars1[names(parstt)=="gamma_A"]
  
  parstt[names(parstt)=="K_N"] <- (alpha)*pars1[names(parstt)=="K_N"] 
  parstt[names(parstt)=="K_A"] <- pars1[names(parstt)=="K_A"]
  
  parstt[names(parstt)=="r_N"] <- pars1[names(parstt)=="r_N"] 
  parstt[names(parstt)=="r_A"] <- phi*pars1[names(parstt)=="r_A"]
  parstt["nu_A"]  <- nu_A
  parstt["nu_N"]  <- nu_N
  parstt["mu"]  <- mu
  CompA_relB <- parstt[names(parstt)=="K_A"]/parstt[names(parstt)=="K_N"]

  # simuate ode add relevant ode state data to merge withExperimental metadata 
  out_N_N <- data.table(ode(y=y0, parms=parstt, times=times, func=dY_fun3,doses=doses))
  out_long0 <- data.table(gather(out_N_N,Variable,y0,-c(time)))
  out_long0[,State:="E"]
  out_long0[Variable%in%c(paste0("N_",1:length(index_N) )),State:="N"];out_long0[Variable%in%c(paste0("Z_N_",1:length(index_N) )),State:="Z_N"];out_long0[Variable%in%c(paste0("S_N_",1:length(index_N) )),State:="S_N"]
  out_long0[Variable%in%c(paste0("A_",1:length(index_N) )),State:="A"];out_long0[Variable%in%c(paste0("Z_A_",1:length(index_A) )),State:="Z_A"];out_long0[Variable%in%c(paste0("S_A_",1:length(index_A) )),State:="S_A"]
  out_long0[, experim_unit := as.numeric(str_extract(Variable, "([0-9]+)"))]
  out_long0$State <- factor(out_long0$State, levels = c("N","A","Z_N", "Z_A", "S_N", "S_A","E"))
  out_long0[,experim_unit:=as.numeric(as.character(experim_unit))]
  out_long0[,Variable:=NULL]
  out_long <- merge(out_long0,Info,by= c("State","experim_unit"))[order(State,experim_unit,time,DoseNum)][Replicate==1]
  out_long$rho <- rho
  out_long$alpha <- alpha
  out_long$phi <- phi
  out_long$mu <- mu
  out_long$nu_A <- nu_A
  out_long$nu_N <- nu_N
  out_long$CompA_relB <- CompA_relB
  out_longsum <- data.table(gather(out_long %>% 
                                     dplyr::select(-c(Initial_Condition, ode_ID)) %>% 
                                     spread( State,y0)%>%group_by(experim_unit, time) %>%
                                     dplyr::mutate(TotN= N + S_N + Z_N, TotA= A + S_A + Z_A ) , State,y0,A:TotA))
  return(out_longsum)
}))

finabund_res2 <- unique(res[time==max(time)][State!="E"]%>%dplyr::select(-c(experim_unit,Replicate)))
finabund_res <- unique(res[Replicate==1][time==max(time)][State!="E"])
finabund_res[,subclonespresent:="both"]
finabund_res[isNpresent==FALSE,subclonespresent:="A"]
finabund_res[isApresent==FALSE,subclonespresent:="N"]
finabund_res[,is_Totstate:=FALSE]
finabund_res[State%in%c("TotA","TotN"),is_Totstate:=TRUE]
finabund_res<-data.table( finabund_res%>%group_by(Replicate,experim_unit,subclonespresent,rho,alpha,phi,mu,nu_A,nu_N,CompA_relB,DoseNum)%>%mutate(Total=sum(y0*TRUE)) )
fit.dd_tt<-fit.dd
fit.dd_tt$State<-fit.dd_tt$State2
qss<-res[Compostition=="polyculture"][rho%in%c(1)][alpha==1][mu==100]%>%spread(State,y0)
qss$qss_E<-(par.vec_mcmc_max["gamma_N"]*(qss$N+qss$S_N)+par.vec_mcmc_max["gamma_A"]*(qss$A+qss$S_A))/par.vec_mcmc_max["delta"]

cell.labs <- c("Resistant", "Sensitive")
names(cell.labs) <- c("TotA","TotN")

ggplot( res[Compostition=="polyculture"][rho%in%c(0.0,0.5,1)][alpha==1][time==max(time)][State%in%c("TotA","TotN")][mu%in%c(1:5)],aes(y=log(1+y0) ,x=DoseNum))+
  geom_smooth(data=res[Compostition=="polyculture"][rho%in%c(1)][alpha==1][time==max(time)][State%in%c("TotA","TotN")][mu%in%c(1:5)][mu==max(mu)],
              se=F,method="gam",formula=y~s(x,k=5),col="black",aes(size=200.5,group=interaction(mu,State)))+ 
  geom_smooth(se=F,method="gam",formula=y~s(x,k=5),aes(col=State , fill=State,size=(mu)/5,group=interaction(mu,State)))+ scale_size(name="Facilitation (%)",  #breaks=c(0,0.5,1),labels = 100*(c(0,0.5,1)), range = c(0.5, 1.4))+
                                                                                                            breaks=c(0,0.2,0.4,0.6,0.8,1),labels = 100*(c(0,0.2,0.4,0.6,0.8,1)), range = c(0.5, 2.04))+
  facet_grid(~State,labeller = labeller(State = cell.labs))+theme_classic(base_size=20)+
  geom_point(data=fit.dd_tt[Day==21][Compostition=="polyculture"], aes(y=log(1+value) ,x=DoseNum ,fill=State2 ),colour="black",pch=21,size=2.5)+
  scale_color_manual(values=c("red","green","blue"),guide = FALSE)+ 
  scale_fill_manual(values=c("red","green","blue"),guide = FALSE)+
  geom_hline(yintercept = log(1+max(y0) ),col="slategrey", linetype="dashed",size=1.5)+theme(aspect.ratio=1)+
  labs(y="Final abundance",x="Drug dose")+scale_y_continuous(labels=
                                                               function(x) format(exp(x)-1, scientific = FALSE)
                                                             ,breaks=log(1+c(0,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8)))

Fig5DModelPredictionData <- res[Compostition=="polyculture"][rho%in%c(0.0,0.5,1)][alpha==1][time==max(time)][State%in%c("TotA","TotN")][mu%in%c(1:5)]
Fig5DObservedData <-fit.dd_tt[Day==21][Compostition=="polyculture"]
#save( Fig5DModelPredictionData, Fig5DObservedData  ,file="/Users/jason/Dropbox/Cancer_pheno_evo/data/Lab Facilitation/Fig5Ddata.RData")

p2<-ggplot( res[Compostition=="polyculture"][rho%in%c(0.0,0.5,1)][alpha==1][time==max(time)][State%in%c("TotA","TotN")][mu%in%c(1:5)],aes(y=log(1+y0) ,x=DoseNum))+
  geom_smooth(data=res[Compostition=="polyculture"][rho%in%c(1)][alpha==1][time==max(time)][State%in%c("TotA","TotN")][mu%in%c(1:5)][mu==max(mu)],
              se=F,method="gam",formula=y~s(x,k=5),col="black",aes(group=interaction(mu,State)),size=2)+ 
  geom_smooth(se=F,method="gam",formula=y~s(x,k=5),aes(col=(mu)/5 , fill=(mu)/5,group=interaction(mu,State)),size=1)+ 
  facet_grid(~State,labeller = labeller(State = cell.labs))+theme_classic(base_size=20)+
  geom_point(data=fit.dd_tt[Day==21][Compostition=="polyculture"], aes(y=log(1+value) ,x=DoseNum ),colour="black",size=2.5)+
  scale_color_viridis_c(name="Facilitation (%)" ,option = "B",end = 0.9, direction=1,breaks=c(0,0.2,0.4,0.6,0.8,1),labels = 100*(c(0,0.2,0.4,0.6,0.8,1)) )+ 
  scale_fill_viridis_c(name="Facilitation (%)" ,option = "B",end = 0.9,direction=1, breaks=c(0,0.2,0.4,0.6,0.8,1),labels = 100*(c(0,0.2,0.4,0.6,0.8,1))  )+
  geom_hline(yintercept = log(1+max(y0) ),col="slategrey", linetype="dashed",size=1.5)+theme(aspect.ratio=1)+
  labs(y="Final abundance",x="Drug dose")+scale_y_continuous(labels=
                                                               function(x) format(exp(x)-1, scientific = FALSE)
                                                             ,breaks=log(1+c(0,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8)))
ggsave(p2,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Fulv block.png")

pleg <- ggpubr::get_legend(p2+theme(legend.text = element_blank(),legend.title = element_blank()))
p2LegBlank <- ggpubr::as_ggplot(pleg)
ggsave(p2LegBlank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Fulv block Legend.png",height=2,width = 2)

p2Blank <- ggplot( res[Compostition=="polyculture"][rho%in%c(0.0,0.5,1)][alpha==1][time==max(time)][State%in%c("TotA","TotN")][mu%in%c(1:5)],aes(y=log(1+y0) ,x=DoseNum))+
  geom_smooth(data=res[Compostition=="polyculture"][rho%in%c(1)][alpha==1][time==max(time)][State%in%c("TotA","TotN")][mu%in%c(1:5)][mu==max(mu)],
              se=F,method="gam",formula=y~s(x,k=5),col="black",aes(group=interaction(mu,State)),size=2)+ 
  geom_smooth(se=F,method="gam",formula=y~s(x,k=5),aes(col=(mu)/5 , fill=(mu)/5,group=interaction(mu,State)),size=1)+ 
  facet_grid(~State,labeller = labeller(State = cell.labs))+theme_classic(base_size=20)+
  geom_point(data=fit.dd_tt[Day==21][Compostition=="polyculture"], aes(y=log(1+value) ,x=DoseNum ),colour="black",size=2.5)+
  scale_color_viridis_c(name="Facilitation (%)" ,option = "B",end = 0.9, direction=1,breaks=c(0,0.2,0.4,0.6,0.8,1),labels = 100*(c(0,0.2,0.4,0.6,0.8,1)) )+ 
  scale_fill_viridis_c(name="Facilitation (%)" ,option = "B",end = 0.9,direction=1, breaks=c(0,0.2,0.4,0.6,0.8,1),labels = 100*(c(0,0.2,0.4,0.6,0.8,1))  )+
  geom_hline(yintercept = log(1+max(y0) ),col="slategrey", linetype="dashed",size=1.5)+theme(aspect.ratio=1)+
  labs(y="Final abundance",x="Drug dose")+scale_y_continuous(labels=
                                                               function(x) format(exp(x)-1, scientific = FALSE)
                                                             ,breaks=log(1+c(0,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8)))+
  theme(axis.title = element_blank(),axis.text =  element_blank(),strip.background = element_blank(),
        strip.text.x = element_blank(),legend.position="none")

ggsave(p2Blank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Fulv block BLANK.png")

