rm(list=ls())
require(deSolve) ;require(ggplot2) ; require(tidyr) ; require(dplyr) ; require(data.table) 
require(Matrix); require(fda); library("readr"); require(boot);
require(stringr)

#Load data 
everolimus_dd <- data.table(read.csv(file = "~/Dropbox/Vince data/processed data/Sorted_everolimus_coculture_4-10-19_5-2-19.csv"))
new_ribo_dd <- data.table(read.csv(file="~/Dropbox/Vince data/processed data/riboadavo/Sorted_190711_start_ribociclib_exp_analyzed_8-13-19.csv"))
ribo_dd <- data.table(read.csv(file="~/Dropbox/Vince data/processed data/riboadavo/Sorted_Ribo_Adavo.csv"))
#calc_params <- data.table(read.csv(file = "~/Dropbox/Vince data/Fitted Model Params/Everolimus Params/Final6Models_params.csv"))
#top5 <- data.table(read.csv(file="~/Dropbox/Vince data/Fitted Model Params/Top5PureCultureParams.csv"))
#happycells_dd <- data.table(read.csv(file="~/Dropbox/Vince data/Fitted Model Params/GrowthFacilitation2.csv"))
# #Function to determine what model will be used in the loop
# Model <- function(model_nm){
#   switch(model_nm, 
#          "CE"=CE_benefits,
#          "CE_GC"=CE_GC_benefits,
#          "CE_alpha"=CE_alpha_benefits,
#          "Growth_Facilitation"=Growth_Facilitation_benefits)
# }

#Function to determine what experimental data will be used
Experiment <- function(type){
  switch(type,
         "ribo" = ribo_dd,
         "everolimus" = everolimus_dd,
         "new_ribo" = new_ribo_dd)
}

#Select Data to use... it is in a semi wide format (not long nor a dense matrix)
desiredexp <- "new_ribo"#"everolimus" #options: "ribo", "everolimus", "new_ribo"
Experiment_Data <- Experiment(desiredexp) ; if(desiredexp == "ribo"){ Experiment_Data <- Experiment_Data[Drug == "ribociclib"] }
Experiment_Data <- Experiment_Data#[Day!=0]
Experiment_Data[,Day:=Day - min(Day)]
Experiment_Data[,Compostition:="polyculture"]
Experiment_Data[Resistance%in%c("100% sensitive","100% resistant"),Compostition:="monoculture"]
#Experiment("ribo")[Drug == "ribociclib"]

# uniquely label each experimental unit (one microcosm of cells)
experim_unit <- unique(Experiment_Data%>%dplyr::select(DoseNum, Replicate ,    Resistance,Compostition))
experim_unit[,experim_unit:=1:nrow(experim_unit)]
Experiment_Data <- merge(Experiment_Data,experim_unit,by=c("DoseNum", "Replicate" ,    "Resistance", "Compostition"))
Experiment_Data <- Experiment_Data[order(experim_unit)]  #plot(Experiment_Data$X,Experiment(desiredexp)$X)
plot(Experiment_Data$DoseNum,ylab="dose")
plot(as.numeric(Experiment_Data$Resistance),ylab="response")
plot(Experiment_Data$Replicate,ylab="replicate")

# extract initial conditions and metadata
lu_table<-data.table(Experiment_Data[Day==0]%>%group_by(Resistance, DoseNum)%>%mutate(N=mean(SCellNum),A=mean(RCellNum),"Z_N_"=0,"Z_A_"=0,"S_N_"=0,"S_A_"=0,E=0))
#setnames(lu_table,old=c("DoseNum","Replicate"),new=c("dose","replicate"))

ndose <- length(unique(Experiment_Data$DoseNum))
doses_seq <- unique(Experiment_Data$DoseNum)#rep(0:1,each=1)          #lu_table[, ID:=(1:nrow(lu_table))]    #lu_table <- lu_table[order(ID)]

# 1008 obs of: 2 states , 3 compositions , 8 doses 3 reps 7 times    = 2*3*8*3*7
# put initial data in long format ... useful for plotting .. contains metadata for only inits
Info <- data.table(gather(lu_table%>% dplyr::select( N,A, E,Z_N_,Z_A_,S_N_,S_A_, DoseNum, Replicate,experim_unit,Compostition),State,Initial_Condition, -c(DoseNum,Replicate,experim_unit,Compostition)) )
Info$State <- factor(Info$State, levels = c("N","A","Z_N_","Z_A_","S_N_","S_A_","E"))
Info <- Info[order(State,experim_unit)]
Info[ ,ode_ID:=1:nrow(Info)]

# ode initial conditions and indexing
inits <- y0 <- c("N_"=lu_table$N,  "A_"=lu_table$A,
                 "Z_N_"=rep(0,length(lu_table$N)) ,
                 "Z_A_"=rep(0,length(lu_table$A)),
                 "S_N_"=rep(0,length(lu_table$N)),
                 "S_A_"=rep(0,length(lu_table$A)),
                 "E_"=lu_table$E )
plot(y0,Info$Initial_Condition);         if(!all(y0==Info$Initial_Condition)){print("DANGER:: Wrong data entry. \n initial conditions not being parsed correctly to stan")}
index_N <- grep("N_",names(y0))[!grep("N_",names(y0))%in%c(grep("Z_N_",names(y0)),grep("S_N_",names(y0)) )]
index_A <- grep("A_",names(y0))[!grep("A_",names(y0))%in%c(grep("Z_A_",names(y0)),grep("S_A_",names(y0)) )]
index_Z_N_ <- grep("Z_N_",names(y0))
index_Z_A_ <- grep("Z_A_",names(y0))
index_S_N_ <- grep("S_N_",names(y0))
index_S_A_ <- grep("S_A_",names(y0))

index_E <- grep("E_",names(y0))          ; if( length(unique(Experiment_Data$experim_unit))!= length(inits)/7){ print("DANGER::incorrect number of initial conditions specified")}

# ode dose vector (length==number of experim_unit)
doses<- lu_table$DoseNum         ;        if(length(doses)!=nrow(experim_unit)){print("DANGER:: Wrong data entry. \n number of ode doses doesnt match number of experim units")}

plot(lu_table$DoseNum,Info[State=="N"]$DoseNum);         if(!all(lu_table$DoseNum==Info[State=="N"]$DoseNum)){print("DANGER:: Wrong data entry. \n dose data not being parsed correctly to stan. \n mismatch identified")}


# Simulation evaluation times (hours) for 21 days
end.day <- max(Experiment_Data$Day)
times   <- seq(0,end.day,by=1)   

### Reformat data for analysis
fit.dd <- data.table(Experiment_Data%>%
                       dplyr::select(X,Day, SCellNum,RCellNum,
                                     Replicate,Compostition,DoseNum,experim_unit)%>%
                       gather(StateOld,value,SCellNum,RCellNum))
#fit.dd$StateOld

# Note that A in the data = sum A+ Z_A
fit.dd[,State:="A"]
fit.dd[StateOld=="SCellNum",State:="N"]
fit.dd <- data.table(merge(fit.dd,Info,by=c("DoseNum","Replicate","experim_unit","Compostition","State")) )
fit.dd <- fit.dd[order(ode_ID,Day)]
stan_Y <-  fit.dd[Day>0]%>%dplyr::select(ode_ID,Day,value)%>%spread(ode_ID,value)%>%select(-Day)
#image( t(as.matrix(log(1+stan_Y[,1:50]))))
if( ncol(stan_Y)!=2*nrow(lu_table)){print("DANGER:: Wrong data entry. \n number of observed ode states doesnt match the metadata")}
if( nrow(stan_Y)!=length(unique(Experiment_Data$Day))-1){print("DANGER:: Wrong data entry. \n number of time points given to ode  doesnt match the number observed after initial conditions")}

ggplot(fit.dd[value>0],aes(y=log(value),x=Day,group=X,col=Compostition))+
  geom_point()+
  facet_grid(State~DoseNum)+theme_classic()


# Stan requires a list of data input objects
obs.times <- unique(fit.dd[Day>0]$Day) # in days
x_iSetts  <- 1:length(index_N)   # index for each populations      #x_rSetts  <- c(Info$dose)        # drug doses of each population
pred_t <- min(obs.times):max(obs.times) # Range of times to generate predictions
pred_t <- 1:max(obs.times) # Range of times to generate predictions

stan_data <- list(nt=length(obs.times),          # Number of days sampled
                  nStates=  length(y0),          # Number of ode states=seen or otherwise
                  nStatesObs=2*length(x_iSetts), # Number of observed states
                  t0=0 ,                         # Start time
                  tObs=obs.times ,               # Observation times       
                  inits=inits,                   # Initial conditions
                  x_r= doses  ,  nSettings=length(doses),
                  x_i=x_iSetts   ,  nSettingsi=length(x_iSetts),
                  Y=  stan_Y,
                  is_non_zero_state=unname(colSums(stan_Y)>0)*1 ,
                  nPred=length(pred_t),            # Prediction variables
                  pred_ts=pred_t )

load(file="Ribociclib facilitation life history FULL NQS.Rdata" )



require(rstan)
extrct_preds<- extract(test,pars="pred_I")$ pred_I  #str(extrct_preds) #[1:mcmc_iter, 1:npredtimes, 1:ntotalstates(unobserved)] 
post_pred <- data.table(
  rbindlist(lapply(1:length(inits),function(state_i){
    s_mcmc_t <-extrct_preds[,,state_i] 
    # set names of cols to be time of prediction
    #colnames(s_mcmc_t)<- pred_t
    x <- data.table(State_numb= names(y0)[state_i] ,s_mcmc_t)
    x[,ode_ID:=state_i]
    x[,mcmc_id:=1:nrow(x)]
    # gather into long format stack times which are columns of data
    
  }))
)

pars.out<-data.table(bind_rows(extract(test,pars=c("r_N","k_N","B_N" ,"c_N" ,  "gamma_N" ,"delta_N","lambda_N","delta_SN",
                                                   "r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta_A","lambda_A","delta_SA","delta","sigma","lp__"))))[lp__ >0]

preds <- data.table(gather(cbind(post_pred,lp=extract(test,pars="lp__")[[1]] ),key,value,-mcmc_id,-State_numb,-lp,-ode_ID))
preds[,State:=(gsub("[[:digit:]]", "", State_numb, perl = TRUE))]
preds[,experim_unit:= gsub(State, "", State_numb),by=c("ode_ID" ,"mcmc_id","State_numb" )]
preds[,State:=sub("_$","",State)]

preds[, timeid:= as.numeric(as.character(str_remove(key, "V")))]
preds2 <- data.table(merge(preds, data.table(Day=pred_t,timeid=1:length(pred_t)),by="timeid"))
preds2[,experim_unit:=as.numeric(as.character(experim_unit))]

preds2_wide <- data.table( preds2 %>%
                             select(-State_numb,-lp,-key,-ode_ID)%>%
                             
                             spread(State,value) ) 
preds2_wide[,Ntot:=N+Z_N+S_N]
preds2_wide[,Atot:=A+Z_A+S_A]
preds2_wide <- preds2_wide%>%select(-c(A,E,N,Z_A,Z_N,S_A,S_N))

setnames(preds2_wide,old=c("Ntot", "Atot"),new=c("N", "A"))
luinfo<-unique(preds2 %>% select(State_numb,lp,key,ode_ID,experim_unit, Day,mcmc_id))

rm(list="pars.out");rm(list="preds")
rm(list="test")
rm(list="preds2")
#rm(list="preds2_wide")

preds3_wide <- merge(preds2_wide, luinfo, by=c("experim_unit", "Day","mcmc_id"))
preds3 <- data.table(gather(preds3_wide, State,value,N : A) )

abc <- merge(preds3,
             unique(fit.dd%>%
                      dplyr::select(experim_unit,DoseNum,Replicate,Compostition,State,Initial_Condition,ode_ID)
             ),
             by=c("ode_ID","State","experim_unit"))


load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/Lab Facilitation/parameterised LH full pars and core data for simulation studies.RData")     

abc[ ,noised:= exp(rnorm(nrow(abc),log(value),sd=mean(pars.out$sigma)))]


mean_credibleregion<-data.table(abc%>%group_by(ode_ID,State,experim_unit, Day,timeid ,State_numb,
                                               DoseNum ,Replicate, Compostition, Initial_Condition)%>%
                                  summarise(lcl=quantile(value, probs = 0.025),
                                            ucl=quantile(value, probs = 0.975),
                                            lcl_sigma=quantile(noised, probs = 0.05,na.rm=T),
                                            ucl_sigma=quantile(noised, probs = 0.95,na.rm=T),
                                            mu=mean(value)))



fit.dd[,StateLab:="Adapted"]
fit.dd[State=="N",StateLab:="Naive"]
mean_credibleregion[,StateLab:="Adapted"]
mean_credibleregion[State=="N",StateLab:="Naive"]

fit.dd[,DoseLab:=paste0(DoseNum," nM")]

fit.dd$DoseLab <- factor(fit.dd$DoseLab, levels = paste0(unique(fit.dd$DoseNum)," nM"))

mean_credibleregion[,DoseLab:=paste0(DoseNum," nM")]
mean_credibleregion$DoseLab <- factor(mean_credibleregion$DoseLab, levels = paste0(unique(fit.dd$DoseNum)," nM"))


fit.dd[,StateLab2:="Sensitive"]
fit.dd[StateLab=="Adapted",StateLab2:="Resistant"]
mean_credibleregion[,StateLab2:="Sensitive"]
mean_credibleregion[StateLab=="Adapted",StateLab2:="Resistant"]
mean_credibleregion$StateLab2 <- factor(mean_credibleregion$StateLab2, levels = c("Resistant","Sensitive"))
fit.dd$StateLab2 <- factor(fit.dd$StateLab2, levels = c("Resistant","Sensitive"))


fit.dd[,Compostition2:="Coculture"]
fit.dd[Compostition=="monoculture",Compostition2:="Monoculture"]
mean_credibleregion[,Compostition2:="Coculture"]
mean_credibleregion[Compostition=="monoculture",Compostition2:="Monoculture"]
mean_credibleregion$Compostition2 <- factor(mean_credibleregion$Compostition2, levels = c("Monoculture","Coculture"))
fit.dd$Compostition2 <- factor(fit.dd$Compostition2, levels = c("Monoculture","Coculture"))



p5 <- ggplot(fit.dd[value>0], aes(y=log(value),x=Day,group=DoseLab))+
  #geom_point(aes(group=X,col=Compostition))+
  facet_grid(StateLab2~Compostition2)+
  theme_classic(base_size=21)+
  #geom_point(data=out_long[y0>0&State!="E"],aes(y=log(1+y0) , x=time))+
  
  geom_ribbon(data=mean_credibleregion[mu>10]
              ,aes(y=log(mu),ymax=log(ucl_sigma),ymin=log(lcl_sigma),x=Day,fill=DoseNum,group=interaction(ode_ID,DoseLab)),colour =NA,alpha=0.1)+
  # geom_ribbon(data=mean_credibleregion[mu>10]
  #             ,aes(y=log(mu),ymax=log(ucl),ymin=log(lcl),x=Day,fill=Compostition,group=interaction(ode_ID)),colour =NA,alpha=0.3)+
  labs(y="Cell count (ln)",x="Day")+
  geom_point(data=fit.dd[value>0],aes(group=X,col=DoseNum))+theme(aspect.ratio = 1)+
  geom_line(data=mean_credibleregion[mu>10],
            aes(y=log(mu),x=Day,col=DoseNum,group=interaction(ode_ID)))+
  
  scale_fill_gradientn(name= "Dose", colours=rainbow(100,end =0.8))+
  scale_color_gradientn(name= "Dose",colours=rainbow(100,end =0.8))

ggsave(p5,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Predictions with uncertainty.png")



p5leg <- ggpubr::get_legend(p5+theme(legend.text = element_blank(),legend.title = element_blank()))
p5LegBlank <- ggpubr::as_ggplot(p5leg)
ggsave(p5LegBlank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Predictions with uncertainty Legend.png",height=2,width = 2)

p5Blank <- ggplot(fit.dd[value>0], aes(y=log(value),x=Day,group=DoseLab))+
  #geom_point(aes(group=X,col=Compostition))+
  facet_grid(StateLab2~Compostition2)+
  theme_classic(base_size=21)+
  #geom_point(data=out_long[y0>0&State!="E"],aes(y=log(1+y0) , x=time))+
  
  geom_ribbon(data=mean_credibleregion[mu>10]
              ,aes(y=log(mu),ymax=log(ucl_sigma),ymin=log(lcl_sigma),x=Day,fill=DoseNum,group=interaction(ode_ID,DoseLab)),colour =NA,alpha=0.1)+
  # geom_ribbon(data=mean_credibleregion[mu>10]
  #             ,aes(y=log(mu),ymax=log(ucl),ymin=log(lcl),x=Day,fill=Compostition,group=interaction(ode_ID)),colour =NA,alpha=0.3)+
  labs(y="Cell count (ln)",x="Day")+
  geom_point(data=fit.dd[value>0],aes(group=X,col=DoseNum))+theme(aspect.ratio = 1)+
  geom_line(data=mean_credibleregion[mu>10],
            aes(y=log(mu),x=Day,col=DoseNum,group=interaction(ode_ID)))+
  
  scale_fill_gradientn(name= "Dose", colours=rainbow(100,end =0.8))+
  scale_color_gradientn(name= "Dose",colours=rainbow(100,end =0.8))+theme(axis.title = element_blank(),axis.text =  element_blank(),strip.background = element_blank(),
        strip.text = element_blank(),legend.position="none")

ggsave(p5Blank,filename ="/Users/jason/Dropbox/Cancer_pheno_evo/images and presentations/Lab Facilitation/Predictions with uncertainty BLANK.png")


