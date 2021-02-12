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
lu_table<-data.table(Experiment_Data[Day==0]%>%group_by(Resistance, DoseNum)%>%mutate(N=mean(SCellNum),A=mean(RCellNum)))
#setnames(lu_table,old=c("DoseNum","Replicate"),new=c("dose","replicate"))

ndose <- length(unique(Experiment_Data$DoseNum))
doses_seq <- unique(Experiment_Data$DoseNum)#rep(0:1,each=1)          #lu_table[, ID:=(1:nrow(lu_table))]    #lu_table <- lu_table[order(ID)]

# 1008 obs of: 2 states , 3 compositions , 8 doses 3 reps 7 times    = 2*3*8*3*7
# put initial data in long format ... useful for plotting .. contains metadata for only inits
Info <- data.table(gather(lu_table%>% dplyr::select( N,A, DoseNum, Replicate,experim_unit,Compostition),State,Initial_Condition, -c(DoseNum,Replicate,experim_unit,Compostition)) )
Info$State <- factor(Info$State, levels = c("N","A"))
Info <- Info[order(State,experim_unit)]
Info[ ,ode_ID:=1:nrow(Info)]

# ode initial conditions and indexing
inits <- y0 <- c("N_"=lu_table$N,  "A_"=lu_table$A )
plot(y0,Info$Initial_Condition);         if(!all(y0==Info$Initial_Condition)){print("DANGER:: Wrong data entry. \n initial conditions not being parsed correctly to stan")}
index_N <- grep("N_",names(y0))
index_A <- grep("A_",names(y0))
if( length(unique(Experiment_Data$experim_unit))!= length(inits)/2){ print("DANGER::incorrect number of initial conditions specified")}

# ode dose vector (length==number of experim_unit)
doses<- lu_table$DoseNum         ;        if(length(doses)!=nrow(experim_unit)){print("DANGER:: Wrong data entry. \n number of ode doses doesnt match number of experim units")}

plot(lu_table$DoseNum,Info[State=="N"]$DoseNum);         if(!all(lu_table$DoseNum==Info[State=="N"]$DoseNum)){print("DANGER:: Wrong data entry. \n dose data not being parsed correctly to stan. \n mismatch identified")}


# Simulation evaluation times (hours) for 21 days
end.day <- max(Experiment_Data$Day)
times <- seq(0,end.day,by=1)   

### Reformat data for analysis
fit.dd <- data.table(Experiment_Data%>%
                       dplyr::select(X,Day, SCellNum,RCellNum,
                                     Replicate,Compostition,DoseNum,experim_unit)%>%
                       gather(StateOld,value,SCellNum,RCellNum))


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
  facet_grid(State~DoseNum)


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



Cancer_mod <- "functions {

real[] dY_funStan(real t, real[] y, real[] params, real[] x_r, int[] x_i) {  
real N; real A;  real doses;

real r_N; real B_N;  real k_N;  real delta_N;  
real r_A; real B_A;  real k_A;  real delta_A;   
real dydt[ 2*size(x_i) ]; int s_i; 

s_i = size(x_i);
r_N = params[1]; B_N = params[2]; k_N = params[3];  delta_N = params[4]; 
r_A = params[5]; B_A = params[6]; k_A = params[7];  delta_A = params[8];

// Rate of change functions: dy/dt
for(i in 1:s_i){ 
N = y[i];   A = y[( s_i + i )]; doses = x_r[i];
dydt[i]         = (r_N*N/( 1 + B_N*doses ))*(1 - N/k_N - A/k_A) - delta_N*N ;
dydt[(s_i + i)] = (r_A*A/( 1 + B_A*doses ))*(1 - N/k_N - A/k_A) - delta_A*A ;

}
return dydt;
}

}
data {
int<lower = 1> nt;          // Number of days sampled
int<lower = 1> nStates;     // number of ODE states
int<lower = 1> nStatesObs;  // number of OBSERVED ODE states
real t0;                    // Initial time point (zero)
real<lower=0> tObs[nt];     // Times of observations
real<lower = 0> inits[nStates]; 
real<lower=0> Y [nt,nStatesObs] ;   // event response variable :Approx Normally distributed data

int<lower = 1> nSettings;   // Number of experimental timing settings
real x_r[nSettings];        // Experimental timing settings
int<lower = 1> nSettingsi;  // Number of integer element indexes for ode solver
int x_i[nSettingsi];        // Integer vector for ode solver
int is_non_zero_state[nStatesObs]; // Binary variable specifying if state is always zero

int<lower = 1> nPred;               // This is to generate predicted/unsampled data
real pred_ts[nPred];                // Time points for predicted/unsampled data
}
transformed data {
real lnY [nt,nStatesObs] ;   // <lower=0>  response variable :Approx Normally distributed data
for(i in 1:nt){for(j in 1:nStatesObs){ lnY[i,j] = log(Y[i,j] + 1.0 ); }}  //for(j in 1:nStatesObs ){   if ( is_non_zero_state[j]==1 ){ lnY[,j] = log(Y[,j])  ;    }} 
//delta_raw= log(1e-25);                                                //delta_N_raw=log(2.0e-01);//delta_A_raw=log(2.0e-01);//c_N_raw= log(1);//c_A_raw= log(1.2);//gamma_N_raw= log(2e-5);//gamma_A_raw= log(4.5e-5);//B_N_raw=log(2.5) ;//B_A_raw=log(5e-2) ; //k_N_raw=log(1.0e+6-1)-10 ; //k_A_raw=log(1.2e+6-1)-10 ;
}
parameters {
real r_N_raw; real k_N_raw; real B_N_raw;  real delta_N_raw; 
real r_A_raw; real k_A_raw; real B_A_raw;  real delta_A_raw; 
real sigma_raw;

}
transformed parameters{
real<lower=0> r_N; real<lower=0> k_N; real<lower=0> B_N; real<lower=0> delta_N;
real<lower=0> r_A; real<lower=0> k_A; real<lower=0> B_A; real<lower=0> delta_A;
real<lower=0> sigma; real params [8];          // parameters need to be passed to ODE integrator as a real vector

real y_hat[nt, nStates];         // Output from the ODE solver

r_N= exp(r_N_raw);
k_N= 1 + exp(10 + k_N_raw);
B_N= exp(5+B_N_raw);
delta_N= exp(delta_N_raw+1); 
r_A= exp(r_A_raw);
k_A= 1 + exp(10 + k_A_raw);
B_A= exp(5+B_A_raw);
delta_A= exp(delta_A_raw+1); 
sigma= exp(sigma_raw);

params[1] = r_N ; params[2] = B_N ; params[3] = k_N ; params[4] = delta_N ; 
params[5] = r_A ; params[6] = B_A ; params[7] = k_A ; params[8] = delta_A ; 
//integrate ode
y_hat = integrate_ode_bdf(dY_funStan, inits, t0, tObs, params, x_r, x_i , 1e-2, 1e-1,1e2);//1e-6, 1e-1, 1e3 ) ;        //integrate_ode_rk45  // pos use:: 1e-2, 1e-2, 5e3// 1e-4,1e-1,5e3//y_hat = integrate_ode_rk45(dY_funStan, inits, t0, tObs, params, x_r, x_i ,1e-6, 1e-3, 0.8e2 ) ;        //integrate_ode_rk45  // pos use:: 1e-2, 1e-2, 5e3// 1e-4,1e-1,5e3//print(y_hat[1:3,1:5]);
}
model {

// Likelihood of model given data
for(i in 1:nt){for(j in 1:nStatesObs){ lnY[i,j] ~ normal( log(y_hat[i,j]  + 1.0 ) , sigma ); }}//for( j in 1:nStatesObs ){   if ( is_non_zero_state[j]==1 ){  lnY[,j] ~ normal(log( 1.0 + to_row_vector(y_hat[,j]) ) , sigma ) ;   }}                     

// Prior on the rate or proliferation of cells to aid identifiability of other fail pars
r_N_raw ~ normal(1e-5, 1);
r_A_raw ~ normal(1e-5, 1);
k_N ~ normal(1.6e+05, 5e4);  
k_A ~ normal(1.6e+05, 5e4);
delta_N ~ normal(0, 0.5);
delta_A ~ normal(0, 0.5);
r_N ~ normal(1e-4, 1);
r_A ~ normal(1e-4, 1);
B_N_raw ~ normal(0, 2);
B_A_raw ~ normal(0, 2);
delta_N_raw ~ normal(0, 2);
delta_A_raw ~ normal(0, 2);
sigma_raw ~ normal(-3, 2);

}
generated quantities {   
//Generate predicted data over the whole time series and pointwise loglik:
vector[nt*nStatesObs] log_lik;
real pred_I[nPred, nStates];
pred_I = integrate_ode_bdf(dY_funStan, inits, t0, pred_ts, params, x_r, x_i, 1e-2, 1e-1,1e2);//1e-6, 1e-1, 2e3);//1e-4,1e-1,5e3);//

for (i in 1:nt) {
for (j in 1:nStatesObs) {
log_lik[j+nStatesObs*(i-1 )] = normal_lpdf(lnY[i,j] |  log(y_hat[i,j] + 1.0 )  , sigma);
}}

}"

require(rstan)
library("bayesplot")
require(parallel)


test <- stan( model_code = Cancer_mod, data = stan_data, init=0,
              chains = 1, cores =1, seed = 1234,
              iter = 150, warmup=100,
              control = list(
                max_treedepth = 2,
                adapt_delta = 0.2 ) )
# test <- stan( model_code = Cancer_mod, data = stan_data, init=0,
#               chains = 2, cores =2, seed = 1234, 
#               iter = 120, warmup=100,
#               control = list(
#                 max_treedepth = 4,
#                 adapt_delta = 0.6 ) )
# 
# test <- stan( model_code = Cancer_mod, data = stan_data, init=0,
#               chains = 3, cores =3, seed = 1234, 
#               iter = 1500, warmup=1000,
#               control = list(
#                 max_treedepth = 7, 
#                 adapt_delta = 0.6 ) )
# 
# test <- stan( model_code = Cancer_mod, data = stan_data, init=0,
#               chains = 3, cores =3, seed = 1234, 
#               iter = 1500, warmup=1000,
#               control = list(
#                 max_treedepth = 4, 
#                 adapt_delta = 0.6 ) )
# 
# test2 <- stan( model_code = Cancer_mod, data = stan_data, init=0,
#                chains = 3, cores =3, seed = 1234, 
#                iter = 500, warmup=200,
#                control = list(
#                  max_treedepth = 3, 
#                  adapt_delta = 0.4 ) )
# test3 <- stan( model_code = Cancer_mod, data = stan_data, init=0,
#                chains = 3, cores =3, seed = 1234, 
#                iter = 1000, warmup=750,
#                control = list(
#                  max_treedepth = 3, 
#                  adapt_delta = 0.8 ) )



test <- stan( model_code = Cancer_mod, data = stan_data, init=0,
              chains = 3, cores =3, seed = 1234, 
              iter = 2500, warmup=2000,
              control = list(
                max_treedepth = 4, 
                adapt_delta = 0.6 ) )
# test2 <- stan( model_code = Cancer_mod, data = stan_data, init=0,
#                chains = 3, cores = 3, seed = 1234, 
#                iter = 3500, warmup=3000,
#                control = list(
#                  max_treedepth = 4, 
#                  adapt_delta = 0.6 ) )
# 
# expose_stan_functions(test2)
require(loo)
log_lik_1=extract_log_lik(test, merge_chains = F)
r_eff_1=relative_eff(log_lik_1)
loo_1 <- waic(log_lik_1)

#save(test,stan_data,fit.dd,Info,doses,doses_seq,experim_unit,Experiment_Data,inits,lu_table,ndose,obs.times,pred_t,times,x_iSetts,y0,
#     index_N,index_A,file="/Users/jason/Downloads/Ribo_fit_Subclonal compModel.Rdata")
load(file="/Users/jason/Dropbox/Vince data/Fitted Model Params/Bayeian Ribo Params/Life history full/Ribo_fit_LifeHistoryModel.Rdata")
pars_guess <- c(
  r_N=exp(0),
  r_A=exp(0),
  K_q_N=exp(3+0),  # B_N in stan
  K_q_A=exp(3+0), # B_N in stan
  c_N=exp(0),#1/2,
  c_A=exp(0),
  K_N =1 + exp(10 + 0) ,
  K_A =1 + exp(10 + 0) ,
  delta_N=exp(0-4),
  delta_A=exp(0-4),
  gamma_N=exp(0-2),
  gamma_A=exp(0-2),
  delta=exp(0),
  lambda_N =exp(0-3),
  lambda_A =exp(0-3) 
)
xx<-c(pars_guess["r_N"],pars_guess["r_A"],
      pars_guess["K_q_N"],pars_guess["K_q_A"],
      pars_guess["c_N"],pars_guess["c_A"],
      pars_guess["K_N"],pars_guess["K_A"],
      pars_guess["gamma_N"],pars_guess["gamma_A"],
      pars_guess["delta_N"],pars_guess["delta_A"],
      pars_guess["delta"],
      pars_guess["lambda_N"],pars_guess["lambda_A"])

index_N <- which(grepl("N",names(inits))&!grepl("Z_N",names(inits)) )
index_A <- which(grepl("A",names(inits))&!grepl("Z_A",names(inits)) )


Life_History <- function(t, y, pars,doses) {
  with( as.list( c(pars, y) ), {
    N <- y[index_N]
    A <- y[index_A]
    Z_N<- y[index_Z_N]
    Z_A<- y[index_Z_A]
    E <- y[index_Eg]
    RC <- (1 - (N + Z_N )/K_N - (A + Z_A )/K_A) # Resource limitation
    
    q_x_N <- (doses)/(K_q_N + (doses)) # Quiesce
    q_x_A <- doses/(K_q_A + doses) # Quiesce
    Gphase_N <- r_N*(1 + E/(1+c_N*E))*RC # Enter G1/S checkpoint phase
    Gphase_A <- r_A*(1 + E/(1+c_A*E))*RC # Enter G1/S checkpoint phase
    
    d_N <- ( Gphase_N * (1 - q_x_N) - Gphase_N*q_x_N - lambda_N )*N     #/(1+B_N*doses)
    d_A <- ( Gphase_A * (1 - q_x_A) - Gphase_A*q_x_A - lambda_A )*A     #/(1+B_A*doses)
    
    d_Z_N <- (lambda_N + Gphase_N*q_x_N)*N - delta_N*Z_N
    d_Z_A <- (lambda_A + Gphase_A*q_x_A)*A - delta_A*Z_A 
    
    d_E <- gamma_N*N +  gamma_A*A  - delta*E   #+ Z_N + Z_A
    list(c(d_N, d_A,d_Z_N, d_Z_A,d_E))
  } )
}
plot(dY_funStan(0 , inits+100, xx, doses,  x_i=1:length(x_iSetts))
     ,Life_History(0, inits+100, pars_guess,doses)[[1]])
length(Life_History(0, y0, pars_guess,doses)[[1]])

#save(test,stan_data,fit.dd,Info,doses,doses_seq,experim_unit,Experiment_Data,inits,lu_table,ndose,obs.times,pred_t,times,x_iSetts,y0,
#     index_N,index_A,index_Z_N_,index_Z_A_,index_S_N_,index_S_A_,index_E,file="/Users/jason/Downloads/Ribo_fit_LifeHistoryModel.Rdata")

#Life_History(0, y0, pars_guess,doses)[[1]]

stan_trace(test,inc_warmup=F,pars=c("r_N","k_N","B_N" ,"c_N" ,  "gamma_N" ,"delta_N","lambda_N","delta_SN",
                                    "r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta_A","lambda_A","delta_SA","delta","sigma","lp__")
)

stan_trace(test,inc_warmup=T,pars=c("k_N"))
stan_trace(test,inc_warmup=T,pars=c("B_A","B_N"))
stan_trace(test,inc_warmup=T,pars=c("r_A","r_N"))

stan_trace(test,inc_warmup=T,pars=c("lambda_N"))
stan_dens(test,inc_warmup=T,pars=c("delta_SA"))
stan_trace(test,inc_warmup=F,pars=c("lp__"))

stan_dens(test,inc_warmup=F,pars=c("k_N"))
stan_trace(test,inc_warmup=T,pars=c("lp__"))

stan_trace(test,inc_warmup=T,pars=c("gamma_A" ,  "gamma_N" ,"sigma"))

stan_trace(test,inc_warmup=F,pars=c("r_N","k_N","B_N" ))
stan_trace(test,inc_warmup=F,pars=c("r_N","k_N","B_N" ,"delta_N",
                                    "r_A","k_A","B_A" ,"delta_A","sigma")
)
stan_trace(test2,inc_warmup=F,pars=c("r_N","k_N","B_N" ,"c_N" ,  "gamma_N" ,"sigma"))
stan_trace(test2,inc_warmup=F,pars=c("r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta"))
stan_dens(test,inc_warmup=F,pars=c("r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta"))

stan_dens(test,inc_warmup=F,pars=c("r","k","B" ,"c" ,  "gamma" ,"sigma"),separate_chains=T)
stan_ac(test,inc_warmup=F,pars=c("r","k","B" ,"c" ,  "gamma" ,"sigma"))
stan_dens(test,inc_warmup=F,pars=c("c"))
stan_plot(test1,inc_warmup=F,pars=c("r","B" ,"c" ,  "gamma" ,"sigma"))
stan_plot(test1,inc_warmup=F,pars=c("k"))
stan_diag(test)
stan_diag(test2,information='stepsize')
stan_diag(test,information='treedepth')
stan_diag(test2,information='divergence')
#,pars=c("r_N","k_N","B_N" ,"c_N" ,  "gamma_N" ,"delta_N","lambda_N",
#       "r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta_A","lambda_A","delta","sigma")
#)
## we can see by plotting param pairs that the r and gamma are strongle correlated, but with slight prior constraints, we are ok.
pairs(test,pars=c("r_N","k_N","B_N" ,"c_N" )) #,  "gamma_N" ,"sigma"))


# Predicting and then plotting estimates
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
pairs(pars.out,cex=0.5)

preds <- data.table(gather(cbind(post_pred,lp=extract(test,pars="lp__")[[1]] ),key,value,-mcmc_id,-State_numb,-lp,-ode_ID))
preds[,State:=(gsub("[[:digit:]]", "", State_numb, perl = TRUE))]
preds[,experim_unit:= gsub(State, "", State_numb),by=c("ode_ID" ,"mcmc_id","State_numb" )]
preds[,State:=sub("_$","",State)]

#preds[, c("State", "experim_unit") := tstrsplit(State_numb, "_", fixed=TRUE)]
preds[, timeid:= as.numeric(as.character(str_remove(key, "V")))]
preds2 <- data.table(merge(preds, data.table(Day=pred_t,timeid=1:length(pred_t)),by="timeid"))
preds2[,experim_unit:=as.numeric(as.character(experim_unit))]

#[mcmc_id==1][Day==1][experim_unit==1]
preds2_wide <- data.table( preds2 %>%
                             select(-State_numb,-lp,-key,-ode_ID)%>%
                             
                             spread(State,value) ) 
preds2_wide[,Ntot:=N+Z_N+S_N]
preds2_wide[,Atot:=A+Z_A+S_A]
preds2_wide <- preds2_wide%>%select(-c(A,E,N,Z_A,Z_N,S_A,S_N))

setnames(preds2_wide,old=c("Ntot", "Atot"),new=c("N", "A"))
luinfo<-unique(preds2 %>% select(State_numb,lp,key,ode_ID,experim_unit, Day,mcmc_id))

preds3_wide<-merge(preds2_wide,luinfo,by=c("experim_unit", "Day","mcmc_id"))
preds3 <- data.table(gather(preds3_wide, State,value,N : A) )
#preds2[mcmc_id==1][Day==1][experim_unit==1] %>% select(timeid,lpmcmc_id,Day,experim_unit,value, State) %>% spread(State,value)

#preds2$State%>%unique()
abc <- merge(preds3,#preds2,
             unique(fit.dd%>%
                      dplyr::select(experim_unit,DoseNum,Replicate,Compostition,State,Initial_Condition,ode_ID)
             ),
             by=c("ode_ID","State","experim_unit"))
abc$State%>%unique()

#save(abc,preds3,preds3_wide,luinfo,preds2_wide,preds2,preds,test,pars.out,extrct_preds,post_pred,stan_data,fit.dd,Info,doses,doses_seq,experim_unit,Experiment_Data,inits,lu_table,ndose,obs.times,pred_t,times,x_iSetts,y0,
#     index_N,index_A,index_Z_N_,index_Z_A_,index_S_N_,index_S_A_,index_E,file="/Users/jason/Downloads/Fitted_Ribo_fit_LifeHistoryModel.Rdata")


ggplot(fit.dd[value>0], aes(y=log(value),x=Day))+
  #geom_point(aes(group=X,col=Compostition))+
  facet_grid(State~DoseNum)+theme_classic(base_size=21)+
  #geom_point(data=out_long[y0>0&State!="E"],aes(y=log(1+y0) , x=time))+
  
  geom_line(data=abc[value>10][mcmc_id%in% round(seq(1,max(mcmc_id),length=200))] #&mcmc_id%in%c(21:40)]
            ,aes(y=log(value),x=Day,col=Compostition,group=interaction(mcmc_id,ode_ID)),alpha=0.1)+
  geom_point(aes(group=X,col=Compostition))+theme(aspect.ratio = 1)+
  labs(y="cell count (ln)",x="Day")


abc[ ,noised:= exp(rnorm(nrow(abc),log(value),sd=mean(pars.out$sigma)))]


mean_credibleregion<-data.table(abc%>%group_by(ode_ID,State,experim_unit, Day,timeid ,State_numb,
                                               DoseNum ,Replicate, Compostition, Initial_Condition)%>%
                                  summarise(lcl=quantile(value, probs = 0.025),
                                            ucl=quantile(value, probs = 0.975),
                                            lcl_sigma=quantile(noised, probs = 0.05,na.rm=T),
                                            ucl_sigma=quantile(noised, probs = 0.95,na.rm=T),
                                            mu=mean(value)))

ggplot(fit.dd[value>0][State=="N"][DoseNum==0], aes(y=log(value),x=Day))+
  #geom_point(aes(group=X,col=Compostition))+
  facet_grid(State~DoseNum)+theme_classic(base_size=21)+
  #geom_point(data=out_long[y0>0&State!="E"],aes(y=log(1+y0) , x=time))+
  
  geom_ribbon(data=mean_credibleregion[mu>10][State=="N"][DoseNum==0]
              ,aes(y=log(mu),ymax=log(ucl_sigma),ymin=log(lcl_sigma),x=Day,fill=Compostition,group=interaction(ode_ID)),colour =NA,alpha=0.4)+
  labs(y="cell count (ln)",x="Day")+
  geom_point(data=fit.dd[value>0][State=="N"][DoseNum==0],aes(group=X,col=Compostition))+theme(aspect.ratio = 1)+
  geom_line(data=mean_credibleregion[mu>10][State=="N"][DoseNum==0]
            ,aes(y=log(mu),x=Day,col=Compostition,group=interaction(ode_ID)))


fit.dd[,StateLab:="Adapted"]
fit.dd[State=="N",StateLab:="Naive"]
mean_credibleregion[,StateLab:="Adapted"]
mean_credibleregion[State=="N",StateLab:="Naive"]

fit.dd[,DoseLab:=paste0(DoseNum," nM")]

fit.dd$DoseLab <- factor(fit.dd$DoseLab, levels = paste0(unique(fit.dd$DoseNum)," nM"))

mean_credibleregion[,DoseLab:=paste0(DoseNum," nM")]
mean_credibleregion$DoseLab <- factor(mean_credibleregion$DoseLab, levels = paste0(unique(fit.dd$DoseNum)," nM"))


p1<-ggplot(fit.dd[value>0], aes(y=log(value),x=Day))+
  #geom_point(aes(group=X,col=Compostition))+
  facet_grid(StateLab~DoseLab)+theme_classic(base_size=21)+
  #geom_point(data=out_long[y0>0&State!="E"],aes(y=log(1+y0) , x=time))+
  
  geom_ribbon(data=mean_credibleregion[mu>10]
              ,aes(y=log(mu),ymax=log(ucl_sigma),ymin=log(lcl_sigma),x=Day,fill=Compostition,group=interaction(ode_ID)),colour =NA,alpha=0.14)+
  # geom_ribbon(data=mean_credibleregion[mu>10]
  #             ,aes(y=log(mu),ymax=log(ucl),ymin=log(lcl),x=Day,fill=Compostition,group=interaction(ode_ID)),colour =NA,alpha=0.3)+
  labs(y="Cell count (ln)",x="Day")+
  geom_point(data=fit.dd[value>0],aes(group=X,col=Compostition))+theme(aspect.ratio = 1)+
  geom_line(data=mean_credibleregion[mu>10],
            aes(y=log(mu),x=Day,col=Compostition,group=interaction(ode_ID)))

p1 # save 10 x 20 pdf

ggplot(abc[value>10],aes(x=lp,y=mcmc_id))+geom_point()

ggplot(pars.out, aes(x= lp__, y= k_N) ) + geom_point()
ggplot(pars.out, aes(x= r_N, y= k_N) ) + geom_point()
# 
# library(MASS)
# library(cluster)
# 
# fit <- cov.mve(samples, quantile.used = nrow(samples) * 0.75)
# 
# ellipse_boundary <- predict(ellipsoidhull( 
#   as.matrix(M)[,3:4][
#   cov.mve(
#     as.matrix(M)[,3:4], quantile.used = round(nrow(M) * 0.75))$best,]
#   ))
# plot(exp(M[,3:4] ))
# lines(exp(ellipse_boundary), col="lightgreen", lwd=3)
# points((fitted_means[3] ),(fitted_means[4] ),cex=2,col="red",pch=16)

# summarise post
M <- log(pars.out%>%dplyr::select(-lp__))
covMat <- cov( M )
corMat <- cor(M, method="pearson")
fitted_means <- exp(colMeans(M))

#save(mean_credibleregion,fitted_means,M,covMat,file="/Users/jason/Downloads/Summarising posterior for Ribo fit Life history.Rdata")
library(ellipse) # References ellipse()

lines( ellipse( corMat, centre = log(fitted_means)) , col='red')
evals <- eigen(corMat)$values
evecs <- eigen(corMat)$vectors
# Angles of a circle
a <- seq(0, 2*pi, len=100)

# Get critical value
c2 <- qchisq(0.95, 2)
c <- sqrt(c2)


corMat
image(corMat)

covMat



dY_fun <- Life_History#function(t, y, pars, index_N, index_E,doses) {
# with( as.list( c(pars, y) ), {
#  N <- y[index_N]
#  A <- y[index_A]
#  E <- y[index_E]
#  d_N <- r_N*N*(1/(1 + B_N*doses))*(1 + E/(1 + c_N*E))*(1 - N/K_N - A/K_A) - deltaN*N # d_N <- r*N*(1/( 1 + B*doses ))*(1 + E/(1 + c*E))*(1 - N/k) #-deltaN*N
#  d_A <- r_A*A*(1/(1 + B_A*doses))*(1 + E/(1 + c_A*E))*(1 - N/K_N - A/K_A) - deltaA*A    
#  d_E <- gamma_N*N + gamma_A*A - delta*E
#  list(c(d_N, d_A, d_E))
#} )
#}

par.vec_mcmc_i <- as.vector(unlist(pars.out[1]))
names(par.vec_mcmc_i) <-colnames(pars.out)

pars1<-par.vec_mcmc_i 
names(pars1)[names(pars1)=="k_N"] <- "K_N" 
names(pars1)[names(pars1)=="k_A"] <- "K_A" 
names(pars1)[names(pars1)=="B_N"] <- "K_q_N" 
names(pars1)[names(pars1)=="B_A"] <- "K_q_A"
#names(pars1)[names(pars1)=="delta_N"] <- "deltaN" 
#names(pars1)[names(pars1)=="delta_A"] <- "deltaA" 
#pars1["delta"]<-1e-25

# Simulation evaluation times (hours) for 21 days

out_N_N <- data.table(ode(y=y0, parms=pars1, times=times, func=dY_fun,doses=doses))
out_long0 <- data.table(gather(out_N_N,Variable,y0,-c(time)))
out_long0[, c("State","ID") := tstrsplit(Variable, "_", fixed=TRUE)]
out_long0$State <- factor(out_long0$State, levels = c("N","A","E"))
out_long0[,experim_unit:=as.numeric(as.character(ID))]

out_long0[,Variable:=NULL]
out_long <- merge(out_long0,Info,by= c("State","experim_unit"))[order(State,experim_unit,time,DoseNum)]


ggplot(data= out_long, aes( y= log(1+y0), x= time , col= DoseNum,group=experim_unit ) )+
  geom_line()+
  theme_classic()+ labs(y= "State", x= "Time (days)")+
  facet_wrap( ~ State,scales="free")#+#, labeller = as_labeller(facet_names))+ #theme(legend.position="none")


ggplot(data= abc, aes( y= lp, x= mcmc_id , col= DoseNum,group=experim_unit ) )+
  geom_line()+
  theme_classic()+ labs(y= "State", x= "Time (days)")+
  facet_wrap( ~ State,scales="free")#+#, labeller = as_labeller(facet_names))+ #theme(legend.position="none")

# &ode_ID%in%1:1500
ggplot(fit.dd[value>0],aes(y=log(value),x=Day))+
  #geom_point(aes(group=X,col=Compostition))+
  facet_grid(State~DoseNum)+theme_classic(base_size=21)+
  #geom_point(data=out_long[y0>0&State!="E"],aes(y=log(1+y0) , x=time))+
  
  geom_line(data=abc[value>10]#&mcmc_id%in%c(21:40)]
            ,aes(y=log(value),x=Day,col=Compostition,group=interaction(mcmc_id,ode_ID)),alpha=0.1)+
  geom_point(aes(group=X,col=Compostition))+theme(aspect.ratio = 1)+
  labs(y="cell count (ln)",x="Day")





#Model comparison
library("loo")

# Extract pointwise log-likelihood
# using merge_chains=FALSE returns an array, which is easier to 
# use with relative_eff()
log_lik_1 <- extract_log_lik(test2, merge_chains = FALSE)
waic1 <- loo::waic(log_lik_1)
# as of loo v2.0.0 we can optionally provide relative effective sample sizes
# when calling loo, which allows for better estimates of the PSIS effective
# sample sizes and Monte Carlo error
r_eff <- relative_eff(exp(log_lik_1), cores = 2) 

# preferably use more than 2 cores (as many cores as possible)
# will use value of 'mc.cores' option if cores is not specified
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)









ggplot(fit.dd[value>0&ode_ID%in%1:1500],aes(y=log(value),x=Day,group=X,col=Compostition))+
  geom_point()+
  facet_grid(State~DoseNum)+theme_classic()+
  geom_line(data=abc[value>10&mcmc_id%in%c(1:54) &ode_ID%in%1:1500],aes(y=log(value),x=Day,group=interaction(mcmc_id,ode_ID),alpha=0.01))

abc[mcmc_id%in%1&DoseNum==10&State=="N"&Compostition=="polyculture"]

ggplot(fit.dd[value>0&ode_ID%in%1:15],aes(y=log(value),x=Day,group=X,col=Compostition))+
  geom_point()+
  facet_grid(State~DoseNum)+
  geom_line(data=abc[mcmc_id==1&ode_ID%in%1:15],aes(y=log(value),x=Day,group=ode_ID,col=Compostition))



ggplot(fit.dd[value>0&ode_ID%in%1:1500],aes(y=log(value),x=Day,group=X))+#,col=Compostition))+
  geom_point()+
  facet_grid(State~DoseNum)+theme_classic()+
  geom_line(data=abc[lp>quantile(lp,probs=0.75)&value>10&mcmc_id%in%c(1:200) &ode_ID%in%1:1500],aes(y=log(value),x=Day,group=interaction(mcmc_id,ode_ID),col=lp))

hist(abc$lp)

ggplot(fit.dd[value>0&ode_ID%in%1:1500],aes(y=log(value),x=Day,group=X,col=Compostition))+
  geom_point()+
  facet_grid(State~DoseNum)+theme_classic()+
  geom_line(data=abc[value>10&mcmc_id%in%c(1:100) &ode_ID%in%1:1500],aes(y=log(value),x=Day,group=interaction(mcmc_id,ode_ID),col=Compostition))






ggplot(fit.dd[value>0],aes(y=log(value),x=Day,group=X,col=Compostition))+
  geom_point()+
  facet_grid(State~DoseNum)




plot(
  preds2[mcmc_id==1&Day==4&State=="N"]$value ,
  fit.dd[Day==4&State=="N"]$value)


plot(
  log(1+preds2[mcmc_id==100&Day==21&State=="N"]$value) ,
  log(1+fit.dd[Day==21&State=="N"]$value) )
abline(0,1)


i=1:3
ggplot(preds2[Rep%in%i&mcmc_id%in%1:20 &State=="N"][value>0],aes(x=time,y=log(1+value)))+
  geom_line(aes(group=interaction(mcmc_id,Rep),col=lp),alpha=0.15)+theme_classic()+
  geom_point(data=fit.dd[ID%in%i&State=="N"&y0>0],aes(y=log(1+y0),x=time),col="blue")


ggplot(preds2[Rep%in%i&id%in%1:20 &State=="N"][value>0],aes(x=time,y=log(1+value)))+
  geom_line(aes(group=interaction(id,Rep),col=lp),alpha=0.15)+theme_classic()+
  geom_point(data=fit.dd[ID%in%i&State=="N"&y0>0],aes(y=log(1+y0),x=time),col="blue")








# ggplot(fit.dd[],aes(y=log(value),x=Day,group=X,col=DoseNum))+
#   geom_point()+
#   facet_grid(Compostition~State)
# ggplot(fit.dd[],aes(y=log(value),x=Day,group=X,col=Compostition))+
#   geom_point()+
#   facet_grid(DoseNum~State)

#fit.dd <- data.table(out_long[State%in%c("N","A") ][ time%in% (c(3,5,8,12,15,19))] %>%

dplyr::select(State,time,y0,dose,ID))[order(State,ID)]
#fit.dd[,Y:=exp(rnorm(length(y0),log(y0),0.1))]
#fit.dd[,State_num:="State1"]
#fit.dd[State=="A",State_num:="State2"]
#fit.dd[State=="E",State_num:="State3"]
#fit.dd[,ode_id:=paste(State_num,sprintf("%09d", ID) ,sep="_")]
ggplot(data= fit.dd, aes( y= log(1+y0), x= time ) ) +
  geom_line(aes(col= dose ,group=interaction(dose,ID,State))) +
  geom_jitter(aes(col= dose ,group=dose),width=0,height=0.15) +
  theme_classic()+ labs(y= "State", x= "Time (days)") +
  geom_point(data=data.table(y=2000,x=0),aes(y=log(y), x=x))+facet_wrap(~State)









# Observed data reformat
#doses <- sort(unique(Experiment_Data$DoseNum))
#doses <- c(0,100,200)
All_Data <- Experiment_Data[order(Resistance,Day,DoseNum,Replicate)]
days <- sort(unique(All_Data$Day))

obsTable_S <- All_Data %>% dplyr::select(Day,DoseNum,Replicate,Resistance,SCellNum) %>% spread(Replicate,SCellNum )
obsTable_R <- All_Data %>% dplyr::select(Day,DoseNum,Replicate,Resistance,RCellNum) %>% spread(Replicate,RCellNum )
obs_S <- as.matrix(obsTable_S %>% dplyr::select(-c(Day, DoseNum,Resistance)))
obs_R <- as.matrix(obsTable_R %>% dplyr::select(-c(Day, DoseNum,Resistance)))

obsMade_S <- is.finite(as.matrix(log(obs_S)))
obsMade_R <- is.finite(as.matrix(log(obs_R)))

StartSR <- obsTable_S[Day==0 ]
doses <-  StartSR$DoseNum #obsTable_S$DoseNum#
# Initial conditions            #All_Data[Day==0 ] %>%dplyr::select(-c(Area,SCellProp.Image.,SCellProp.FACS.,Dose))
inits1 <- c(Nu=rowMeans( obsTable_S[Day==0 ] %>% dplyr::select(-c(Day, DoseNum,Resistance))), 
            Ni=rep(0,nrow(StartSR) ),
            Au=rowMeans( obsTable_R[Day==0 ] %>% dplyr::select(-c(Day, DoseNum,Resistance))), 
            Ai=rep(0,nrow(StartSR)),
            
            Eg=rep(0,nrow(StartSR) ))

#indexing state variables
index_Nu <- which(grepl("Nu",names(inits1)))
index_Ni <- which(grepl("Ni",names(inits1)))
index_Au <- which(grepl("Au",names(inits1)))
index_Ai <- which(grepl("Ai",names(inits1)))
index_Eg <- which(grepl("Eg",names(inits1)))

#nrow(obsTable_S[Day==0 ])
treatment_lookup <- data.table(Variable=names(inits1),
                               celltype=rep(c("Nu", "Ni","Au", "Ai", "Eg"), each=nrow(StartSR)  ),
                               dose=rep(StartSR$DoseNum,length(c("Nu", "Ni","Au", "Ai", "Eg"))),
                               Resistance=StartSR$Resistance)

