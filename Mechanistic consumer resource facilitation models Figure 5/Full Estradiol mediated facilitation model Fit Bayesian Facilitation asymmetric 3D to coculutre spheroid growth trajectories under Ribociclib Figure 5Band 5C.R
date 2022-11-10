rm(list=ls())
require(deSolve) ;require(ggplot2) ; require(tidyr) ; require(dplyr) ; require(data.table) 
require(Matrix); require(fda); library(readr); require(boot);
require(stringr)
require(rstan)
library(bayesplot)
require(parallel)
require(loo)
library(ellipse) 

#Load data 
Experiment_Data <- data.table(read.csv(file="~/Dropbox/Vince data/processed data/riboadavo/Sorted_190711_start_ribociclib_exp_analyzed_8-13-19.csv"))
Experiment_Data[, Day:= Day - min(Day)]
Experiment_Data[, Compostition:= "polyculture"]
Experiment_Data[Resistance %in% c("100% sensitive", "100% resistant"), Compostition:= "monoculture"]

# uniquely label each experimental unit (one microcosm of cells)
experim_unit <- unique(Experiment_Data %>% dplyr::select(DoseNum, Replicate ,    Resistance, Compostition))
experim_unit[, experim_unit:= 1:nrow(experim_unit)]
Experiment_Data <- merge(Experiment_Data,experim_unit,by=c("DoseNum", "Replicate" ,    "Resistance", "Compostition"))
Experiment_Data <- Experiment_Data[order(experim_unit)]  #plot(Experiment_Data$X,Experiment(desiredexp)$X)
plot(Experiment_Data$DoseNum,ylab="dose")
plot(as.numeric(Experiment_Data$Resistance),ylab="response")
plot(Experiment_Data$Replicate,ylab="replicate")

# extract initial conditions and metadata
lu_table<-data.table(Experiment_Data[Day==0] %>% group_by(Resistance, DoseNum) %>% mutate(N=mean(SCellNum), A=mean(RCellNum), "Z_N_"=0, "Z_A_"=0, "S_N_"=0, "S_A_"=0, E=0))
ndose <- length(unique(Experiment_Data$DoseNum))
doses_seq <- unique(Experiment_Data$DoseNum)#rep(0:1,each=1)        

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

# index parameters
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


## Stan model : described in Methods section and reported in Supplementary information
Cancer_mod <- "functions {

real[] dY_funStan(real t, real[] y, real[] params, real[] x_r, int[] x_i) {  
real N; real A;    real E; real Z_N; real Z_A; real S_N; real S_A;     real doses;

real r_N; real B_N; real c_N; real k_N; real gamma_N; real delta_N;  real lambda_N; real delta_SN;
real r_A; real B_A; real c_A; real k_A; real gamma_A; real delta_A;  real lambda_A; real delta_SA;
real delta;
real dydt[ 7*size(x_i) ]; int s_i;  real RC; real q_x_N; real q_x_A; real Gphase_N; real Gphase_A;

s_i = size(x_i);
r_N = params[1]; B_N = params[3];c_N = params[5]; k_N = params[7]; gamma_N = params[9]; delta_N = params[11]; 
r_A = params[2]; B_A = params[4];c_A = params[6]; k_A = params[8]; gamma_A = params[10]; delta_A = params[12];
delta = params[13];
lambda_N = params[14];
lambda_A = params[15];
delta_SN = params[16];
delta_SA = params[17];

// Rate of change functions: dy/dt
for(i in 1:s_i){ 
N = y[i];   A = y[( s_i + i )]; Z_N = y[( 2*s_i + i )]; Z_A = y[( 3*s_i + i )];
S_N = y[( 4*s_i + i )]; S_A = y[( 5*s_i + i )]; E = y[( 6*s_i + i )];
doses = x_r[i];

RC = (1 - (N + Z_N +S_N )/k_N - (A + Z_A + S_A )/k_A);  // Resource limitation
q_x_N = doses/(B_N + doses)   ; // Quiesce
q_x_A = doses/(B_A + doses)   ; // Quiesce

Gphase_N = r_N*(1 + E/(1+c_N*E))*RC ;// Enter G1/S checkpoint phase
Gphase_A = r_A*(1 + E/(1+c_A*E))*RC ;// Enter G1/S checkpoint phase

// Adapted and non-adapted cell state
dydt[i] = ( Gphase_N * (1 - q_x_N) - Gphase_N*q_x_N - lambda_N )*N          ; 
dydt[(s_i + i)] = ( Gphase_A * (1 - q_x_A) - Gphase_A*q_x_A - lambda_A )*A  ; 

// Quiescent adapted and non-adapted cell Be careful that order is consistent
dydt[(2*s_i + i)] =  (lambda_N + Gphase_N*q_x_N)*N - delta_N*Z_N;
dydt[(3*s_i + i)] =  (lambda_A + Gphase_A*q_x_A)*A - delta_A*Z_A ;

//Scenescent
dydt[(4*s_i + i)] = delta_N*Z_N - delta_SN*S_N;
dydt[(5*s_i + i)] = delta_A*Z_A - delta_SA*S_A;

// 6th block of states is E
dydt[(6*s_i + i)] =  gamma_N*N +  gamma_A*A - delta*E;
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
for(i in 1:nt){for(j in 1:nStatesObs){ lnY[i,j] = log(Y[i,j] + 1.0 ); }}  
}
parameters {
real r_N_raw; real k_N_raw; real B_N_raw;  real gamma_N_raw; real delta_N_raw; real delta_SN_raw;
real r_A_raw; real k_A_raw; real B_A_raw;  real gamma_A_raw; real delta_A_raw; real delta_SA_raw;
real c_raw; 

real sigma_raw;
real lambda_N_raw;
real lambda_A_raw;
real delta_raw;                                                       
}
transformed parameters{
real<lower=0> r_N; real<lower=0> k_N; real<lower=0> B_N; real<lower=0> c_N; real<lower=0> gamma_N; real<lower=0> delta_N;real<lower=0> delta_SN;
real<lower=0> r_A; real<lower=0> k_A; real<lower=0> B_A; real<lower=0> c_A; real<lower=0> gamma_A; real<lower=0> delta_A;real<lower=0> delta_SA;
real<lower=0> delta; real<lower=0> lambda_N; real<lower=0> lambda_A;
real<lower=0> sigma; real params [17];          // parameters need to be passed to ODE integrator as a real vector

real y_hat[nt, nStates];         // Output from the ODE solver

r_N= exp(r_N_raw);
k_N= 1 + exp(10 + k_N_raw);
B_N= exp(5+B_N_raw);
c_N= exp(c_raw);
gamma_N= exp(gamma_N_raw-4); 
delta_N= exp(delta_N_raw+1); 
lambda_N=exp(lambda_N_raw-3);
delta_SN= exp(delta_SN_raw+1); 
r_A= exp(r_A_raw);
k_A= 1 + exp(10 + k_A_raw);
B_A= exp(5+B_A_raw);
c_A= exp(c_raw);
gamma_A= exp(gamma_A_raw-4); 
delta_A= exp(delta_A_raw+1); 
lambda_A=exp(lambda_A_raw-3);
delta_SA= exp(delta_SA_raw+1); 
delta= exp(delta_raw);
sigma= exp(sigma_raw);

params[1] = r_N ; params[3] = B_N ; params[5] = c_N ; params[7] = k_N ; params[9] = gamma_N ;  params[11] = delta_N ; 
params[2] = r_A ; params[4] = B_A ; params[6] = c_A ; params[8] = k_A ; params[10] = gamma_A ; params[12] = delta_A ; 
params[13] = delta ; params[14] = lambda_N ; params[15] = lambda_A ;
params[16] = delta_SN ; params[17] = delta_SA ;
//integrate ode
y_hat = integrate_ode_bdf(dY_funStan, inits, t0, tObs, params, x_r, x_i , 1e-2, 1e-1,1e2);
}
model {


// Likelihood of model given data
for(i in 1:nt){for(j in 1:nStatesObs){ lnY[i,j] ~ normal( log(y_hat[i,j] + y_hat[i,j+nStatesObs] + y_hat[i,j+2*nStatesObs] + 1.0 ) , sigma ); }}
// Prior on the rate or proliferation of cells to aid identifiability of other fail pars
r_N_raw ~ normal(1e-5, 1);
r_A_raw ~ normal(1e-5, 1);
k_N ~ normal(1.6e+05, 5e4);  
k_A ~ normal(1.6e+05, 5e4);
lambda_N ~ normal(0, 0.5);
lambda_A ~ normal(0, 0.5);
delta_SN ~ normal(0, 0.5);
delta_SA ~ normal(0, 0.5);
gamma_N ~ normal(0, 0.5);
gamma_A ~ normal(0, 0.5);
delta ~ normal(5, 5);
c_N ~ normal(0, 1);
delta_N ~ normal(0, 0.5);
delta_A ~ normal(0, 0.5);

r_N ~ normal(1e-4, 1);
r_A ~ normal(1e-4, 1);

delta_raw ~ normal(0, 2);
B_N_raw ~ normal(0, 2);
B_A_raw ~ normal(0, 2);
c_raw ~ normal(0, 2);
gamma_N_raw ~ normal(0, 2);
gamma_A_raw ~ normal(0, 2);
delta_N_raw ~ normal(0, 2);
delta_A_raw ~ normal(0, 2);
lambda_N_raw ~ normal(0, 2);
lambda_A_raw ~ normal(0, 2);
delta_SN_raw ~ normal(0, 2);
delta_SA_raw ~ normal(0, 2);
sigma_raw ~ normal(-3, 2);

}
generated quantities {   
//Generate predicted data over the whole time series and pointwise loglik:
vector[nt*nStatesObs] log_lik;
real pred_I[nPred, nStates];
pred_I = integrate_ode_bdf(dY_funStan, inits, t0, pred_ts, params, x_r, x_i, 1e-2, 1e-1,1e2);

for (i in 1:nt) {
for (j in 1:nStatesObs) {
log_lik[j+nStatesObs*(i-1 )] = normal_lpdf(lnY[i,j] |  log(y_hat[i,j] + y_hat[i,j+nStatesObs]+ y_hat[i,j+2*nStatesObs] + 1.0 )  , sigma);
}}

}"

# Fit model to data
test <- stan( model_code = Cancer_mod, data = stan_data, init=0,
              chains = 3, cores =3, seed = 1234, 
              iter = 5000, warmup=4000,
              control = list(
                max_treedepth = 4, 
                adapt_delta = 0.6 ) )

# expose_stan_functions(test)

# Extract log likelihood for WAIC calculation
log_lik_1=extract_log_lik(test, merge_chains = F)
r_eff_1=relative_eff(log_lik_1)
loo_1 <- waic(log_lik_1)

#save(test,stan_data,fit.dd,Info,doses,doses_seq,experim_unit,Experiment_Data,inits,lu_table,ndose,obs.times,pred_t,times,x_iSetts,y0,
#     index_N,index_A,index_Z_N_,index_Z_A_,index_S_N_,index_S_A_,index_E,file="/Users/jason/Downloads/Ribociclib facilitation life history FULL NQS.Rdata")
load("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/Lab Facilitation/Ribociclib facilitation life history FULL NQS.Rdata")

pars_guess <- c(
  r_N=exp(0),
  r_A=exp(0),
  K_q_N=exp(3+0),  # B_N in stan
  K_q_A=exp(3+0), # B_N in stan
  c_N=exp(0),
  c_A=exp(0),
  K_N =1 + exp(10 + 0) ,
  K_A =1 + exp(10 + 0) ,
  delta_N=exp(0-4),
  delta_A=exp(0-4),
  gamma_N=exp(0-2),
  gamma_A=exp(0-2),
  delta=exp(0),
  lambda_N =exp(0-3),
  lambda_A =exp(0-3) ,
  delta_SN= 1,
  delta_SA= 1
)
xx<-c(pars_guess["r_N"],pars_guess["r_A"],
      pars_guess["K_q_N"],pars_guess["K_q_A"],
      pars_guess["c_N"],pars_guess["c_A"],
      pars_guess["K_N"],pars_guess["K_A"],
      pars_guess["gamma_N"],pars_guess["gamma_A"],
      pars_guess["delta_N"],pars_guess["delta_A"],
      pars_guess["delta"],
      pars_guess["lambda_N"],pars_guess["lambda_A"],
      pars_guess["delta_SN"],pars_guess["delta_SA"])

index_N <- which(grepl("N",names(inits))&!grepl("Z_N",names(inits)) &!grepl("S_N",names(inits)) )
index_A <- which(grepl("A",names(inits))&!grepl("Z_A",names(inits)) &!grepl("S_A",names(inits)))
index_Z_N <- which(grepl("Z_N",names(inits)))
index_Z_A <- which(grepl("Z_A",names(inits)))
index_Eg <- which(grepl("E",names(inits)))


# Translation of full model into R ode function
Life_History <- function(t, y, pars,doses) {
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
    Gphase_N <- r_N*(1 + E/(1+c_N*E))*RC # Enter G1/S checkpoint phase
    Gphase_A <- r_A*(1 + E/(1+c_A*E))*RC # Enter G1/S checkpoint phase
    
    d_N <- ( Gphase_N * (1 - q_x_N) - Gphase_N*q_x_N - lambda_N )*N     
    d_A <- ( Gphase_A * (1 - q_x_A) - Gphase_A*q_x_A - lambda_A )*A     
    
    d_Z_N <- (lambda_N + Gphase_N*q_x_N)*N - delta_N*Z_N
    d_Z_A <- (lambda_A + Gphase_A*q_x_A)*A - delta_A*Z_A 
    
    d_S_N <- delta_N*Z_N - delta_SN*S_N
    d_S_A <- delta_A*Z_A - delta_SA*S_A 
    
    d_E <- gamma_N*N +  gamma_A*A  - delta*E   
    list(c(d_N, d_A,d_Z_N, d_Z_A,d_S_N, d_S_A,d_E))
  } )
}
plot(dY_funStan(0 , inits+100, xx, doses,  x_i=1:length(x_iSetts))
     ,Life_History(0, inits+100, pars_guess,doses)[[1]])
length(Life_History(0, y0, pars_guess,doses)[[1]])

#save(test,stan_data,fit.dd,Info,doses,doses_seq,experim_unit,Experiment_Data,inits,lu_table,ndose,obs.times,pred_t,times,x_iSetts,y0,
#     index_N,index_A,index_Z_N_,index_Z_A_,index_S_N_,index_S_A_,index_E,file="/Users/jason/Downloads/Ribo_fit_LifeHistoryModel.Rdata")


## Investigation of Bayesian posteriors
stan_trace(test,inc_warmup=F,pars=c("r_N","k_N","B_N" ,"c_N" ,  "gamma_N" ,"delta_N","lambda_N","delta_SN",
                                    "r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta_A","lambda_A","delta_SA","delta","sigma","lp__")
)
stan_trace(test,inc_warmup=T,pars=c("r_N","k_N","B_N" ,"c_N" ,  "gamma_N" ,"delta_N","lambda_N","delta_SN",
                                    "r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta_A","lambda_A","delta_SA","delta","sigma","lp__")
)
stan_dens(test,inc_warmup=F,pars=c("r_N","k_N","B_N" ,"c_N" ,  "gamma_N" ,"delta_N","lambda_N","delta_SN",
                                   "r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta_A","lambda_A","delta_SA","delta","sigma","lp__")
)

stan_diag(test)
stan_diag(test,information='stepsize')
stan_diag(test,information='treedepth')
stan_diag(test,information='divergence')

## we can see by plotting param pairs that the r and gamma are strongle correlated, but with slight prior constraints, we are ok.
pairs(test) 

# Predicting spheroid states over time and parameter estimates
extrct_preds<- extract(test,pars="pred_I")$ pred_I  
post_pred <- data.table(
  rbindlist(lapply(1:length(inits),function(state_i){
    s_mcmc_t <-extrct_preds[,,state_i] 
    # set names of cols to be time of prediction
    x <- data.table(State_numb= names(y0)[state_i] ,s_mcmc_t)
    x[,ode_ID:=state_i]
    x[,mcmc_id:=1:nrow(x)]
    # gather into long format stack times which are columns of data
  }))
)

pars.out <- data.table(bind_rows(extract(test,pars=c("r_N","k_N","B_N" ,"c_N" ,  "gamma_N" ,"delta_N","lambda_N","delta_SN",
                                                     "r_A","k_A","B_A" ,"c_A" ,  "gamma_A" ,"delta_A","lambda_A","delta_SA","delta","sigma","lp__"))))[lp__ >0]
#pairs(pars.out,cex=0.5)

# Extract the log likelihood estimates
preds <- data.table(gather(cbind(post_pred,lp=extract(test,pars="lp__")[[1]] ),key,value,-mcmc_id,-State_numb,-lp,-ode_ID))
preds[,State:=(gsub("[[:digit:]]", "", State_numb, perl = TRUE))]
preds[,experim_unit:= gsub(State, "", State_numb),by=c("ode_ID" ,"mcmc_id","State_numb" )]
preds[,State:=sub("_$","",State)]
preds[, timeid:= as.numeric(as.character(str_remove(key, "V")))]

# Merge in day information
preds2 <- data.table(merge(preds, data.table(Day=pred_t,timeid=1:length(pred_t)),by="timeid"))
preds2[,experim_unit:=as.numeric(as.character(experim_unit))]

# Restructure and calculate total resistant and sensitive cells across states 
preds2_wide <- data.table( preds2 %>%
                             dplyr::select(-State_numb,-lp,-key,-ode_ID)%>%
                             spread(State,value) ) 
preds2_wide[,Ntot:=N+Z_N+S_N]
preds2_wide[,Atot:=A+Z_A+S_A]
preds2_wide <- preds2_wide %>% dplyr::select(-c(A,E,N,Z_A,Z_N,S_A,S_N))
setnames(preds2_wide,old=c("Ntot", "Atot"),new=c("N", "A"))

# Restructure further by adding in meta data again and recondensing into long format
luinfo <- unique(preds2 %>% select(State_numb, lp, key, ode_ID,experim_unit, Day,mcmc_id))
preds3_wide<-merge(preds2_wide,luinfo,by=c("experim_unit", "Day", "mcmc_id"))
preds3 <- data.table(gather(preds3_wide, State,value, N : A) )

# Merge predicted state and observed state data
abc <- merge(preds3,#preds2,
             unique(fit.dd%>%
                      dplyr::select(experim_unit,DoseNum,Replicate,Compostition,State,Initial_Condition,ode_ID)
             ),
             by=c("ode_ID","State","experim_unit"))
abc$State%>%unique()

#save(abc,preds3,preds3_wide,luinfo,preds2_wide,preds2,preds,test,pars.out,extrct_preds,post_pred,stan_data,fit.dd,Info,doses,doses_seq,experim_unit,Experiment_Data,inits,lu_table,ndose,obs.times,pred_t,times,x_iSetts,y0,
#     index_N,index_A,index_Z_N_,index_Z_A_,index_S_N_,index_S_A_,index_E,file="/Users/jason/Downloads/Fitted_Ribo_fit_LifeHistoryModel.Rdata")
load("/Users/jason/Jason Griffiths Dropbox/jason griffiths/Vince data/Fitted Model Params/Bayeian Ribo Params/Life history full/Fitted_Ribo_fit_LifeHistoryModel.Rdata")
abc[ ,noised:= exp(rnorm(nrow(abc),log(value),sd=mean(pars.out$sigma)))]

# plot observed and expected timecourse
ggplot(fit.dd[value>0], aes(y=log(value),x=Day))+
  facet_grid(State~DoseNum)+theme_classic(base_size=21)+
  geom_line(data=abc[value>10][mcmc_id%in% round(seq(1,max(mcmc_id),length=200))] 
            ,aes(y=log(value),x=Day,col=Compostition,group=interaction(mcmc_id,ode_ID)),alpha=0.1)+
  geom_point(aes(group=X,col=Compostition))+theme(aspect.ratio = 1)+
  labs(y="cell count (ln)",x="Day")

# get mean and credicble regions
mean_credibleregion<-data.table(abc%>%group_by(ode_ID,State,experim_unit, Day,timeid ,State_numb,
                                               DoseNum ,Replicate, Compostition, Initial_Condition)%>%
                                  summarise(lcl=quantile(value, probs = 0.025),
                                            ucl=quantile(value, probs = 0.975),
                                            lcl_sigma=quantile(noised, probs = 0.05,na.rm=T),
                                            ucl_sigma=quantile(noised, probs = 0.95,na.rm=T),
                                            mu=mean(value)))

ggplot(fit.dd[value>0][State=="N"][DoseNum==0], aes(y=log(value),x=Day))+
  facet_grid(State~DoseNum)+theme_classic(base_size=21)+
  geom_ribbon(data=mean_credibleregion[mu>10][State=="N"][DoseNum==0]
              ,aes(y=log(mu),ymax=log(ucl_sigma),ymin=log(lcl_sigma),x=Day,fill=Compostition,group=interaction(ode_ID)),colour =NA,alpha=0.4)+
  labs(y="cell count (ln)",x="Day")+
  geom_point(data=fit.dd[value>0][State=="N"][DoseNum==0],aes(group=X,col=Compostition))+theme(aspect.ratio = 1)+
  geom_line(data=mean_credibleregion[mu>10][State=="N"][DoseNum==0]
            ,aes(y=log(mu),x=Day,col=Compostition,group=interaction(ode_ID)))

# Edit the labels to enhance the plotting 
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

p0 <- ggplot(fit.dd[value>0], aes(y=log(value),x=Day,group=DoseLab))+
  facet_grid(StateLab2~Compostition2)+
  theme_classic(base_size=21)+
  geom_ribbon(data=mean_credibleregion[mu>10]
              ,aes(y=log(mu),ymax=log(ucl_sigma),ymin=log(lcl_sigma),x=Day,fill=DoseNum,group=interaction(ode_ID,DoseLab)),colour =NA,alpha=0.1)+
  labs(y="Cell count (ln)",x="Day")+
  geom_point(data=fit.dd[value>0],aes(group=X,col=DoseNum))+theme(aspect.ratio = 1)+
  geom_line(data=mean_credibleregion[mu>10],
            aes(y=log(mu),x=Day,col=DoseNum,group=interaction(ode_ID)))+
  scale_fill_gradientn(name= "Dose", colours=rainbow(100,end =0.8))+
  scale_color_gradientn(name= "Dose",colours=rainbow(100,end =0.8))


p1<-ggplot(fit.dd[value>0], aes(y=log(value),x=Day))+
  facet_grid(StateLab~DoseLab)+theme_classic(base_size=24)+
  geom_ribbon(data=mean_credibleregion[mu>10]
              ,aes(y=log(mu),ymax=log(ucl_sigma),ymin=log(lcl_sigma),x=Day,fill=Compostition,group=interaction(ode_ID)),colour =NA,alpha=0.14)+
  labs(y="Cell count (ln)",x="Day")+
  geom_point(data=fit.dd[value>0],aes(group=X,col=Compostition))+theme(aspect.ratio = 1)+
  geom_line(data=mean_credibleregion[mu>10],
            aes(y=log(mu),x=Day,col=Compostition,group=interaction(ode_ID)))
p1 # save 10 x 20 pdf

## Exploration of posterior parameter estimates
ggplot(abc[value>10],aes(x=lp,y=mcmc_id))+geom_point()
ggplot(pars.out, aes(x= lp__, y= k_N) ) + geom_point()
ggplot(pars.out, aes(x= r_N, y= k_N) ) + geom_point()
ggplot(pars.out, aes(x= 1/k_A, y= 1/k_N) ) + geom_point()

# Calc ratio of resistant and sensitive cell parameter values
pars.out[,r_RbyN:=r_A/r_N ]
pars.out[,gamma_RbyN:=gamma_A/gamma_N ]
pars.out[,delta_RbyN:=delta_A/delta_N ]
pars.out[,lambda_RbyN:=lambda_A/lambda_N ]
pars.out[,delta_SRbySN:=delta_SA/delta_SN ]
pars.out[,B_RbyN:=(1/B_A)/(1/B_N) ]
pars.out[,k_RbyN:=(k_A)/(k_N) ]

pars.outlong <- data.table(gather(pars.out,par,val,r_RbyN:k_RbyN))
pars.outlong$par <- factor(pars.outlong$par, levels = rev(c("gamma_RbyN","k_RbyN","r_RbyN","lambda_RbyN","B_RbyN","delta_SRbySN","delta_RbyN")))
ggplot(pars.outlong, aes(y=log(val),x= par ,col=par,fill=par) ) + geom_violin()+theme_classic()

fig5cData<- pars.outlong[!par%in%c("delta_SRbySN","delta_RbyN")]
ggplot(pars.outlong[!par%in%c("delta_SRbySN","delta_RbyN")], aes(y=log(val),x= par ,col=par,fill=par) ) + geom_violin(scale="width",bw=.015)+theme_classic(base_size = 19)+
  theme(aspect.ratio=1,legend.position = "none")+coord_flip()+
  labs(x="Process",y="Resistant cell performance \n (relative to sensitive cells)")+
  geom_hline(yintercept=0,linetype=2)+
  scale_x_discrete(labels=rev(c("Facilitation \n contribution","Competition \n effect","Division \n (Baseline)","Quiescence \n (Baseline)","Drug sensitivity")))+
  scale_y_continuous(breaks=log(c(0.125,0.25,0.5,1,2,4,8)),labels=c(0.125,0.25,0.5,1,2,4,8))

ggplot(pars.out, aes(x= gamma_A/gamma_N ) ) + geom_histogram()+theme_classic()
ggplot(pars.out, aes(x= r_A/r_N ) ) + geom_histogram()+theme_classic()
ggplot(pars.out, aes(x= delta_SA/delta_SN ) ) + geom_histogram()+theme_classic()
ggplot(pars.out, aes(x= lambda_A/lambda_N ) ) + geom_histogram()+theme_classic()
ggplot(pars.out, aes(x= delta_A/delta_N ) ) + geom_histogram()+theme_classic()
ggplot(pars.out, aes(x= B_A/B_N ) ) + geom_histogram()+theme_classic()
ggplot(pars.out, aes(x= k_A/k_N) ) + geom_histogram()+theme_classic()
ggplot(pars.out, aes(x= c_N) ) + geom_histogram()+theme_classic()
ggplot(pars.out, aes(x= gamma_A ) ) + geom_histogram()+ geom_histogram(aes(x= gamma_N ))

# summarise post
M <- log(pars.out%>%dplyr::select(r_N:delta,sigma))#-lp__))
covMat <- cov( M )
corMat <- cor(M, method="pearson")
fitted_means <- exp(colMeans(M))
#save(mean_credibleregion,fitted_means,M,covMat,file="/Users/jason/Downloads/Summarising posterior for Ribo fit Life history.Rdata")

plot(0,0,pch="",xlim=c(-2,2),ylim=c(-2,2))
plot( ellipse( corMat, centre = log(fitted_means)) , col='red')
evals <- eigen(corMat)$values
evecs <- eigen(corMat)$vectors
# Angles of a circle
a <- seq(0, 2*pi, len=100)
# Get critical value
c2 <- qchisq(0.95, 2)
c <- sqrt(c2)
corMat
image(corMat)

dY_fun <- Life_History
par.vec_mcmc_i <- as.vector(unlist(pars.out[1]))
names(par.vec_mcmc_i) <-colnames(pars.out)
pars1<-par.vec_mcmc_i 
names(pars1)[names(pars1)=="k_N"] <- "K_N" 
names(pars1)[names(pars1)=="k_A"] <- "K_A" 
names(pars1)[names(pars1)=="B_N"] <- "K_q_N" 
names(pars1)[names(pars1)=="B_A"] <- "K_q_A"

# Simulation evaluation times (hours) for 21 days
index_Z_N<-index_Z_N_
index_Z_A<-index_Z_A_
index_S_N<-index_S_N_
index_S_A<-index_S_A_
index_Eg<-index_E
length( dY_fun(0,y0,pars1,doses ) [[1]])
out_N_N <- data.table(ode(y=y0, parms=pars1, times=times, func=dY_fun,doses=doses))
out_long0 <- data.table(gather(out_N_N,Variable,y0,-c(time)))
out_long0[,State:="E"]
out_long0[Variable%in%c(paste0("N_",1:length(index_N) )),State:="N"];out_long0[Variable%in%c(paste0("Z_N_",1:length(index_N) )),State:="Z_N"];out_long0[Variable%in%c(paste0("S_N_",1:length(index_N) )),State:="S_N"]
out_long0[Variable%in%c(paste0("A_",1:length(index_N) )),State:="A"];out_long0[Variable%in%c(paste0("Z_A_",1:length(index_A) )),State:="Z_A"];out_long0[Variable%in%c(paste0("S_A_",1:length(index_A) )),State:="S_A"]
out_long0[, experim_unit := as.numeric(str_extract(Variable, "([0-9]+)"))] 
out_long0$State <- factor(out_long0$State, levels = c("N","A","Z_N", "Z_A", "S_N", "S_A","E"))
out_long0[,experim_unit:=as.numeric(as.character(experim_unit))]
out_long0[,Variable:=NULL]
Info[, isNpresent:=TRUE];Info[, isApresent:=TRUE]
Info[experim_unit%in%which(y0[index_N]==0), isNpresent:=FALSE]
Info[experim_unit%in%which(y0[index_A]==0), isApresent:=FALSE]
out_long <- merge(out_long0,Info,by= c("State","experim_unit"))[order(State,experim_unit,time,DoseNum)]

ggplot(data= out_long, aes( y= log(1+y0), x= time , col= experim_unit,group=experim_unit ,linetype=interaction(isNpresent,isApresent)) )+
  geom_line()+
  theme_classic()+ labs(y= "State", x= "Time (days)")+
  facet_wrap(Compostition ~ State,scales="free")#+#, labeller = as_labeller(facet_names))+ #theme(legend.position="none")

#### 
#save(index_N, index_A, index_Z_N, index_Z_A, index_S_N, index_S_A, index_Eg, times, y0, doses, Info, dY_fun, Life_History, pars.out, 
#     file="/Users/jason/Downloads/parameterised LH full pars and core data for simulation studies.RData")
load(file="/Users/jason/Dropbox/Cancer_pheno_evo/data/Lab Facilitation/parameterised LH full pars and core data for simulation studies.RData")     
save( fig5cData, file="/Users/jason/Dropbox/Cancer_pheno_evo/data/Lab Facilitation/Fig5Cdata.RData")

#Model comparison test for this model. This approach will be applied to compare all fitted models in WAIC analysis script
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
