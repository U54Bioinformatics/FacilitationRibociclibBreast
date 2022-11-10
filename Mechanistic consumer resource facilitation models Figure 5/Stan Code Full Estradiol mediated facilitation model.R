## Stan model of Estradiol mediated facilitation of sensitive by resistant cancer cells under ribociclib treatment
# ODE model structure described in Methods section and Bayesian inference reported in Supplementary information
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