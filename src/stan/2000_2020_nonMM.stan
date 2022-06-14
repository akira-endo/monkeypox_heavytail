Model4 = "
functions{
// discrete weiull distr.
 real discrete_weibull_lpdf(real y, real x, real alpha, real gamma){
    real lp =  y * log(exp(-gamma*(x^alpha -1)/alpha) - exp(-gamma*((x+1)^alpha -1)/alpha));
    return lp;
  }
  
}
data {
int N;
real<lower=0> x_contact[N];
real<lower=0> y_case00[N];
}
parameters{
real<lower=0> alpha;
real<lower=0> beta;
}

model{ 
alpha ~ normal(0,10);
beta ~ normal(0,10);
for(n in 1:N)
target += discrete_weibull_lpdf(y_case00[n] | x_contact[n], alpha, beta) ;//1/log(N) * discrete_weibull_lpdf(y_case00[n] | x_contact[n], alpha, beta);
}

generated quantities{
vector[N] llk;
for(n in 1:N)
llk[n] = discrete_weibull_lpdf(y_case00[n] | x_contact[n], alpha, beta);
}
"