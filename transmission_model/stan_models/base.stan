data 
{
  int <lower=1> M; 
  int <lower=1> N0; 
  int<lower=1> N[M]; 
  int<lower=1> N2; 
  int deaths[N2, M]; 
  matrix[N2, M] f; 
  int EpidemicStart[M];
  real pop[M];
  int W; 
  int week_index[M,N2];
  real SI[N2]; 
  real AR_SD_MEAN;
  int CasesStart;
  int T2; 
  real par; 
  int phylo_N_len;
  int phylo_N[phylo_N_len]; 
  int phylo_PSamples[phylo_N_len];
  int phylo_NSamples[phylo_N_len];
  matrix[N2, M] PCR_pos_prob; 
  matrix[N2, M] seroconv_cdf; 
  matrix[N2, M] serorev_surv; 
}

parameters 
{
  real<lower=0> R_difference; 
  real<lower=1> y_v1[M];
  real<lower=1> y_v2[M];
  real<lower=0> phi;
  real<lower=0,upper=1> cross;
  real<lower=0> tau;
  real <lower=0, upper=100> ifr1[M];
  real <lower=0> RR[M];
  matrix[W+1,M] weekly_effect;
  real<lower=0, upper=1> weekly_rho;
  real<lower=0, upper=1> weekly_rho1;
  real<lower=0> weekly_sd;
}

transformed parameters 
{
  matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
  matrix[N2, M] prediction = rep_matrix(0,N2,M);
  matrix[N2, M] E_deaths  = rep_matrix(0,N2,M);
  matrix[N2,M] immune = rep_matrix(0,N2,M);
  matrix[N2, M] prediction_v1 = rep_matrix(0,N2,M);
  matrix[N2, M] prediction_v2 = rep_matrix(0,N2,M);
  matrix[N2, M] E_deaths_v1  = rep_matrix(0,N2,M);
  matrix[N2, M] E_deaths_v2  = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_v1 = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_v2 = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_adj_immune_v1 = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_adj_immune_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] cumm_sum_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] cumm_sum_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] immune_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] immune_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] alpha_sus1 = rep_matrix(0,N2,M);
  matrix[N2,M] alpha_sus2 = rep_matrix(0,N2,M);
  matrix[N2,M] n1 = rep_matrix(0,N2,M);
  matrix[N2,M] n2 = rep_matrix(0,N2,M);
  matrix[N2,M] seroconv_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] seroconv_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] pcr_pos_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] pcr_pos_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] E_fraction = rep_matrix(0,N2,M);
  real <lower=0> ifr2[M];
  
  for (m in 1:M)
  {
    for (i in 2:N0)
    {
      cumm_sum_v1[i,m] = cumm_sum_v1[i-1,m] + y_v1[m];
    }
    
    for ( i in ( T2+1 ):( T2+1 ) ) 
    {
      cumm_sum_v2[i,m] = cumm_sum_v2[i-1,m] + y_v2[m];
    }
    
    prediction_v1[1:N0,m] = rep_vector(y_v1[m],N0); 
    prediction_v2[T2:(T2+N0-1),m] = rep_vector(y_v2[m],N0); 
    Rt_v1[,m] = 3.28 * 2 * inv_logit( - weekly_effect[week_index[m],m] );
    Rt_v2[T2:N2,m] = Rt_v1[T2:N2,m] * (R_difference); 
    
    Rt_adj_immune_v1[1:N0,m] = Rt_v1[1:N0,m]; 
    for (i in (N0+1):N2) 
    {
      real convolution_v1 = 1e-15;
      for (j in 1:(i-1)) 
      {
        convolution_v1 += prediction_v1[j, m] * SI[i-j];
        immune_v1[i,m] += prediction_v1[j, m] * exp( - 0.5 * (i-j) * (i-j) / (par * par) );
      }
      if ( i > (T2) ) 
      {
        real convolution_v2 = 1e-15;
        for(j in T2:(i-1))  // start the v2 convolution at T2
        {
          convolution_v2 += prediction_v2[j, m] * SI[i-j];
          immune_v2[i,m] += prediction_v2[j, m] * exp( - 0.5 * (i-j) * (i-j) / ( par * par) );
        }
        alpha_sus2[i,m] = (1 - cross) * immune_v2[i,m] / pop[m];
        n2[i,m] = immune_v2[i,m] + cross * (immune_v1[i,m] * ( 1 - alpha_sus2[i,m] ));
        Rt_adj_immune_v2[i,m] = ( 1 - n2[i,m] / pop[m]) * Rt_v2[i,m]; 
        prediction_v2[i, m] = Rt_adj_immune_v2[i,m] * convolution_v2;
        cumm_sum_v2[i,m]  = cumm_sum_v2[i-1,m] +  prediction_v2[i,m];
      }
      alpha_sus1[i,m] = (1 - cross) * immune_v1[i,m] / pop[m];
      n1[i,m] = immune_v1[i,m] + cross * (immune_v2[i,m] * ( 1 - alpha_sus1[i,m] ));
      
      Rt_adj_immune_v1[i,m] = ( 1 - n1[i,m] / pop[m]) * Rt_v1[i,m];
      prediction_v1[i, m] = Rt_adj_immune_v1[i,m] * convolution_v1;
      cumm_sum_v1[i,m] = cumm_sum_v1[i-1,m] +  prediction_v1[i,m];
      
      cumm_sum[i,m] = cumm_sum_v1[i,m] + cumm_sum_v2[i,m];
      prediction[i, m] = prediction_v1[i, m] + prediction_v2[i, m];
      immune[i, m] = immune_v1[i, m] + immune_v2[i, m];
    }
    E_deaths_v1[1, m]= 1e-15 * prediction_v1[1,m];
    E_deaths_v2[1, m]= 1e-15 * prediction_v2[1,m];
    E_deaths[1, m]= 1e-15 * (prediction_v1[1,m] + prediction_v2[1,m]);
    for (i in 2:N2)
    {
      for(j in 1:(i-1))
      {	
        pcr_pos_v1[i,m] += prediction_v1[j,m] * PCR_pos_prob[i-j,m];
        seroconv_v1[i,m] += prediction_v1[j,m] * seroconv_cdf[i-j,m];
        E_deaths_v1[i,m] += prediction_v1[j,m] * f[i-j,m] * ifr1[m];
      }
      
      if (i > T2)  // start accumulating deaths after T2
      {
        for(j in T2:(i-1)) 
        {
          pcr_pos_v2[i,m] += prediction_v2[j,m] * PCR_pos_prob[i-j,m];
          seroconv_v2[i,m] += prediction_v2[j,m] * seroconv_cdf[i-j,m];
          ifr2[m] = ifr1[m] * RR[m];
          E_deaths_v2[i,m] += prediction_v2[j,m] * f[i-j,m] * ifr2[m];
        }
      }
      E_fraction[i,m] = pcr_pos_v2[i,m]/(pcr_pos_v1[i,m] + pcr_pos_v2[i,m]);	
      E_deaths[i,m] = E_deaths_v1[i,m] + E_deaths_v2[i,m];
    }
  }
}

model 
{
  ifr1 ~ normal(0.32,0.1);
  cross ~ beta(2,1);
  RR ~ lognormal(0,0.5);
  R_difference ~ normal(1,1);
  tau ~ exponential(0.03);
  weekly_sd ~ normal(0, AR_SD_MEAN);
  weekly_rho ~ normal(0.8, 0.05);
  weekly_rho1 ~ normal(0.1, 0.05);
  for (m in 1:M)
  {
    y_v1[m] ~ exponential(1/tau);
    y_v2[m] ~ normal(1,1);
    weekly_effect[3:(W+1), m] ~ normal( weekly_effect[2:W,m]* weekly_rho + weekly_effect[1:(W-1),m]* weekly_rho1,weekly_sd *sqrt(1-pow(weekly_rho,2)-pow(weekly_rho1,2) - 2 * pow(weekly_rho,2) * weekly_rho1/(1-weekly_rho1)));
    
  }
  weekly_effect[2, ] ~ normal(0,weekly_sd *sqrt(1-pow(weekly_rho,2)-pow(weekly_rho1,2) - 2 * pow(weekly_rho,2) * weekly_rho1/(1-weekly_rho1)));
  weekly_effect[1, ] ~ normal(0, 0.01);
  phi ~ normal(0,5);
  for(m in 1:M)
  {
    deaths[EpidemicStart[m]:(N[m]), m] ~ neg_binomial_2(E_deaths[EpidemicStart[m]:(N[m]), m], phi);
    for ( i in 1:phylo_N_len )
    {
      phylo_PSamples[i] ~ binomial(phylo_NSamples[i]+phylo_PSamples[i],E_fraction[phylo_N[i],m]);
    }
  }
}

generated quantities 
{
  real RR_prior = lognormal_rng(0,0.5);
  real cross_prior = beta_rng(2,1);
  real R_difference_prior = normal_rng(1,1);
  real ifr1_prior = normal_rng(0.32,0.1);
}

