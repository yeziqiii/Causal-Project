library(rstan)
library(V8)
library(DIRECT)
library(loo)
source("determine_K.R")
source("initialization.R")
mcmc_clust <- function(
    bx, by, byse,
    etas = NULL,
    num_null = 0,
    num_junk = 0,
    Q1 = 0,
    Q2 = 0,
    df = 4
){
  
  
  theta_estimates <- by / bx
  sigma_estimates <- byse / abs(bx)
  
  # Deal with null cluster in the synthetic data
  
  if (num_null == 1 && Q1 > 0) {
    theta_null <- rnorm(Q1, 0, 1)
    sigma_null <- rep(1, Q1)
    
    theta_estimates <- c(theta_null, theta_estimates)
    sigma_estimates <- c(sigma_null, sigma_estimates)
  }
  
  # Deal with junk cluster in the synthetic data
  
  if (num_junk == 1 && Q2 > 0) {
    theta_junk <- rt(Q2, df)
    sigma_junk <- rep(sqrt(2), Q2)
    
    theta_estimates <- c(theta_estimates, theta_junk)
    sigma_estimates <- c(sigma_estimates, sigma_junk)
  }
  
  J <- length(theta_estimates)
  
  
  # If we do not provide eta, then we estimate the number of clusters and use uniform etas
  
  if (is.null(etas)) {
    
    K <- determine_K(bx, by, byse, num_null, num_junk)
    
    etas <- generate_priors(
      K = K,
      J = J,
      num_null = num_null,
      num_junk = num_junk,
      mode = "uniform"
    )
    
  } else {
    
    K <- ncol(etas) - num_null - num_junk
    
    if (nrow(etas) != J) {
      stop("etas must have J rows")
    }
  }
  

  rho <- rep(0, K)
  phi <- rep(1, K)
  
  data_list <- list(
    J = J,
    K = K,
    rho = rho,
    phi = phi,
    theta_estimates = theta_estimates,
    sigma_estimates = sigma_estimates,
    num_null = num_null,
    num_junk = num_junk,
    deg_freedom = df,
    etas = etas
  )
  
  
  model <- stan_model(file = "genetic_variant_modelv2.stan")
  
  base_init <- initialize_centers_weighted(
    theta_estimates,
    sigma_estimates,
    K
  )
  
  init_fun <- function() {
    list(cluster_center = base_init + rnorm(K, 0, 1e-3))
  }
  
  
  fit <- rstan::sampling(
    model,
    init = init_fun,
    data = data_list,
    iter = 2000,
    warmup = 1000,
    chains = 4
  )
  
  
  log_lik <- loo::extract_log_lik(fit, merge_chains = FALSE)
  waic_result <- loo::waic(log_lik, pointwise = TRUE)
  waic_values <- waic_result$estimates["waic", "Estimate"]
  
 
  
  fit_summary <- summary(fit)$summary
  
  center_estimates <- fit_summary[
    paste0("cluster_center[", 1:K, "]"),
    "mean"
  ]
  
  # pi_estimates <- fit_summary[
  #   paste0("pi[", 1:(K + num_null + num_junk), "]"),
  #   "mean"
  # ]
  
  return(list(
    center_estimates = center_estimates,
   # pi_estimates = pi_estimates,
    waic = waic_values
  ))
}