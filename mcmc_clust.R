library(rstan)
library(V8)
library(DIRECT)
library(loo)
source("determine_K.R")
source("initialization.R")

mcmc_clust <- function(
    bx, by, byse,
    eta = NULL,
    num_null = 0,
    num_junk = 0,
    df = 4
){
  theta_estimates <- by / bx
  sigma_estimates <- byse / abs(bx)
  J <- length(theta_estimates)
  
  if (is.null(eta)) {
    cat("--------------------- Estimating K (number of regular clusters) ---------------------")
    K_res <- determine_K(
      bx,
      by,
      byse,
      num_null,
      num_junk
    )
    
    K <- K_res$best_K
    cat("Estimated K =", K, "\n\n")
    eta <- generate_priors(
      K = K,
      J = J,
      num_null = num_null,
      num_junk = num_junk,
      mode = "uniform"
    )
  } else {
    K <- ncol(eta) - num_null - num_junk
    if (nrow(eta) != J) {
      stop("eta must have J rows")
    }
  }
  
  rho <- rep(0, K)
  phi <- rep(1, K)
  gamma <- 1
  kappa <- 2
  
  data_list <- list(
    J = J,
    K = K,
    rho = rho,
    phi = phi,
    gamma = gamma,
    kappa = kappa,
    theta_estimates = theta_estimates,
    sigma_estimates = sigma_estimates,
    num_null = num_null,
    num_junk = num_junk,
    deg_freedom = df,
    eta = eta
  )
  
  model <- stan_model(file = "genetic_variant_modelv2.stan")
  
  base_init <- initialize_centers_eta(
    theta_estimates = theta_estimates,
    sigma_estimates = sigma_estimates,
    eta = eta,
    K = K,
    num_null = num_null,
    num_junk = num_junk
  )
  
  if (length(base_init) != K) {
    stop("Length of base_init must be equal to K.")
  }
  
  init_fun <- function() {
    x <- sort(base_init)
    for (i in 2:K) {
      if (x[i] <= x[i - 1]) {
        x[i] <- x[i - 1] + 1e-6
      }
    }
    
    list(
      cluster_center = x
    )
  }
  
  fit <- rstan::sampling(
    model,
    init = init_fun,
    data = data_list,
    iter = 2000,
    warmup = 1000,
    chains = 4
  )
  
  log_lik <- loo::extract_log_lik(
    fit,
    merge_chains = FALSE
  )
  
  waic_result <- loo::waic(
    log_lik,
    pointwise = TRUE
  )
  
  waic_values <- waic_result$estimates["waic", "Estimate"]
  
  fit_summary <- summary(fit)$summary
  
  center_estimates <- fit_summary[
    paste0("cluster_center[", 1:K, "]"),
    "mean"
  ]
  
  post <- rstan::extract(fit)
  
  pi_estimates <- apply(
    post$pi,
    c(2, 3),
    mean
  )
  
  cluster_assignment <- apply(
    pi_estimates,
    1,
    which.max
  )
  
  return(list(
    center_estimates = center_estimates,
    pi_estimates = pi_estimates,
    cluster_assignment = cluster_assignment,
    waic = waic_values
  ))
}