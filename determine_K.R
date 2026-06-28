determine_K <- function(bx, by, byse, num_null, num_junk, K_max = 6, df = 4) {
  
  J <- length(bx)
  
  waic_vals <- rep(NA, K_max)
  
  for (k in 2:K_max) {
    
    eta_uniform <- generate_priors(
      K = k,
      J = J,
      num_null,
      num_junk,
      mode = "uniform"
    )
    
    fit_res <- mcmc_clust(
      bx = bx,
      by = by,
      byse = byse,
      eta = eta_uniform,
      num_null = num_null,
      num_junk = num_junk,
      df = df
    )
    
    waic_vals[k] <- fit_res$waic
  }
  
 
  best_k <- which.min(waic_vals[2:K_max]) + 1
  
  return(list(
    best_K = best_k,
    waic_all = waic_vals[2:K_max]
  ))
}