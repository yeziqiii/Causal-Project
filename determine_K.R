determine_K <- function(bx, by, byse, num_null, num_junk, K_max = 5) {
  
  J <- length(bx)
  
  
  waic_vals <- rep(NA, K_max)
  
  for (k in 2:K_max) {
    
    etas_uniform <- generate_priors(
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
      etas = etas_uniform,
      K = k,
      num_null,
      num_junk
    )
    
    waic_vals[k] <- fit_res$waic
  }
  
 
  best_k <- which.min(waic_vals[2:K_max]) + 1
  
  return(list(
    best_K = best_k,
    waic_all = waic_vals[2:K_max]
  ))
}