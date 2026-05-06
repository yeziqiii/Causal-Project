experiments <- function(
    N_list, p, a, b, c, true_theta,
    Q1, Q2, df,
    etas,
    rep,
    outputfile,
    raw_outputfile
) {
  
  
  
  J_regular <- sum(sapply(a, length))
  J <- J_regular + Q1 + Q2
  
  num_null <- as.integer(Q1 > 0)
  num_junk <- as.integer(Q2 > 0)
  
  K <- ncol(etas) - num_null - num_junk
  
  
  true_centers <- true_theta + c / b
  ord <- order(true_centers)
  true_centers <- true_centers[ord]
  
  total_cols <- 2 + 4 * K
  
  results <- matrix(
    NA,
    nrow = length(N_list),
    ncol = total_cols
  )
  
  colnames(results) <- c(
    "N",
    paste0("mean_center", 1:K),
    paste0("var_center", 1:K),
    paste0("MSE", 1:K),
    paste0("true_center", 1:K),
    "WAIC"
  )
  
  
  
  raw_results <- data.frame()
  row_idx <- 1
  
  
  for (N in N_list) {
    
    cat("Running N =", N, "\n")
    
    center_estimates <- matrix(
      NA,
      nrow = rep,
      ncol = K
    )
    
    waic_values <- numeric(rep)
    
    for (r in 1:rep) {
      
      cat("  rep =", r, "\n")
      
      samples <- generate_data(
        N, p, b, a, true_theta, c,
        0.95, 0.95, 0.95
      )
      
      estimates <- estimate_theta(samples)
      
      bx   <- estimates$bx
      by   <- estimates$by
      byse <- estimates$byse
      
      fit_res <- mcmc_clust(
        bx = bx,
        by = by,
        byse = byse,
        etas = etas,
        num_null = num_null,
        num_junk = num_junk,
        Q1 = Q1,
        Q2 = Q2,
        df = df
      )
      
      center_estimates[r, ] <- fit_res$center_estimates
      waic_values[r] <- fit_res$waic
      
      
      raw_row <- data.frame(
        N = N,
        rep = r
      )
      
      for (j in 1:K) {
        raw_row[[paste0("center", j)]] <-
          fit_res$center_estimates[j]
      }
      
      raw_results <- rbind(raw_results, raw_row)
    }
    
    
    
    waic_mean <- mean(waic_values)
    
    mse_per_center <- sapply(1:K, function(j) {
      mean(
        (center_estimates[, j] - true_centers[j])^2
      )
    })
    
    results[row_idx, ] <- c(
      N,
      signif(colMeans(center_estimates), 3),
      signif(apply(center_estimates, 2, var), 3),
      signif(mse_per_center, 3),
      true_centers,
      signif(waic_mean, 4)
    )
    
    row_idx <- row_idx + 1
  }
  
  
  # Save summary results
  
  write.table(
    results,
    file = outputfile,
    sep = ",",
    row.names = FALSE
  )
  
  
  # Save raw estimates
  
  write.csv(
    raw_results,
    file = raw_outputfile,
    row.names = FALSE
  )
}