initialize_centers_eta <- function(
    theta_estimates,
    sigma_estimates,
    eta,
    K,
    num_null = 0,
    num_junk = 0
) {
  # Find cols for the substantive cluster
  start_col <- 1 + num_null
  end_col <- num_null + K
  
  eta_regular <- eta[, start_col:end_col, drop = FALSE]
  
  centers <- numeric(K)
  
  for (k in 1:K) {
    w <- eta_regular[, k] / sigma_estimates^2
    centers[k] <- sum(w * theta_estimates) / sum(w)
  }
  
  centers <- sort(centers)
  
  return(centers)
}