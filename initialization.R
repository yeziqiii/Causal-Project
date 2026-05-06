initialize_centers_weighted <- function(theta, sigma, K, jitter_sd = 1e-3) {
  w <- 1 / (sigma^2 + 1e-12)
  probs <- seq(0, 1, length.out = K + 2)[2:(K + 1)]
  centers <- as.numeric(weighted_quantile(theta, w, probs))
  
  # Aligh with ordered center in the stan model
  for (k in 2:K) {
    if (centers[k] <= centers[k - 1]) {
      centers[k] <- centers[k - 1] + 1e-4
    }
  }
  
  centers
}

weighted_quantile <- function(x, w, probs) {
  o <- order(x)
  x <- x[o]
  w <- w[o] / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) {
    idx <- which(cw >= p)[1]
    x[idx]
  })
}