generate_priors <- function(
    K, J,
    num_null = 0,
    num_junk = 0,
    Q1 = NULL,
    Q2 = NULL,
    true_cluster_weight = 0.9,
    a = NULL,
    b = NULL,
    c = NULL,
    true_theta = NULL,
    mode = c("true", "wrong", "uniform")
) {
  
  mode <- match.arg(mode)
  M <- K + num_null + num_junk
  
  is_simulation <- !is.null(a) && !is.null(b) && !is.null(c) && !is.null(true_theta)
  
  if (!is_simulation && mode != "uniform") {
    stop("Real data must use mode 'uniform'")
  }
  
  # Uniform prior for the real data
  if (!is_simulation) {
    return(matrix(1 / M, nrow = J, ncol = M))
  }
  
  if (is.null(Q1)) Q1 <- 0
  if (is.null(Q2)) Q2 <- 0
  
  
  true_centers <- true_theta + c / b
  center_order <- order(true_centers)   
  
  
  labels <- c()
  
  # null
  if (num_null == 1) {
    labels <- c(labels, rep(1, Q1))
  }
  
  # regular clusters
  for (k in seq_along(a)) {
    cluster_id <- which(center_order == k) + num_null
    labels <- c(labels, rep(cluster_id, length(a[[k]])))
  }
  
  # junk
  if (num_junk == 1) {
    labels <- c(labels, rep(M, Q2))
  }
  
  if (length(labels) != J) {
    stop("label length != J")
  }
  
  
  false_weight <- (1 - true_cluster_weight) / (M - 1)
  
  eta <- matrix(false_weight, nrow = J, ncol = M)
  
  for (j in 1:J) {
    eta[j, labels[j]] <- true_cluster_weight
  }
  
  # Random permute regular cluster weights
  
  if (mode == "wrong") {
    
    perm <- seq_len(M)
    
    start <- 1 + num_null
    end   <- M - num_junk
    
    perm[start:end] <- sample(perm[start:end])
    
    eta <- eta[, perm]
  }
  
  return(eta)
}