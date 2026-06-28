generate_priors <- function(
    K, J,
    num_null = 0,
    num_junk = 0,
    Q1 = NULL,
    Q2 = NULL,
    true_cluster_weight = 0.9,
    wrong_fraction = 0.4,
    wrong_rows = NULL,
    a = NULL,
    b = NULL,
    c = NULL,
    true_theta = NULL,
    mode = c("true", "wrong", "uniform")
) {
  
  mode <- match.arg(mode)
  
  M <- K + num_null + num_junk
  
  is_simulation <- !is.null(a) &&
    !is.null(b) &&
    !is.null(c) &&
    !is.null(true_theta)
  
  if (!is_simulation && mode != "uniform") {
    stop("Real data must use mode 'uniform'.")
  }
  
  # Uniform prior, used for real data or simulation
  if (mode == "uniform") {
    return(matrix(1 / M, nrow = J, ncol = M))
  }
  
  if (!is_simulation) {
    stop("Modes 'true' and 'wrong' require simulation parameters.")
  }
  
  if (is.null(Q1)) Q1 <- 0
  if (is.null(Q2)) Q2 <- 0
  
  true_centers <- true_theta + c / b
  center_order <- order(true_centers)
  
  labels <- c()
  
  # Null cluster label
  if (num_null == 1) {
    labels <- c(labels, rep(1, Q1))
  }
  
  # Regular cluster labels
  # Because cluster_center is ordered in Stan, labels are assigned
  # according to the order of true cluster centers.
  for (k in seq_along(a)) {
    cluster_id <- which(center_order == k) + num_null
    labels <- c(labels, rep(cluster_id, length(a[[k]])))
  }
  
  # Junk cluster label
  if (num_junk == 1) {
    labels <- c(labels, rep(M, Q2))
  }
  
  if (length(labels) != J) {
    stop("label length != J.")
  }
  
  false_weight <- (1 - true_cluster_weight) / (M - 1)
  
  eta <- matrix(false_weight, nrow = J, ncol = M)
  
  # Correct-confident prior by default
  for (j in 1:J) {
    eta[j, labels[j]] <- true_cluster_weight
  }
  
  # Wrong-confident prior:
  # For selected rows, move the largest weight from the true cluster
  # to an incorrect cluster.
  if (mode == "wrong") {
    
    if (is.null(wrong_rows)) {
      wrong_rows <- sample(
        1:J,
        size = ceiling(wrong_fraction * J)
      )
    }
    
    for (j in wrong_rows) {
      
      true_col <- labels[j]
      
      wrong_col <- sample(
        setdiff(1:M, true_col),
        size = 1
      )
      
      eta[j, ] <- false_weight
      eta[j, wrong_col] <- true_cluster_weight
    }
  }
  
  return(eta)
}
diffuse_eta <- function(
    eta,
    new_weight = 0.6,
    diffuse_rows = NULL,
    diffuse_fraction = 0.4
) {
  
  J <- nrow(eta)
  M <- ncol(eta)
  
  if (is.null(diffuse_rows)) {
    diffuse_rows <- sample(
      1:J,
      size = ceiling(diffuse_fraction * J)
    )
  }
  
  eta_new <- eta
  
  for (j in diffuse_rows) {
    
    max_col <- which.max(eta_new[j, ])
    
    false_weight <- (1 - new_weight) / (M - 1)
    
    eta_new[j, ] <- false_weight
    eta_new[j, max_col] <- new_weight
  }
  
  return(eta_new)
}
