generate_data <- function(N, p, b, a, true_theta, c, h2_m = 0.5, h2_x = 0.5, h2_y = 0.5) {
  d <- length(a)  # Num of mechanisms
  J <- sum(sapply(a, length))  # Num of genetic variants
  
  
  var_G <- 2 * p * (1 - p)  # Var(S_j) = 2p(1-p)
  Vg_m_list <- numeric(d)
  for (k in 1:d) {
    L_k <- length(a[[k]])
    sum_sq <- sum(a[[k]]^2)
    Vg_m_list[k] <- sum_sq * var_G
  }
  Vg_m <- mean(Vg_m_list)
  
  
  Vg_x <- 0
  for (k in 1:d) {
    Vg_x <- Vg_x + b[k]^2 * Vg_m_list[k]
  }
  
  
  Vg_y <- true_theta^2 * Vg_x
  for (k in 1:d) {
    Vg_y <- Vg_y + c[k]^2 * Vg_m_list[k]
  }
  
  
  Vem <- (1 - h2_m) / h2_m * Vg_m
  Vex <- (1 - h2_x) / h2_x * Vg_x
  Vey <- (1 - h2_y) / h2_y * Vg_y
  
  
  sd_em <- sqrt(Vem)
  sd_ex <- sqrt(Vex)
  sd_ey <- sqrt(Vey)
  
  
  M <- matrix(0, nrow = N, ncol = d)
  G <- matrix(0, nrow = N, ncol = J)
  X <- numeric(N)
  Y <- numeric(N)
  
  for (i in 1:N) {
    
    #G_sample <- rnorm(J, mean = 2*p, sd = sqrt(2*p*(1-p)))  # Approximate the binomial dist
    G_sample <- matrix(rbinom(J, size = 2, prob = p), ncol = J)
    G[i, ] <- G_sample
    
    
    M_sample <- numeric(d)
    start_idx <- 1
    for (k in 1:d) {
      L_k <- length(a[[k]])
      S_k <- G_sample[start_idx:(start_idx + L_k - 1)]
      start_idx <- start_idx + L_k
      
      
      M_sample[k] <- sum(a[[k]] * S_k) + rnorm(1, sd = sd_em)
    }
    M[i, ] <- M_sample
    
    
    genetic_signal_X <- sum(b * M_sample)
    X[i] <- genetic_signal_X + rnorm(1, sd = sd_ex)
    
    
    genetic_signal_Y <- true_theta * X[i] + sum(c * M_sample)
    Y[i] <- genetic_signal_Y + rnorm(1, sd = sd_ey)
  }
  
  return(list(G = G, M = M, X = X, Y = Y))
}


estimate_theta <- function(samples){
  N <- length(samples$X)
  J <- ncol(samples$G)  
  bx <- numeric(J)  
  by <- numeric(J)  
  theta_estimates <- numeric(J)  
  sigma_estimates <- numeric(J)  
  bxse <- numeric(J)  
  byse <- numeric(J)  
  
  X_sample <- samples$X  
  Y_sample <- samples$Y 
  
  for (j in 1:J) {
    G_jth_sample <- samples$G[, j]  
    
    # X ~ G
    X_G <- data.frame(X = X_sample, G = G_jth_sample)
    model_X <- lm(X ~ G, data = X_G)
    
    # Y ~ G
    Y_G <- data.frame(Y = Y_sample, G = G_jth_sample)
    model_Y <- lm(Y ~ G, data = Y_G)
    
    beta_x <- coef(model_X)["G"]
    beta_y <- coef(model_Y)["G"]
    bx[j] <- beta_x
    by[j] <- beta_y
    
    theta_hat <- beta_y / beta_x
    theta_estimates[j] <- theta_hat
    
    
    bxse[j] <- summary(model_X)$coefficients["G", "Std. Error"]
    byse[j] <- summary(model_Y)$coefficients["G", "Std. Error"]
    sigma_estimates[j] <- abs(byse[j] / bx[j])
  }
  
  return(list(
    bx = bx,
    by = by,
    theta_estimates = theta_estimates,
    sigma_estimates = sigma_estimates,
    bxse = bxse,  
    byse = byse   
  ))
}