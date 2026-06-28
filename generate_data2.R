generate_data <- function(
    N,
    p,
    b,
    a,
    true_theta,
    c,
    h2_m = 0.5,
    h2_x = 0.5,
    h2_y = 0.5,
    rho = 0.2
) {
  
  D <- length(a)                        # number of mechanisms
  J <- sum(sapply(a, length))           # number of genetic variants
  
  var_G <- 2 * p * (1 - p)              # Gj ~ Binomial(2, p)
  
  # Genetic variance components for mechanisms
  V_M_g <- numeric(D)
  
  for (i in 1:D) {
    V_M_g[i] <- var_G * sum(a[[i]]^2)
  }
  
  # Residual variance for mechanisms
  V_M_noise <- ((1 - h2_m) / h2_m) * V_M_g
  
  # Genetic variance component for X
  V_X_g <- 0
  
  for (i in 1:D) {
    V_X_g <- V_X_g + b[i]^2 * V_M_g[i]
  }
  
  # Residual variance for X
  V_X_noise <- ((1 - h2_x) / h2_x) * V_X_g
  
  # Genetic variance component for Y
  V_Y_g <- true_theta^2 * V_X_g
  
  for (i in 1:D) {
    V_Y_g <- V_Y_g + c[i]^2 * V_M_g[i]
  }
  
  # Residual variance for Y
  V_Y_noise <- ((1 - h2_y) / h2_y) * V_Y_g
  
  # Standard deviations
  sd_M_noise <- sqrt(V_M_noise)
  sd_X_noise <- sqrt(V_X_noise)
  sd_Y_noise <- sqrt(V_Y_noise)
  
  M <- matrix(0, nrow = N, ncol = D)
  G <- matrix(0, nrow = N, ncol = J)
  
  X <- numeric(N)
  Y <- numeric(N)
  
  # Generate samples
  for (n in 1:N) {
    
    # Generate genetic variants
    G_sample <- rbinom(
      J,
      size = 2,
      prob = p
    )
    
    G[n, ] <- G_sample
    
    # Generate mechanisms
    M_sample <- numeric(D)
    
    start_idx <- 1
    
    for (i in 1:D) {
      
      L_i <- length(a[[i]])
      
      s_i <- G_sample[
        start_idx:(start_idx + L_i - 1)
      ]
      
      start_idx <- start_idx + L_i
      
      M_sample[i] <- sum(a[[i]] * s_i) +
        rnorm(1, mean = 0, sd = sd_M_noise[i])
    }
    
    M[n, ] <- M_sample
    
    # Generate correlated residuals for X and Y
    z1 <- rnorm(1)
    z2 <- rnorm(1)
    
    epsilon_X <- sd_X_noise * z1
    
    epsilon_Y <- rho * sd_Y_noise * z1 +
      sqrt(1 - rho^2) * sd_Y_noise * z2
    
    # Generate exposure
    X_signal <- sum(b * M_sample)
    
    X[n] <- X_signal + epsilon_X
    
    # Generate outcome
    Y_signal <- true_theta * X[n] +
      sum(c * M_sample)
    
    Y[n] <- Y_signal + epsilon_Y
  }
  
  return(list(
    G = G,
    M = M,
    X = X,
    Y = Y
  ))
}

estimate_theta <- function(samples) {
  
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
    
    X_G <- data.frame(
      X = X_sample,
      G = G_jth_sample
    )
    
    model_X <- lm(
      X ~ G,
      data = X_G
    )
    
    Y_G <- data.frame(
      Y = Y_sample,
      G = G_jth_sample
    )
    
    model_Y <- lm(
      Y ~ G,
      data = Y_G
    )
    
    # Regression coefficients
    beta_x <- coef(model_X)["G"]
    beta_y <- coef(model_Y)["G"]
    
    bx[j] <- beta_x
    by[j] <- beta_y
    
    # Ratio estimate
    theta_hat <- beta_y / beta_x
    
    theta_estimates[j] <- theta_hat
    
    # Standard errors
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