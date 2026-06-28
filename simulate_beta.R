# Simulate beta data function
simulate_beta <- function(theta_data, N = 1000, tau) {
  
  J <- length(theta_data)
  
  beta_xj_hat_list <- numeric(J)
  beta_yj_hat_list <- numeric(J)
  mu_beta_xj_list <- numeric(J)
  MAF_j_list <- numeric(J)
  se_beta_yj_hat_list <- numeric(J)
  theta_j_hat_list <- numeric(J)
  sigma_j_hat_list <- numeric(J)
  
  for (j in 1:J) {
    # Do calculation for each j
    mu_beta_xj <- rnorm(1, 0, 1)
    MAF_j <- runif(1, 0.05, 0.5)
    beta_xj_hat <- rnorm(1, mu_beta_xj, sqrt(1 / (N * MAF_j * (1 - MAF_j))))
    beta_yj_hat <- rnorm(1, theta_data[j] * beta_xj_hat, sqrt(tau / (N * MAF_j * (1 - MAF_j))))
    theta_j_hat <- beta_yj_hat / beta_xj_hat
    se_beta_yj_hat <- sqrt(1 / (N * MAF_j * (1 - MAF_j)))
    sigma_j_hat <- se_beta_yj_hat / abs(beta_xj_hat)
    
    # Append values to lists
    beta_xj_hat_list[j] <- beta_xj_hat
    beta_yj_hat_list[j] <- beta_yj_hat
    mu_beta_xj_list[j] <- mu_beta_xj
    MAF_j_list[j] <- MAF_j
    se_beta_yj_hat_list[j] <- se_beta_yj_hat
    theta_j_hat_list[j] <- theta_j_hat
    sigma_j_hat_list[j] <- sigma_j_hat
  }
  
  iteration_dict <- list(
    j = 1:J,
    beta_xj_hat = beta_xj_hat_list,
    beta_yj_hat = beta_yj_hat_list,
    mu_beta_xj = mu_beta_xj_list,
    MAF_j = MAF_j_list,
    se_beta_yj_hat = se_beta_yj_hat_list,
    theta_j_hat = theta_j_hat_list,
    sigma_j_hat = sigma_j_hat_list
  )
  
  data_frame <- as.data.frame(iteration_dict)
  return(data_frame)
}