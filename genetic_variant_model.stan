data {
  int J;  
  int K;  
  vector[K] rho;
  vector[K] phi;
  real<lower = 0> kappa;
  real<lower = 0> ga;
  vector[J] theta_estimates;  
  vector[J] sigma_estimates;  
  int has_null;
  int has_junk;
  real deg_freedom;
  matrix[J, K + has_null + has_junk] eta;
}


parameters {
  ordered[K] cluster_center;
  real<lower=0> alpha;
  simplex[K + has_null + has_junk] pi[J];
}


model {
  int num_of_clusters = K + has_null + has_junk;
  
  cluster_center ~ normal(rho, phi);
  
  alpha ~ gamma(ga, kappa);
  
  
  if (num_of_clusters > 1) {
    for (j in 1:J) {
      vector[num_of_clusters] diri_param = to_vector(eta[j]);
      pi[j] ~ dirichlet(alpha * diri_param);
    }
  }
  
  for (j in 1:J) {
    vector[num_of_clusters] likelihood = log(pi[j]);
    real abs_sigma_hat = abs(sigma_estimates[j]);
    for (num in 1:num_of_clusters) {
      if (has_null == 1){  // There is a null cluster
        if (num == 1){
          likelihood[num] += normal_lpdf(theta_estimates[j] | 0, abs_sigma_hat);
        }
        else{
          if (has_junk == 1){  // There is a junk cluster
            if (num < num_of_clusters){
                likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num - 1], abs_sigma_hat);
            }
            else{
              likelihood[num] += student_t_lpdf(theta_estimates[j] | deg_freedom, theta_estimates[j], abs_sigma_hat);
            }
          }
          else { // There isn’t a junk cluster
            likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num - 1], abs_sigma_hat);
          }
        }
      }
      else{ // There isn't a null cluster
        if (has_junk == 1){ // There is a junk cluster
          if (num < num_of_clusters){
              likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num], abs_sigma_hat);
          }
          else{ // There isn’t a junk cluster
            likelihood[num] += student_t_lpdf(theta_estimates[j] | deg_freedom, theta_estimates[j], abs_sigma_hat);
          }
        }
          
        else {// There isn’t a junk cluster
          likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num], abs_sigma_hat);
        }
      }
    }
  target += log_sum_exp(likelihood);
  }
}



