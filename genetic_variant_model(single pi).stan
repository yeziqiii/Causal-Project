data {
  int J;  
  int K;  
  real rho;
  real phi;

  vector[J] theta_estimates;  
  vector[J] sigma_estimates;  

  int num_null;
  int num_junk;

  real deg_freedom;
}

parameters {
  real cluster_center;  
  real<lower = 0> alpha;
  simplex[K + num_null + num_junk] pi;
}

model {

  int num_of_clusters = 1 + num_null + num_junk;

  cluster_center ~ normal(rho, phi);
  alpha ~ lognormal(0, 0.2);

  pi ~ dirichlet(alpha * rep_vector(1, num_of_clusters));

  for (j in 1:J) {

    vector[num_of_clusters] likelihood = log(pi);
    real abs_sigma_hat = abs(sigma_estimates[j]);

    for (num in 1:num_of_clusters) {

      // null
      if (num_null == 1 && num == 1) {

        likelihood[num] += normal_lpdf(theta_estimates[j] | 0, abs_sigma_hat);

      }

      // junk
      else if (num_junk == 1 && num == num_of_clusters) {

        likelihood[num] += student_t_lpdf(theta_estimates[j] | deg_freedom, 0, 1);

      }

      // main cluster
      else {

        likelihood[num] += normal_lpdf(
          theta_estimates[j] | cluster_center,
          abs_sigma_hat
        );
      }
    }

    target += log_sum_exp(likelihood);
  }
}

generated quantities {
  vector[J] log_lik;
  int num_of_clusters = 1 + num_null + num_junk;

  for (j in 1:J) {

    vector[num_of_clusters] likelihood = log(pi);
    real abs_sigma_hat = abs(sigma_estimates[j]);

    for (num in 1:num_of_clusters) {

      if (num_null == 1 && num == 1) {

        likelihood[num] += normal_lpdf(theta_estimates[j] | 0, abs_sigma_hat);

      } else if (num_junk == 1 && num == num_of_clusters) {

        likelihood[num] += student_t_lpdf(theta_estimates[j] | deg_freedom, 0, 1);

      } else {

        likelihood[num] += normal_lpdf(
          theta_estimates[j] | cluster_center,
          abs_sigma_hat
        );
      }
    }

    log_lik[j] = log_sum_exp(likelihood);
  }
}


