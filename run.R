source("mcmc_clust.R")
source("data_simulation.R")
source("determine_K.R")
library(mrclust)

res <- determine_K(bx, by, byse, num_null, num_junk, K_max = 5)
print(res)
K <- res$best_K

# counts <- generate_uniform_counts(K, J)
# etas_uniform <- generate_priors(
#   true_cluster_weight = 1 / K,
#   K = K,
#   J = J,
#   counts = counts
# )

res_em <- mr_clust_em_fixed_count(
  theta = theta_estimates,
  theta_se = sigma_estimates,
  bx = bx,
  by = by,
  bxse = bxse,
  byse = byse,
  num_substantive_clusters = 2
)
cat("====================================\n")
cat("MR-Clust results\n")
cat("====================================\n")

mr_centers <- sort(
  unique(
    round(res_em$results$all$cluster_mean, 3)
  )
)

print(mr_centers)
mcmc_clust(bx, by, byse, eta = eta_wrong, num_null = num_null, num_junk = num_junk, true_labels = true_labels)
mcmc_clust(bx, by, byse, eta = eta_true, num_null = num_null, num_junk = num_junk, true_labels = true_labels)

