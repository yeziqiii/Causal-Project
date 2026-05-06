source("mcmc_clust.R")
source("data_simulation.R")
source("determine_K.R")
source("generate_uniform_counts.R")
res <- determine_K(bx, by, byse, num_null, num_junk)
print(res)
K <- res$best_K

# counts <- generate_uniform_counts(K, J)
# etas_uniform <- generate_priors(
#   true_cluster_weight = 1 / K,
#   K = K,
#   J = J,
#   counts = counts
# )



mcmc_clust(bx, by, byse, etas = etas_true, num_null = num_null, num_junk = num_junk)

