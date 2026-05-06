source("generate_data2.R")
source("generate_priors.R")
source("density_plot.R")


#N_list <- seq(from = 2000, to = 60000, length.out = 30)
N_list <- list(100, 500, 1000, 3000, 5000, 10000, 20000, 50000, 100000)
#N_list <- list(5000)
p <- 0.3

a <- list(c(rep(1,40)),c(rep(-1.4,50)))
b <- c(1, 1)
true_theta <- c(-2)
c <- c(4.2, -6)
K <- length(a)
J_regular <- sum(sapply(a, length))
Q1 <- 0
Q2 <- 0
J <- Q1 + Q2 + J_regular
df <- 4
num_null <- (Q1 > 0)
num_junk <- (Q2 > 0)

counts <- unlist(lapply(a, length))
if (Q1 != 0) {
  counts <- c(Q1, counts)
}

if (Q2 != 0) {
  counts <- c(counts, Q2)
}

rho <- rep(0, K)
phi <- rep(1, K)

num_of_clusters <- K + as.integer(num_null) + as.integer(num_junk)
rep <- 100

# samples <- generate_data(50000, p, b, a, true_theta, c, h2_m = 0.5, h2_x = 0.5, h2_y = 0.5)
# estimates <- estimate_theta(samples)
# bx <- estimates$bx
# by <- estimates$by
# byse <- estimates$byse



etas_true <- generate_priors(
  K = K,
  J = J,
  num_null = num_null,
  num_junk = num_junk,
  Q1 = Q1,
  Q2 = Q2,
  true_cluster_weight = 0.95,
  a = a,
  b = b,
  c = c,
  true_theta = true_theta,
  mode = "true"
)


etas_uniform <- generate_priors(
  K = K,
  J = J,
  num_null = num_null,
  num_junk = num_junk,
  mode = "uniform"
)


etas_wrong <- generate_priors(
  K = K,
  J = J,
  num_null = num_null,
  num_junk = num_junk,
  Q1 = Q1,
  Q2 = Q2,
  true_cluster_weight = 0.95,
  a = a,
  b = b,
  c = c,
  true_theta = true_theta,
  mode = "wrong"
)

