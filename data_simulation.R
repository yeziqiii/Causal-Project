source("generate_data2.R")
source("generate_priors.R")
source("density_plot.R")
source("true_labels.R")


#N_list <- seq(from = 2000, to = 60000, length.out = 30)
N_list <- list(100, 500, 1000, 3000, 5000, 10000, 20000, 50000, 100000)
#N_list <- list(5000)
p <- 0.3

#a <- list(c(rep(1,40)),c(rep(-1.4,50)))
a <- list(c(rep(2,25)),c(rep(-3.4,25)), c(rep(-0.7,25)), c(rep(-5.8,25)))
#b <- c(1, 1)
b <- c(1, 1, 1, 1)
true_theta <- c(-2)
#c <- c(0.6, 5.8)
c <- c(-6, -1, 1, 4.2)
K <- length(a)
J_regular <- sum(sapply(a, length))
Q1 <- 10
Q2 <- 10
J <- Q1 + Q2 + J_regular
df <- 4
num_null <- as.integer(Q1 > 0)
num_junk <- as.integer(Q2 > 0)

rho <- rep(0, K)
phi <- rep(1, K)
gamma <- 1
kappa <- 2
num_of_clusters <- K + num_null + num_junk
rep <- 100

samples <- generate_data(50000, p, b, a, true_theta, c, h2_m = 0.5, h2_x = 0.5, h2_y = 0.5)
estimates <- estimate_theta(samples)
theta_estimates <- estimates$theta_estimates
sigma_estimates <- estimates$sigma_estimates
bx <- estimates$bx
by <- estimates$by
bxse <- estimates$bxse
byse <- estimates$byse



set.seed(1)

wrong_rows <- sample(
  1:J,
  size = ceiling(0.4 * J)
)

set.seed(1)

wrong_rows <- sample(
  1:J,
  size = ceiling(0.4 * J)
)

eta_true <- generate_priors(
  K = K,
  J = J,
  num_null = num_null,
  num_junk = num_junk,
  Q1 = Q1,
  Q2 = Q2,
  true_cluster_weight = 0.9,
  a = a,
  b = b,
  c = c,
  true_theta = true_theta,
  mode = "true"
)

eta_uniform <- generate_priors(
  K = K,
  J = J,
  num_null = num_null,
  num_junk = num_junk,
  mode = "uniform"
)

eta_wrong <- generate_priors(
  K = K,
  J = J,
  num_null = num_null,
  num_junk = num_junk,
  Q1 = Q1,
  Q2 = Q2,
  true_cluster_weight = 0.9,
  wrong_fraction = 0.4,
  wrong_rows = wrong_rows,
  a = a,
  b = b,
  c = c,
  true_theta = true_theta,
  mode = "wrong"
)

eta_true_diffuse <- diffuse_eta(
  eta = eta_true,
  new_weight = 0.6,
  diffuse_rows = 1:J
)

eta_wrong_diffuse <- diffuse_eta(
  eta = eta_wrong,
  new_weight = 0.6,
  diffuse_rows = 1:J
)
true_labels <- true_labels(a, b, c , true_theta, Q1, Q2)
