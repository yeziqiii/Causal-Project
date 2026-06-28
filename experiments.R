library(rstan)
library(V8)
library(DIRECT)
library(loo)
#source("data_simulation.R")
source("generate_priors.R")
source("experiment_design.R")
source("initialization.R")

source("./mrclust-master/R/mr_clust_em_fixed_count.R")

N_list = c(10000, 25000, 50000, 100000)
N_list = c(50000)
p <- 0.3
a <- list(c(rep(2,25)),c(rep(-3.4,25)), c(rep(-0.7,25)), c(rep(-5.8,25)))
b <- c(1, 1, 1, 1)
true_theta <- c(-2)
c <- c(-6, -1, 1, 4.2)

# a <- list(c(rep(1,40)),c(rep(-1.4,50)))
# b <- c(1, 1)
# c <- c(0.6, 5.8)

Q1 <- 10
Q2 <- 10
df <- 4
K <- length(a)
rho <- rep(0, K)
phi <- rep(1, K)
gamma <- 1
kappa <- 2
rep <- 2

res1 <- experiment1(
  N_list = N_list,
  p = p,
  a = a,
  b = b,
  c = c,
  true_theta = true_theta,
  Q1 = Q1,
  Q2 = Q2,
  df = df,
  rep = rep,
  outputfile = "exp1_summary.csv",
  raw_outputfile = "exp1_raw.csv"
)

res2 <- experiment2(
  N_list = N_list,
  p = p,
  a = a,
  b = b,
  c = c,
  true_theta = true_theta,
  Q1 = Q1,
  Q2 = Q2,
  df = df,
  rep = rep,
  outputfile = "exp2_summary.csv",
  raw_outputfile = "exp2_raw.csv"
)

settings_exp3 <- list(
  # Centers: -6.5, 1.5
  K2 = list(
    a = list(
      rep(1.6, 25),
      rep(-2.8, 25)
    ),
    b = c(1, 1),
    c = c(-4.5, 3.5),
    true_theta = -2
  ),
  
  # Centers: -8, -3, -1, 2.2
  K4 = list(
    a = list(
      rep(2, 25),
      rep(-3.4, 25),
      rep(-0.7, 25),
      rep(-5.8, 25)
    ),
    b = c(1, 1, 1, 1),
    c = c(-6, -1, 1, 4.2),
    true_theta = -2
  ),
  
  #Centers: -7.2, -4.5, -2, 0.5, 2.8, 5
  K6 = list(
    a = list(
      rep(1.2, 25),
      rep(-2.1, 25),
      rep(3.3, 25),
      rep(-4.6, 25),
      rep(0.8, 25),
      rep(-5.2, 25)
    ),
    b = c(1, 1, 1, 1, 1, 1),
    c = c(-5.2, -2.5, 0, 2.5, 4.8, 7),
    true_theta = -2
  )
)

res3 <- experiment3(
  settings = settings_exp3,
  N = 50000,
  p = 0.3,
  Q1 = 10,
  Q2 = 10,
  df = 4,
  rep = 100,
  K_max = 6,
  rho = 0.2,
  outputfile = "exp3.csv"
)


