library(rstan)
library(V8)
library(DIRECT)
library(loo)
source("data_simulation.R")
source("experiment_design.R")
source("initialization.R")
source("mcmc_clust.R")
experiments(N_list, p, a, b, c, true_theta,Q1, Q2, df, etas_true, 
            rep, "results_true1.txt", "raw_true1.csv")

experiments(N_list, p, a, b, c, true_theta, Q1, Q2, df, etas_uniform, 
            rep, "results_uniform1.txt", "raw_uniform2.csv")

experiments(N_list, p, a, b, c, true_theta, Q1, Q2, df, etas_wrong, 
            rep, "results_wrong1.txt", "raw_wrong3.csv")

