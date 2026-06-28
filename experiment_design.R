source("generate_data2.R")
source("rand.R")
source("determine_K.R")
source("mcmc_clust.R")
source("true_labels.R")
experiment1 <- function(
    N_list, p, a, b, c, true_theta,
    Q1 = 10,
    Q2 = 10,
    df = 4,
    rep = 100,
    outputfile = "exp1_summary.csv",
    raw_outputfile = "exp1_raw.csv",
    true_weight = 0.9,
    diffuse_weight = 0.6,
    wrong_fraction = 0.4,
    h2_m = 0.5,
    h2_x = 0.5,
    h2_y = 0.5,
    rho = 0.2
) {
  
  J_regular <- sum(sapply(a, length))
  J <- J_regular + Q1 + Q2
  
  num_null <- as.integer(Q1 > 0)
  num_junk <- as.integer(Q2 > 0)
  
  K <- length(a)
  
  true_centers <- true_theta + c / b
  true_centers <- sort(true_centers)
  
  true_lab <- true_labels(a, b, c, true_theta, Q1, Q2)
  
  prior_names <- c(
    "True-confident",
    "True-diffuse",
    "Uniform",
    "Wrong-confident",
    "Wrong-diffuse"
  )
  
  raw_results <- data.frame()
  
  for (N in N_list) {
    cat("N =", N, "\n")
    
    for (r in 1:rep) {
      
      cat("  rep =", r, "\n")
      samples <- generate_data(
        N = N,
        p = p,
        b = b,
        a = a,
        true_theta = true_theta,
        c = c,
        h2_m = h2_m,
        h2_x = h2_x,
        h2_y = h2_y, 
        rho = rho
      )
      
      estimates <- estimate_theta(samples)
      
      bx_regular <- estimates$bx
      by_regular <- estimates$by
      byse_regular <- estimates$byse
      
      # null, ensure theta samples ~ N(0, 1)
      bx_null <- rep(1, Q1)
      by_null <- rnorm(Q1, 0, 1)
      byse_null <- rep(1, Q1)
      
      # junk, ensure theta samples ~ T4
      bx_junk <- rep(1, Q2)
      by_junk <- rt(Q2, df)
      byse_junk <- rep(sqrt(2), Q2)
      
      # Combine
      bx <- c(bx_null, bx_regular, bx_junk)
      by <- c(by_null, by_regular, by_junk)
      byse <- c(byse_null, byse_regular, byse_junk)
      
      if (length(bx) != J) {
        stop("Length of bx does not match J")
      }
      
      wrong_rows <- sample(
        1:J,
        size = ceiling(wrong_fraction * J)
      )
      
      eta_true <- generate_priors(
        K = K,
        J = J,
        num_null = num_null,
        num_junk = num_junk,
        Q1 = Q1,
        Q2 = Q2,
        true_cluster_weight = true_weight,
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
        Q1 = Q1,
        Q2 = Q2,
        mode = "uniform"
      )
      
      eta_wrong <- generate_priors(
        K = K,
        J = J,
        num_null = num_null,
        num_junk = num_junk,
        Q1 = Q1,
        Q2 = Q2,
        true_cluster_weight = true_weight,
        wrong_fraction = wrong_fraction,
        wrong_rows = wrong_rows,
        a = a,
        b = b,
        c = c,
        true_theta = true_theta,
        mode = "wrong"
      )
      
      eta_true_diffuse <- diffuse_eta(
        eta = eta_true,
        new_weight = diffuse_weight,
        diffuse_rows = 1:J
      )
      
      eta_wrong_diffuse <- diffuse_eta(
        eta = eta_wrong,
        new_weight = diffuse_weight,
        diffuse_rows = 1:J
      )
      
      eta_list <- list(
        "True-confident" = eta_true,
        "True-diffuse" = eta_true_diffuse,
        "Uniform" = eta_uniform,
        "Wrong-confident" = eta_wrong,
        "Wrong-diffuse" = eta_wrong_diffuse
      )
      
      for (prior_name in prior_names) {
        
        cat("    prior =", prior_name, "\n")
        
        fit_res <- mcmc_clust(
          bx = bx,
          by = by,
          byse = byse,
          eta = eta_list[[prior_name]],
          num_null = num_null,
          num_junk = num_junk,
          df = df
        )
        
        ri_value <- rand_mcmc(
          true_labels = true_lab,
          pred_labels = fit_res$cluster_assignment,
          K = K,
          num_null = num_null,
          num_junk = num_junk
        )
        
        
        row <- data.frame(
          N = N,
          rep = r,
          prior = prior_name,
          RI = ri_value,
          WAIC = fit_res$waic
        )
        
        for (k in 1:K) {
          row[[paste0("center", k)]] <- round(
            fit_res$center_estimates[k],
            3
          )
        }
        # Store data in raw results by combining results in each loop
        raw_results <- rbind(
          raw_results,
          row
        )
      }
    }
  }
  
  # Store raw data
  write.csv(
    raw_results,
    file = raw_outputfile,
    row.names = FALSE
  )
  
  summary_results <- data.frame()
  
  for (N in N_list) {
    for (prior_name in prior_names) {
      # Consider the sub dataset:
      # Results of all repetitions for fixed N and prior
      data <- raw_results[
        raw_results$N == N & raw_results$prior == prior_name,
      ]
      
      center_mat <- as.matrix(
        data[, paste0("center", 1:K)]
      )
      
      mean_centers <- colMeans(center_mat)
      
      # Find variance for each substantive center estimates
      var_centers <- apply(center_mat, 2, var)
      
      # Calculate squared error for each substantive center estimates
      sq_error_mat <- sweep(
        center_mat,
        MARGIN = 2,
        STATS = true_centers,
        FUN = "-"
      )^2
      
      mse_centers <- colMeans(sq_error_mat)
      
      summary_row <- data.frame(
        N = N,
        Prior_specification = prior_name
      )
      
      for (k in 1:K) {
        summary_row[[paste0("estimated_center", k)]] <- signif(
          mean_centers[k],
          4
        )
      }
      
      for (k in 1:K) {
        summary_row[[paste0("variance_center", k)]] <- signif(
          var_centers[k],
          4
        )
      }
      
      for (k in 1:K) {
        summary_row[[paste0("mse_center", k)]] <- signif(
          mse_centers[k],
          4
        )
      }
      
      summary_row$RI <- signif(
        mean(data$RI, na.rm = TRUE),
        4
      )
      
      summary_results <- rbind(
        summary_results,
        summary_row
      )
    }
  }
  
  write.csv(
    summary_results,
    file = outputfile,
    row.names = FALSE
  )
  
  return(list(
    summary = summary_results,
    raw = raw_results
  ))
}

extract_mrclust_centers <- function(res_em, K) {
  
  best <- res_em$results$best
  
  center_df <- aggregate(
    cluster_mean ~ cluster_class,
    data = best,
    FUN = mean
  )
  
  centers <- center_df$cluster_mean[
    match(1:K, center_df$cluster_class)
  ]
  
  centers <- sort(centers)
  
  return(centers)
}

experiment2 <- function(
    N_list, p, a, b, c, true_theta,
    Q1 = 10,
    Q2 = 10,
    df = 4,
    rep = 100,
    outputfile = "exp2_summary.csv",
    raw_outputfile = "exp2_raw.csv",
    true_weight = 0.9,
    diffuse_weight = 0.6,
    wrong_fraction = 0.4,
    h2_m = 0.5,
    h2_x = 0.5,
    h2_y = 0.5,
    rho = 0.2
) {
  J_regular <- sum(sapply(a, length))
  J <- J_regular + Q1 + Q2
  
  num_null <- as.integer(Q1 > 0)
  num_junk <- as.integer(Q2 > 0)
  
  K <- length(a)
  
  true_centers <- true_theta + c / b
  true_centers <- sort(true_centers)
  
  true_lab <- true_labels(a, b, c, true_theta, Q1, Q2)
  
  method_names <- c(
    "Proposed",
    "MR-Clust"
  )
  
  raw_results <- data.frame()
  
  for (N in N_list) {
    cat("N =", N, "\n")
    for (r in 1:rep) {
      cat("  rep =", r, "\n")
      samples <- generate_data(
        N = N,
        p = p,
        b = b,
        a = a,
        true_theta = true_theta,
        c = c,
        h2_m = h2_m,
        h2_x = h2_x,
        h2_y = h2_y,
        rho = rho
      )
      
      estimates <- estimate_theta(samples)
      
      # Regular variants
      bx_regular <- estimates$bx
      by_regular <- estimates$by
      bxse_regular <- estimates$bxse
      byse_regular <- estimates$byse
      
      if (length(bx_regular) != J_regular) {
        stop("Length of bx_regular does not match J_regular.")
      }
      
      # Null variants
      bx_null <- rep(1, Q1)
      by_null <- rnorm(Q1, mean = 0, sd = 1)
      bxse_null <- rep(0, Q1)
      byse_null <- rep(1, Q1)
      
      # Junk variants
      bx_junk <- rep(1, Q2)
      by_junk <- rt(Q2, df = df)
      bxse_junk <- rep(0, Q2)
      byse_junk <- rep(sqrt(2), Q2)
      
      # Full data: null + regular + junk
      bx <- c(
        bx_null,
        bx_regular,
        bx_junk
      )
      
      by <- c(
        by_null,
        by_regular,
        by_junk
      )
      
      bxse <- c(
        bxse_null,
        bxse_regular,
        bxse_junk
      )
      
      byse <- c(
        byse_null,
        byse_regular,
        byse_junk
      )
      
      if (length(bx) != J) {
        stop("Length of bx does not match J.")
      }
      
      
      theta_full <- by / bx
      theta_se_full <- byse / abs(bx)
      
      obs_names <- paste0("element", 1:J)
      
      # Proposed method use true-diffuse
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
      
      eta_right_diffuse <- diffuse_eta(
        eta = eta_true,
        new_weight = diffuse_weight,
        diffuse_rows = 1:J
      )
      
      fit_mcmc <- mcmc_clust(
        bx = bx,
        by = by,
        byse = byse,
        eta = eta_right_diffuse,
        num_null = num_null,
        num_junk = num_junk,
        df = df
      )
      
      ri_mcmc <- rand_mcmc(
        true_labels = true_lab,
        pred_labels = fit_mcmc$cluster_assignment,
        K = K,
        num_null = num_null,
        num_junk = num_junk
      )
      
      row_mcmc <- data.frame(
        N = N,
        rep = r,
        method = "Proposed",
        RI = ri_mcmc
      )
      
      for (k in 1:K) {
        row_mcmc[[paste0("center", k)]] <- round(fit_mcmc$center_estimates[k], 3)
      }
      
      raw_results <- rbind(
        raw_results,
        row_mcmc
      )
      
      # MR-Clust with fixed K substantive clusters
      res_em <- mr_clust_em_fixed_count(
        theta = theta_full,
        theta_se = theta_se_full,
        bx = bx,
        by = by,
        bxse = bxse,
        byse = byse,
        obs_names = obs_names,
        num_substantive_clusters = K,
        plot_results = NULL
      )
      
      mr_centers <- extract_mrclust_centers(
        res_em = res_em,
        K = K
      )
      
      mr_labels <- res_em$results$best$cluster_class
      
      ri_mr <- rand_mrclust(
        true_labels = true_lab,
        pred_labels = mr_labels,
        K = K,
        num_null = num_null,
        num_junk = num_junk
      )
      
      row_mr <- data.frame(
        N = N,
        rep = r,
        method = "MR-Clust",
        RI = ri_mr
      )
      
      for (k in 1:K) {
        row_mr[[paste0("center", k)]] <- mr_centers[k]
      }
      
      raw_results <- rbind(
        raw_results,
        row_mr
      )
    }
  }
  
  write.csv(
    raw_results,
    file = raw_outputfile,
    row.names = FALSE
  )
  
  summary_results <- data.frame()
  
  for (N in N_list) {
    for (method_name in method_names) {
      data <- raw_results[
        raw_results$N == N &
          raw_results$method == method_name,
      ]
      
      center_mat <- as.matrix(
        data[, paste0("center", 1:K)]
      )
      
      mean_centers <- colMeans(
        center_mat,
      )
      
      summary_row <- data.frame(
        N = N,
        Method = method_name,
        RI = signif(
          mean(data$RI, na.rm = TRUE),
          4
        )
      )
      
      for (k in 1:K) {
        summary_row[[paste0("estimated_center", k)]] <- signif(
          mean_centers[k],
          4
        )
      }
      
      summary_results <- rbind(
        summary_results,
        summary_row
      )
    }
  }
  
  write.csv(
    summary_results,
    file = outputfile,
    row.names = FALSE
  )
  
  return(list(
    summary = summary_results,
    raw = raw_results
  ))
}

experiment3 <- function(
    settings,
    N,
    p = 0.3,
    Q1 = 10,
    Q2 = 10,
    df = 4,
    rep = 100,
    K_max = 6,
    outputfile = "exp3_K_selection_counts.csv",
    h2_m = 0.5,
    h2_x = 0.5,
    h2_y = 0.5,
    rho = 0.2
) {
  
  num_null <- as.integer(Q1 > 0)
  num_junk <- as.integer(Q2 > 0)
  
  summary_results <- data.frame()
  
  for (setting_name in names(settings)) {
    
    setting <- settings[[setting_name]]
    
    a <- setting$a
    b <- setting$b
    c <- setting$c
    true_theta <- setting$true_theta
    
    K0 <- length(a)
    J_regular <- sum(sapply(a, length))
    J <- J_regular + Q1 + Q2
    
    cat("True K0 =", K0, "\n")
    
    estimated_K_vec <- numeric(rep)
    
    for (r in 1:rep) {
      
      cat("  rep =", r, "\n")
      
      samples <- generate_data(
        N = N,
        p = p,
        b = b,
        a = a,
        true_theta = true_theta,
        c = c,
        h2_m = h2_m,
        h2_x = h2_x,
        h2_y = h2_y,
        rho = rho
      )
      
      estimates <- estimate_theta(samples)
      
      bx_regular <- estimates$bx
      by_regular <- estimates$by
      byse_regular <- estimates$byse
      
      if (length(bx_regular) != J_regular) {
        stop("Length of bx_regular does not match J_regular.")
      }
      
      # Null
      bx_null <- rep(1, Q1)
      by_null <- rnorm(Q1, mean = 0, sd = 1)
      byse_null <- rep(1, Q1)
      
      # Junk
      bx_junk <- rep(1, Q2)
      by_junk <- rt(Q2, df = df)
      byse_junk <- rep(sqrt(2), Q2)
      
      # Full data
      bx <- c(
        bx_null,
        bx_regular,
        bx_junk
      )
      
      by <- c(
        by_null,
        by_regular,
        by_junk
      )
      
      byse <- c(
        byse_null,
        byse_regular,
        byse_junk
      )
      
      if (length(bx) != J) {
        stop("Length of bx does not match J.")
      }
      
      fit_k <- determine_K(
        bx = bx,
        by = by,
        byse = byse,
        num_null = num_null,
        num_junk = num_junk,
        K_max = K_max,
        df = df
      )
      
      estimated_K_vec[r] <- fit_k$best_K
    }
    
    estimated_counts <- table(
      factor(
        estimated_K_vec,
        levels = 2:K_max
      )
    )
    
    summary_row <- data.frame(
      setting = setting_name,
      K0 = K0
    )
    
    for (k in 2:K_max) {
      summary_row[[paste0("estimated_K", k)]] <- as.numeric(
        estimated_counts[as.character(k)]
      )
    }
    
    summary_row$accuracy <- mean(estimated_K_vec == K0)
    
    summary_results <- rbind(
      summary_results,
      summary_row
    )
  }
  
  write.csv(
    summary_results,
    file = outputfile,
    row.names = FALSE
  )
  
  return(summary_results)
}