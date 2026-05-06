plot_segmented_density <- function(theta_estimates, counts,
                                   lwd = 2,
                                   xlab = "Theta estimates",
                                   ylab = "Density",
                                   main = "Cluster density for theta estimates") {
  
  cols <- rainbow(length(counts))
  
  
  
  breaks <- cumsum(counts)                 
  starts <- c(1, breaks[-length(breaks)] + 1)
  
  
  all_dens <- list()
  ymax <- 0
  
  for (i in seq_along(counts)) {
    start_idx <- starts[i]
    end_idx <- breaks[i]
    seg_theta <- theta_estimates[start_idx:end_idx]
    
    
    
    
    d <- density(seg_theta)
    all_dens[[i]] <- d  
    
    ymax <- max(ymax, max(d$y, na.rm = TRUE))
  }

  
  plot(NULL,
       xlim = range(theta_estimates, na.rm = TRUE),
       ylim = c(0, ymax * 1.1),
       xlab = xlab,
       ylab = ylab,
       main = main)
  
  for (i in seq_along(all_dens)) {
    d <- all_dens[[i]]
    lines(d$x, d$y,
          col = cols[i],
          lwd = lwd)
  }
  
  legend("topright",
         legend = paste("Cluster", seq_along(all_dens)),
         col = cols[seq_along(all_dens)],
         lwd = lwd)
}