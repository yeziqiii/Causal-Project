true_labels <- function(a, b, c, true_theta, Q1 = 0, Q2 = 0
){
  centers <- true_theta + c / b
  center_order <- order(centers)
  labels <- c()
  num_null <- as.integer(Q1 > 0)
  num_junk <- as.integer(Q2 > 0)
  # null cluster
  
  if(num_null == 1){
    labels <- c(
      labels,
      rep(1,Q1)
    )
  }
  
  # regular clusters
  
  for(k in seq_along(a)){
    cluster_id <-
      which(center_order == k) + num_null
    
    labels <- c(
      labels,
      rep(
        cluster_id,
        length(a[[k]])
      )
    )
  }
  
  # junk cluster
  
  if(num_junk == 1){
    
    labels <- c(
      labels,
      rep(
        length(a) + num_null + 1,
        Q2
      )
    )
  }
  
  return(labels)
}



