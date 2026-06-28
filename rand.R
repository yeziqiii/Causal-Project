rand_index <- function(
    true_labels,
    pred_labels
) {
  J <- length(true_labels)
  
  total_pairs <- choose(J, 2)
  agree <- 0
  
  for (i in 1:(J - 1)) {
    for (j in (i + 1):J) {
      same_true <- true_labels[i] == true_labels[j]
      same_pred <- pred_labels[i] == pred_labels[j]
      if (same_true == same_pred) {
        agree <- agree + 1
      }
    }
  }
  
  rand_index <- agree / total_pairs
  
  return(rand_index)
}


rand_mcmc <- function(
    true_labels,
    pred_labels,
    K,
    num_null = 1,
    num_junk = 1
) {
  
  if (length(true_labels) != length(pred_labels)) {
    stop("true_labels and pred_labels must have the same length.")
  }
  
  true_new <- true_labels
  pred_new <- pred_labels
  
  # Only combine null and junk when both are present
  if (num_null == 1 && num_junk == 1) {
    
    null_label <- 1
    junk_label <- K + num_null + num_junk
    
    true_new[true_new == junk_label] <- null_label
    pred_new[pred_new == junk_label] <- null_label
  }
  
  return(rand_index(true_new, pred_new))
}



rand_mrclust <- function(
    true_labels,
    pred_labels,
    K,
    num_null = 1,
    num_junk = 1
) {
  
  if (length(true_labels) != length(pred_labels)) {
    stop("true_labels and pred_labels must have the same length.")
  }
  
  true_new <- true_labels
  pred_new <- pred_labels
  
  # Combine null and junk when both are present
  if (num_null == 1 && num_junk == 1) {
    
    # True-label convention:
    # null: 1
    # substantive clusters: 2,...,K+1
    # junk: K + 2
    true_null_label <- 1
    true_junk_label <- K + 2
    
    true_new[true_new == true_junk_label] <- true_null_label
    
    
    # MR-Clust convention:
    # null: K+1 
    # substantive clusters: 1,...,K
    # junk: K+2
    pred_null_label <- K + 1
    pred_junk_label <- K + 2
    
    pred_new[pred_new == pred_junk_label] <- pred_null_label
  }
  
  return(rand_index(true_new, pred_new))
}