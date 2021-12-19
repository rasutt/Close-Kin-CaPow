# Summarise data for POPAN model
FindPopSum <- function() {
  # Get capture histories as matrix
  cap.hists <- as.matrix(pop.cap.hist[, 4:(3 + k)])
  
  # Create objects to hold summaries for each capture history
  non.caps.mat <- survive.mat <- matrix(0, nrow = n.cap.hists, ncol = k)
  first.obs <- last.obs <- numeric(n.cap.hists)
  
  # Loop over capture histories
  for(cp.hst.ind in 1:n.cap.hists){
    # Find first and last captures
    which.1 <- which(cap.hists[cp.hst.ind, ] == 1)
    first <- min(which.1)
    last <- max(which.1)
    first.obs[cp.hst.ind] <- first
    last.obs[cp.hst.ind] <- last
    
    # Find non-captures and survival between first and last captures
    if(last > first){
      non.caps.mat[cp.hst.ind, first:last] <- 
        1 - cap.hists[cp.hst.ind, first:last]
      survive.mat[cp.hst.ind, first:(last - 1)] <- 1
    }
  }
  
  # Find and return sufficient statistics
  data.frame(
    first.tab = as.vector(table(factor(first.obs, levels=1:k))), 
    caps = colSums(cap.hists), 
    non.caps = colSums(non.caps.mat), 
    survives = c(colSums(survive.mat)[1:(k-1)], NA),
    last.tab = as.vector(table(factor(last.obs, levels=1:k)))
  )
}