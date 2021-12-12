# Functions for POPAN-lambda model

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

# Function to find proportions of superpopulation entering population in each
# survey
pent_func <- function(lambda.gaps, phi.gaps){
  pent_vec <- 
    c(1, (lambda.gaps - phi.gaps) * cumprod(c(1, lambda.gaps[1:(k - 2)])))
  pent_vec / sum(pent_vec)
}

# Vectorised pent function, as above but for vectors of parameter estimates,
# returns matrix with entry proportions as rows
pent_vec_func <- function(lmbd.gps.mat, phi.gps.mat){
  pent_mat <- 
    cbind(1, (lmbd.gps.mat - phi.gps.mat) * 
            t(apply(cbind(1, lmbd.gps.mat[, 1:(k - 2)]), 1, cumprod)))
  pent_mat / rowSums(pent_mat)
}

# Function to find final population size for each of multiple sets of parameter
# estimates
N_fin_vec_func <- function(Ns, lambda, phi){
  lmbd.gps.mat <- outer(lambda, srvy.gaps, "^")
  phi.gps.mat <- outer(phi, srvy.gaps, "^")
  rowSums(pent_vec_func(lmbd.gps.mat, phi.gps.mat) * Ns * 
            apply(phi.gps.mat, 1, prod) / 
            t(apply(cbind(1, phi.gps.mat), 1, cumprod)))
}

# Function to find superpopulation size for each of multiple sets of parameter
# estimates
Ns_vec_func <- function(N_fin, lambda, phi){
  lmbd.gps.mat <- outer(lambda, srvy.gaps, "^")
  phi.gps.mat <- outer(phi, srvy.gaps, "^")
  N_fin / lambda^stdy.len / pent_vec_func(lmbd.gps.mat, phi.gps.mat)[, 1]
}

# Function to find probability that animal alive in population and not
# previously observed at each survey
psi_func <- function(pent_vec, phi.gaps, p_vec) {
  psi_vec <- numeric(k)
  psi_vec[1] <- pent_vec[1]
  for(t in 2:k) {
    psi_vec[t] <- psi_vec[t - 1] * (1 - p_vec[t - 1]) * phi.gaps[t - 1] + 
      pent_vec[t]
  }
  psi_vec
}

# Function to find probability that animal not consequenty observed given alive
# in population at each survey
chi_func <- function(phi.gaps, p_vec) {
  chi_vec <- numeric(k)
  chi_vec[k] <- 1
  for(t in (k - 1):1) {
    chi_vec[t] <- 1 - phi.gaps[t] + phi.gaps[t] * (1 - p_vec[t + 1]) * 
      chi_vec[t + 1]
  }
  chi_vec
}

# Function to find probability that animal captured at least once in study
p_theta_func <- function(pent_vec, p_vec, chi_vec) {
  1 - sum(pent_vec * (1 - p_vec) * chi_vec)
}

# Function to calculate POPAN-lambda model negative log likelihood
popan_nll <- function(pars) {
  # print(pars)
  
  # Unpack parameters and transform according to gaps between surveys
  rho <- pars[1]
  phi <- pars[2]
  Ns <- pars[3]
  p_vec <- pars[4:(3 + k)]
  lambda <- rho + phi
  lambda.gaps <- lambda^srvy.gaps
  phi.gaps <- phi^srvy.gaps
  
  # Find entry proportions
  pent_vec <- pent_func(lambda.gaps, phi.gaps)
  
  # Find psi, chi, and ptheta
  psi_vec <- psi_func(pent_vec, phi.gaps, p_vec)
  chi_vec <- chi_func(phi.gaps, p_vec)
  p_theta <- p_theta_func(pent_vec, p_vec, chi_vec)
  
  # Find negative log-likelihood
  -sum(
    # Binomial coefficients
    lgamma(Ns + 1),
    # -lgamma(n.cap.hists + 1),
    -lgamma(Ns - n.cap.hists + 1),
    
    # Capture histories including undetected
    (Ns - n.cap.hists) * log(1 - p_theta),
    pop.sum$first.tab * log(psi_vec),
    pop.sum$caps * log(p_vec),
    pop.sum$non.caps * log(1 - p_vec),
    pop.sum$survives[-k] * log(phi.gaps),
    pop.sum$last.tab * log(chi_vec)
  )
}
