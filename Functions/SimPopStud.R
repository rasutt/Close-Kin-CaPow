# Function to simulate a population of animals over time and a mark-recapture
# study of it.

# Returns: Dataframe with IDs, parents (NA for initial population), and
# capture histories for captured animals, with parameters, implied birthrate
# beta, and population trajectory N.t.vec attached.
SimPopStud <- function(phi, lambda, N.init, hist.len, srvy.yrs, k, f.year, p) {
  # Record start-time
  s.time <- proc.time()
  
  # Find implied birth rate for mature females surviving to birth year
  beta <- 2 * (lambda / phi - 1) * (lambda / phi)^alpha
  if (beta < 0 | beta > 1)
    cat("Implied birthrate for mature females:", round(beta, 3), "\n")
  if (beta < 0) stop("Negative birth rates impossible")
  if (beta > 1) stop("Maximum one calf at a time")
  
  # Set birthyears, parents, sexes, times when animals last had a calf, and life
  # statuses for first animals
  init.ages <- rgeom(n = N.init, prob = 1 - phi / lambda)
  b.year <- sort(1 - init.ages)
  mum <- rep(NA, N.init)
  dad <- rep(NA, N.init)
  female <- rbinom(N.init, 1, 0.5)
  t.lst.clf <- rep(NA, N.init)
  alive <- rep(T, N.init)
  
  # Create vector for population size and enter first value
  N.t.vec <- numeric(hist.len)
  N.t.vec[1] <- N.init
  
  # Create lists for life and calving statuses of animals in survey years
  alv.srvy <- vector("list", k)
  clvng.srvy <- vector("list", k)
  
  # Set survey counter to zero
  srvy.cnt <- 0
  
  # If survey year increment counter and add life statuses of animals to list
  if ((f.year - hist.len + 1) %in% srvy.yrs) {
    srvy.cnt <- srvy.cnt + 1
    alv.srvy[[srvy.cnt]] <- alive
  }
  
  # Display progress
  # cat('Initialized. Year =', f.year - hist.len + 1)
  
  # Loop over remaining years generating data based on the previous year
  for (t in 2:hist.len) {
    # Display progress
    # if (t %% 20 == 1) cat('', f.year - hist.len + t)
    
    # Find animals alive and mature last year
    mature <- alive & t - 1 - b.year >= alpha
    
    # Find survivors to current year
    alive[alive] <- as.logical(
      rbinom(N.t.vec[t - 1], 1, phi * (1 - pmt.emgn * !female[alive]))
    )
    
    # Find possible parents (mothers must survive to give birth)
    mums.poss <- which(alive & mature & female)
    dads.poss <- which(mature & !female)
    
    # If there was at least one possible father find number of calves
    if (length(dads.poss) > 0) 
      # Either stochastic or deterministic numbers of births
      if (stch.bths) n.calves <- rbinom(1, length(mums.poss), beta)
      else round(N.t.vec[t - 1] * lambda - sum(alive))
    else n.calves <- 0

    # If there are calves
    if (n.calves > 0) {
      # Order possible mothers by time since last calving
      mums.poss.ord <- mums.poss[order(t.lst.clf[mums.poss], na.last = F)]
      
      # Find possible mothers tied for longest time since last calving
      lgst <- mums.poss.ord == mums.poss.ord[n.calves]
      
      # If more than one then order randomly
      if (sum(lgst) > 1) mums.poss.ord[lgst] <- sample(mums.poss.ord[lgst])
      
      # Find mothers
      mums.new <- mums.poss.ord[1:n.calves]
      mum <- c(mum, mums.new)
      t.lst.clf[mums.new] <- t
      
      # If more than one possible father then choose randomly for each calf
      if (length(dads.poss) > 1) 
        dad <- c(dad, sample(dads.poss, n.calves, replace = T))
      else dad <- c(dad, rep(dads.poss, n.calves))
      
      # Add data for calves
      b.year <- c(b.year, rep(t, n.calves))
      female <- c(female, rbinom(n.calves, 1, 0.5))
      t.lst.clf <- c(t.lst.clf, rep(NA, n.calves))
      alive <- c(alive, rep(T, n.calves))
      
    } # End of if there are calves
    
    # Record population size
    N.t.vec[t] <- sum(alive)

    # If survey year add life and calving statuses of animals to lists
    if ((f.year - hist.len + t) %in% srvy.yrs) {
      srvy.cnt <- srvy.cnt + 1
      alv.srvy[[srvy.cnt]] <- alive
      clvng.srvy[[srvy.cnt]] <- t.lst.clf == t
    }
    
  } # End of loop over times
  
  # Display progress
  # cat('', f.year - hist.len + t, '\n')
  
  # Adjust birth years with respect to final year
  b.year <- b.year + f.year - hist.len
  
  # Create matrices and enter life and calving statuses of animals in survey
  # years
  cap.hists <- clvng.hists <- matrix(F, length(alive), k)
  for (srvy.ind in 1:k) {
    n.alv.srvy <- length(alv.srvy[[srvy.ind]])
    cap.hists[1:n.alv.srvy, srvy.ind] <- alv.srvy[[srvy.ind]]
    clvng.hists[1:n.alv.srvy, srvy.ind] <- clvng.srvy[[srvy.ind]]
  }
  clvng.hists[is.na(clvng.hists)] <- F
  mode(clvng.hists) <- "integer"
  
  # Find superpopulation size of study
  Ns <- sum(rowSums(cap.hists) > 0)
  
  # Find numbers calving in survey years
  ns.clvng <- colSums(clvng.hists)
  
  # Change life statuses to capture histories
  # The capture probability depends on male temporary emigration, and calving
  cap.hists[cap.hists] <- 
    rbinom(sum(cap.hists), 1, p * (1 - tmp.emgn * !rep(female, k)[cap.hists]) + 
             clvng.hists[cap.hists] * clvng.p)
  
  # Label capture and calving history columns by year
  colnames(cap.hists) <- paste0("C", srvy.yrs)
  colnames(clvng.hists) <- paste0("Cvg", srvy.yrs)
  
  # Make animal IDs
  ID <- seq_along(alive)
  
  # Combine results for captured animals
  pop.hist <- 
    data.frame(ID, mum, dad, cap.hists, clvng.hists)[rowSums(cap.hists) > 0, ]
  
  # Attach parameters and implied birthrate
  attributes(pop.hist)$beta <- beta
  attributes(pop.hist)$N.t.vec <- N.t.vec
  attributes(pop.hist)$ns.caps <- colSums(cap.hists)
  attributes(pop.hist)$Ns <- Ns
  attributes(pop.hist)$ns.clvng <- ns.clvng
  
  # Display runtime
  # cat("Final population size:", tail(N.t.vec, 1), "\n")
  # cat("Sim took", (proc.time() - s.time)["user.self"], "seconds \n")
  
  # Return population history data 
  pop.hist
}