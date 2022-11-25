# Function to simulate a population of animals over time and a mark-recapture
# study of it.

# Returns: Data frame with IDs, parents (NA for initial population), and
# capture histories for captured animals, with parameters, implied birthrate
# beta, and population trajectory N.t.vec attached.
SimPopStud <- function(
  phi, lambda, N.init, hist.len, srvy.yrs, k, f.year, p, L, clvng.p, tmp.emgn,
  alpha, clvng.ints
) {
  # Record start-time
  s.time <- proc.time()
  
  # Find implied birth rate for mature females surviving to birth year
  # beta <- 2 * (lambda / phi - 1) * (lambda / phi)^alpha
  beta <- 2 * (1 - phi / lambda) * (lambda / phi)^alpha
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
  
  # Set initial genotypes
  gts = array(as.integer(runif(2 * L * N.init) > 0.5), c(2, L, N.init))
  
  # Create vectors for population size and observed survival rates and enter
  # first value
  N.t.vec <- numeric(hist.len)
  N.t.vec[1] <- N.init
  phi.obs = numeric(hist.len - 1)

  # Create lists for life and calving statuses of animals in survey years
  alv.srvy <- clvng.srvy <- vector("list", k)
  
  # I was keeping track of life statuses in all years for testing derivations,
  # just commenting out in case I want to again later
  # alv.all = vector("list", hist.len)
  # alv.all[[1]] = alive

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
    
    # Find animals mature this year
    mature <- t - b.year >= alpha
    
    # Find possible dads (alive last year and would be mature this year)
    dads.poss <- which(alive & mature & !female)

    # Find survivors to current year
    alive[alive] <- as.logical(
      # rbinom(N.t.vec[t - 1], 1, phi * (1 - pmt.emgn * !female[alive]))
      rbinom(N.t.vec[t - 1], 1, phi)
    )
    phi.obs[t - 1] = sum(alive) / N.t.vec[t - 1]

    # Find possible mums (alive and mature this year)
    mums.poss <- which(alive & mature & female)
    
    # If there was at least one possible father find number of calves
    if (length(dads.poss) > 0) 
      # Either stochastic or deterministic numbers of births
      if (stch.bths) n.calves <- rbinom(1, length(mums.poss), beta)
      else round(N.t.vec[t - 1] * lambda - sum(alive))
    else n.calves <- 0

    # If there are calves
    if (n.calves > 0) {
      # If calving intervals requested
      if (clvng.ints) {
        # Order possible mothers by time since last calving
        mums.poss.ord <- mums.poss[order(t.lst.clf[mums.poss], na.last = F)]
        
        # Find possible mothers tied for longest time since last calving
        lgst <- mums.poss.ord == mums.poss.ord[n.calves]
        
        # If more than one then order randomly
        if (sum(lgst) > 1) mums.poss.ord[lgst] <- sample(mums.poss.ord[lgst])
      } else {
        # Randomly select mothers
        if (length(mums.poss) > 1) mums.poss.ord <- sample(mums.poss)
        else mums.poss.ord <- mums.poss
      }
      
      # Find mothers
      mums.new <- mums.poss.ord[1:n.calves]
      mum <- c(mum, mums.new)
      t.lst.clf[mums.new] <- t
      
      # If more than one possible father then choose randomly for each calf
      if (length(dads.poss) > 1) 
        dads.new <- sample(dads.poss, n.calves, replace = T)
      else dads.new <- rep(dads.poss, n.calves)
      dad <- c(dad, dads.new)
      
      # Add data for calves
      b.year <- c(b.year, rep(t, n.calves))
      female <- c(female, rbinom(n.calves, 1, 0.5))
      t.lst.clf <- c(t.lst.clf, rep(NA, n.calves))
      alive <- c(alive, rep(T, n.calves))
      
      # Add genotypes for calves. Arithmetic seems to be a tiny bit faster. Keep
      # both versions to help with trying in pytorch. If change should also
      # change type of initial genotypes.
      gt.slct = array(runif(2 * L * n.calves) > 0.5, c(2, L, n.calves))
      # gt.clvs = array(F, c(2, L, n.calves))
      # gt.clvs[aperm(array(
      #   c(
      #     (gt.slct[1, , ] & gt[1, , mums.new]) |
      #       (!gt.slct[1, , ] & gt[2, , mums.new]),
      #     (gt.slct[2, , ] & gt[1, , dads.new]) |
      #       (!gt.slct[2, , ] & gt[2, , dads.new])
      #   ), c(L, n.calves, 2)
      # ), c(3, 1, 2))] = T
      gt.clvs = aperm(array(
        c(
          gt.slct[1, , ] * gts[1, , mums.new] +
            (!gt.slct[1, , ]) * gts[2, , mums.new],
          gt.slct[2, , ] * gts[1, , dads.new] +
            (!gt.slct[2, , ]) * gts[2, , dads.new]
        ), c(L, n.calves, 2)
      ), c(3, 1, 2))
      
      gts = array(c(gts, gt.clvs), c(2, L, length(alive)))
      
    } # End of if there are calves
    
    # Record population size
    N.t.vec[t] <- sum(alive)
    
    # Record life statuses
    # alv.all[[t]] = alive
    
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
  alv.s.yrs <- cap.hists <- clvng.hists <- matrix(F, length(alive), k)
  for (srvy.ind in 1:k) {
    n.alv.srvy <- length(alv.srvy[[srvy.ind]])
    alv.s.yrs[1:n.alv.srvy, srvy.ind] <- alv.srvy[[srvy.ind]]
    clvng.hists[1:n.alv.srvy, srvy.ind] <- clvng.srvy[[srvy.ind]]
  }
  clvng.hists[is.na(clvng.hists)] <- F
  mode(clvng.hists) <- "integer"

  # Make matrix for life histories 
  # alv.mat = matrix(F, length(alive), hist.len)
  # for (t in 1:hist.len) {
  #   alv.mat[1:length(alv.all[[t]]), t] = alv.all[[t]]
  # }

  # Find super-population size of study
  Ns <- sum(rowSums(alv.s.yrs) > 0)
  
  # Change life statuses to capture histories
  # The capture probability depends on male temporary emigration, and calving
  cap.hists[alv.s.yrs] <- 
    rbinom(sum(alv.s.yrs), 1, p * (1 - tmp.emgn * !rep(female, k)[alv.s.yrs]) + 
             clvng.hists[alv.s.yrs] * clvng.p)
  
  # Find numbers calving in survey years
  ns.clvng <- colSums(clvng.hists)
  
  # Label capture and calving history columns by year
  colnames(cap.hists) <- paste0("C", srvy.yrs)
  colnames(clvng.hists) <- paste0("Cvg", srvy.yrs)
  
  # Make animal IDs
  ID <- seq_along(alive)
  
  # Combine results for captured animals
  pop.hist <- 
    data.frame(ID, mum, dad, cap.hists, clvng.hists)[rowSums(cap.hists) > 0, ]
  
  # rownames(alv.mat) <- ID
  
  # Attach parameters and implied birthrate
  attributes(pop.hist)$avg.phi.obs <- mean(phi.obs)
  attributes(pop.hist)$beta <- beta
  attributes(pop.hist)$N.t.vec <- N.t.vec
  attributes(pop.hist)$ns.caps <- colSums(cap.hists)
  attributes(pop.hist)$Ns <- Ns
  attributes(pop.hist)$ns.clvng <- ns.clvng
  # attributes(pop.hist)$alv.mat <- alv.mat
  attributes(pop.hist)$alive <- alive
  attributes(pop.hist)$alv.s.yrs <- alv.s.yrs
  attributes(pop.hist)$f.age <- f.year - b.year
  attributes(pop.hist)$mum <- mum
  attributes(pop.hist)$dad <- dad
  attributes(pop.hist)$ID <- ID
  attributes(pop.hist)$unq.smp.gts <- gts[, , rowSums(cap.hists) > 0]
  
  # Display runtime
  # cat("Final population size:", tail(N.t.vec, 1), "\n")
  # cat("Sim took", (proc.time() - s.time)["user.self"], "seconds \n")
  
  # Return population history data 
  pop.hist
}