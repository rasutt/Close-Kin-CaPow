# Find close kin model negative log likelihood
CloseKinNLL <- function(params, k, f.year, srvy.yrs, ns.kps.lst, ns.caps) {
  # print(params)
  
  # Unpack parameters
  rho <- params[1]
  phi <- params[2]
  N.fin <- params[3]
  lambda <- rho + phi
  
  # Rho must be strictly greater than zero for the closed form series below
  # to hold, otherwise nlminb gives warnings
  if (rho == 0) return(Inf)
  
  # Set negative log likelihood to zero
  nll <- 0
  
  # Parent-offspring pairs within samples
  
  # Create vector for probability of POPs within samples
  prb.POPs.wtn.vec <- numeric(k)
  
  # Loop over number of surveys
  for (srvy.ind in 1:k) {
    # Find the expected number alive at the survey year
    exp.N.srvy.yr <- N.fin / lambda^(f.year - srvy.yrs[srvy.ind])
    
    # If exp.N.srvy.yr < 5 some of the kin pair probabilities can be greater
    # than one, or less than zero, producing tonnes of warnings in nlminb
    if (exp.N.srvy.yr < 5) return(Inf)
    
    # Probability of POPs within sample
    prb.POPs.wtn <- 4 / (exp.N.srvy.yr - 1) * phi * rho / (lambda - phi^2)
    
    # Store probability of POPs between samples
    prb.POPs.wtn.vec[srvy.ind] <- prb.POPs.wtn
    
    # Probability of HSPs within sample
    prb.HSPs.wtn <- 4 * rho / (exp.N.srvy.yr - 1) * (1 - phi / lambda) * 
      (lambda / phi)^alpha * 
      (lambda * (phi + lambda) / (lambda - phi^2)^2 -
         4 * rho * (lambda / phi)^alpha * phi^2 * lambda / (lambda - phi^3)^2)
    
    # Find number of POPs, HSPs, and non-POPs-non-HSPs within sample, assuming
    # no mixed kinships
    POPs.wtn <- ns.kps.lst$ns.POPs.wtn[srvy.ind]
    HSPs.wtn <- ns.kps.lst$ns.HSPs.wtn[srvy.ind]
    non.POPs.HSPs.wtn <- choose(ns.caps[srvy.ind], 2) - POPs.wtn - HSPs.wtn
    
    # Add negative log likelihood from number of POPs within sample
    nll <- nll - POPs.wtn * log(prb.POPs.wtn) - HSPs.wtn * log(prb.HSPs.wtn) -
      non.POPs.HSPs.wtn * log(1 - prb.POPs.wtn - prb.HSPs.wtn)
  } 
  
  # Self and parent-offspring pairs between samples
  
  # Set pair counter to zero
  pr.cnt <- 0
  
  # Loop over all but last survey
  for (srvy.ind.1 in 1:(k - 1)) {
    # Find first survey year
    srvy.yr.1 <- srvy.yrs[srvy.ind.1]
    
    # Loop over surveys with greater indices than first
    for (srvy.ind.2 in (srvy.ind.1 + 1):k) {
      # Increment pair counter
      pr.cnt <- pr.cnt + 1
      
      # Find second survey year
      srvy.yr.2 <- srvy.yrs[srvy.ind.2]
      
      # Find gap between surveys.  Remember not necessarily consecutive
      srvy.gap <- srvy.yr.2 - srvy.yr.1
      
      # Find the expected numbers alive in survey years
      exp.N.srvy.yr.1 <- N.fin / lambda^(f.year - srvy.yr.1)
      exp.N.srvy.yr.2 <- N.fin / lambda^(f.year - srvy.yr.2)
      
      # Probability of SPs between samples
      prb.SPs.btn <- phi^srvy.gap / exp.N.srvy.yr.2
      
      # Probability of POPs where one born between samples
      prb.POPs.brn.btn <- 2 * rho / lambda / exp.N.srvy.yr.1 *
        sum((phi / lambda)^pmax((srvy.gap - 1):0, srvy.gap - alpha - 1))
      
      # Add modified probability from POPS that existed in the first survey.
      # Strangely simplified with self-pair probability lol
      prb.POPs.btn <- prb.POPs.brn.btn + prb.POPs.wtn.vec[srvy.ind.1] * 
        prb.SPs.btn * (exp.N.srvy.yr.1 - 1)

      # Find numbers of POPs, SPs, and other pairs between samples
      SPs.btn <- ns.kps.lst$ns.SPs.btn[pr.cnt]
      POPs.btn <- ns.kps.lst$ns.POPs.btn[pr.cnt]
      non.POPs.SPs.btn <- ns.caps[srvy.ind.1] * ns.caps[srvy.ind.2] - 
        SPs.btn - POPs.btn
      
      # Add negative log likelihood from numbers of SPs and POPs observed
      nll <- nll - SPs.btn * log(prb.SPs.btn) - POPs.btn * log(prb.POPs.btn) -
        non.POPs.SPs.btn * log(1 - prb.POPs.btn - prb.SPs.btn)
    }
  }
  
  # Return negative log likelihood
  nll
}