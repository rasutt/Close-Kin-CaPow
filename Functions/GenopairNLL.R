# Find genopair model negative log likelihood
GenopairNLL <- function(params, k, f.year, srvy.yrs, ns.kps.lst, ns.caps) {
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
  
  # Create vector for probability of POPs and SPs over survey-pairs
  prb.POPs.mat <- prb.SPs.mat <- matrix(0, k, k)
  
  # Loop over number of surveys
  for (srvy.ind in 1:k) {
    # Find the expected number alive at the survey year
    exp.N.srvy.yr <- N.fin / lambda^(f.year - srvy.yrs[srvy.ind])
    
    # If exp.N.srvy.yr < 5 some of the kin pair probabilities can be greater
    # than one, or less than zero, producing tonnes of warnings in nlminb
    if (exp.N.srvy.yr < 5) return(Inf)
    
    # Probability of POPs within sample
    prb.POPs.mat[srvy.ind, srvy.ind] <- 
      2 / (exp.N.srvy.yr - 1) * rho * (1 + phi) / (lambda - phi^2)
  } 
  
  # Self and parent-offspring pairs between samples
  
  # Loop over all but last survey
  for (srvy.ind.1 in 1:(k - 1)) {
    # Find first survey year
    srvy.yr.1 <- srvy.yrs[srvy.ind.1]
    
    # Loop over surveys with greater indices than first
    for (srvy.ind.2 in (srvy.ind.1 + 1):k) {
      # Find second survey year
      srvy.yr.2 <- srvy.yrs[srvy.ind.2]
      
      # Find gap between surveys.  Remember not necessarily consecutive
      srvy.gap <- srvy.yr.2 - srvy.yr.1
      
      # Find the expected numbers alive in survey years
      exp.N.srvy.yr.1 <- N.fin / lambda^(f.year - srvy.yr.1)
      exp.N.srvy.yr.2 <- N.fin / lambda^(f.year - srvy.yr.2)
      
      # Probability of SPs between samples
      prb.SPs.mat[srvy.ind.1, srvy.ind.2] <- phi^srvy.gap / exp.N.srvy.yr.2
      
      # Probability of POPs between samples
      s.yrs.1.p.a = s.yrs.1 + alpha
      
      prb.POPs.mat[srvy.ind.1, srvy.ind.2] <- 
        2 * exp.ns.POPs.wtn[srvy.ind.1] * p.t.s.gs +
        2 * exp.N.s.yrs.2 * (1 - p.o.l) * p.o.l^s.yrs.2 *
        ((l.o.p^(s.yrs.1 + 1) - l.o.p^(pmin(s.yrs.1.p.a, s.yrs.2) + 1)) /
           (1 - l.o.p) +
           ifelse(
             s.yrs.1.p.a < s.yrs.2,
             (s.yrs.2 - s.yrs.1.p.a) * l.o.p^s.yrs.1.p.a,
             0
           ))
      
    }
  }
  
  # Return negative log likelihood
  nll
}