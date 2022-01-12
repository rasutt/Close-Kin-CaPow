# Function to find expected numbers of kinpairs 
FindExpNsKPs <- function(
  k, n.srvy.prs, exp.N.fin, lambda, f.year, srvy.yrs, phi, rho, ns.caps, alpha
) {
  # Expected numbers of pairs within each survey and between each pair of
  # surveys
  exp.N.srvy.yrs = exp.N.fin / lambda^(f.year - srvy.yrs)
  exp.ns.APs.wtn = choose(exp.N.srvy.yrs, 2)
  exp.ns.APs.btn = combn(exp.N.srvy.yrs, 2, function(x) x[1] * x[2])
  
  # Create vectors for POPs within each sample, and POPs and SPs between each
  # pair of samples
  exp.ns.POPs.wtn <- exp.ns.HSPs.wtn <- exp.ns.SMPs.wtn <- numeric(k)
  exp.ns.POPs.btn <- exp.ns.SPs.btn <- numeric(n.srvy.prs)
  
  # Parent-offspring pairs within sample
  
  # Loop over number of surveys
  for (srvy.ind in 1:k) {
    # Find the expected number alive at the survey year
    exp.N.srvy.yr <- exp.N.fin / lambda^(f.year - srvy.yrs[srvy.ind])
    
    # Probability of POPs within sample
    prb.POPs.wtn <- 4 / (exp.N.srvy.yr - 1) * phi * rho / (lambda - phi^2)
    
    # Find expected number of POPs within sample
    exp.ns.POPs.wtn[srvy.ind] <- prb.POPs.wtn * choose(ns.caps[srvy.ind], 2)
    
    # Probability of HSPs within sample
    prb.HSPs.wtn <- 4 * rho / (exp.N.srvy.yr - 1) * (1 - phi / lambda) * 
      (lambda / phi)^alpha * 
      (lambda * (phi + lambda) / (lambda - phi^2)^2 -
         4 * rho * (lambda / phi)^alpha * phi^2 * lambda / (lambda - phi^3)^2)
    
    # Find expected number of HSPs within sample
    exp.ns.HSPs.wtn[srvy.ind] <- prb.HSPs.wtn * choose(ns.caps[srvy.ind], 2)
  
    # Same-mother pairs
    exp.ns.SMPs.wtn[srvy.ind] <- 4 / exp.N.srvy.yr * lambda * 
      (1 - phi / lambda)^2 * 
      (lambda / phi)^alpha * lambda * phi / (lambda - phi^2)^2 *
      choose(ns.caps[srvy.ind], 2)
  } 
  
  # Self and parent-offspring pairs between samples
  
  # Set sample pair count to zero
  smp.pr.cnt <- 0
  
  # Loop over all but last survey
  for (srvy.ind.1 in 1:(k - 1)) {
    # Find first survey year
    srvy.yr.1 <- srvy.yrs[srvy.ind.1]
    
    # Loop over surveys with greater indices than first
    for (srvy.ind.2 in (srvy.ind.1 + 1):k) {
      # Increment sample pair count
      smp.pr.cnt <- smp.pr.cnt + 1
      
      # Find second survey year
      srvy.yr.2 <- srvy.yrs[srvy.ind.2]
      
      # Find gap between surveys
      srvy.gap <- srvy.yr.2 - srvy.yr.1
      
      # Find the expected number alive at second sample year
      exp.N.srvy.yr <- exp.N.fin / lambda^(f.year - srvy.yr.2)
      
      # Probability of SPs between samples
      prb.SPs.btn <- phi^srvy.gap / exp.N.srvy.yr
      
      # Number of pairs between samples
      n.prs.btn <- ns.caps[srvy.ind.1] * ns.caps[srvy.ind.2]
      
      # Find expected number of SPs between samples
      exp.ns.SPs.btn[smp.pr.cnt] <- prb.SPs.btn * n.prs.btn
      
      # Probability of POPs where one born between samples
      prb.POPs.brn.btn <- 2 * rho / lambda *
        sum((phi / lambda)^pmax((srvy.gap - 1):0, srvy.gap - alpha - 1)) /
        (exp.N.srvy.yr / lambda^srvy.gap)
      
      # Find expected number of POPs between samples
      exp.ns.POPs.btn[smp.pr.cnt] <- exp.ns.POPs.wtn[srvy.ind.1] * 2 * 
        phi^srvy.gap + prb.POPs.brn.btn * n.prs.btn
    }
  }
  
  # Return expected numbers of kinpairs
  list(
    exp.ns.APs.wtn = exp.ns.APs.wtn, 
    exp.ns.APs.btn = exp.ns.APs.btn, 
    exp.ns.POPs.wtn = exp.ns.POPs.wtn, 
    exp.ns.HSPs.wtn = exp.ns.HSPs.wtn,
    exp.ns.POPs.btn = exp.ns.POPs.btn,
    exp.ns.SPs.btn = exp.ns.SPs.btn,
    exp.ns.SMPs.wtn = exp.ns.SMPs.wtn
  )
}