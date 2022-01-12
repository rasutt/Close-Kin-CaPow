FindNsKinPairs <- function(k, n.srvy.prs, pop.cap.hist, ns.alv.srvy.yrs) {
  # Find total numbers of pairs within each survey year, and between each pair
  # of survey years
  ns.APs.wtn = choose(ns.alv.srvy.yrs, 2)
  ns.APs.btn = combn(ns.alv.srvy.yrs, 2, function(x) x[1] * x[2])
  
  # Create vectors for POPs within each sample, and POPs and SPs between each
  # pair of samples
  ns.POPs.wtn <- ns.HSPs.wtn <- ns.SMPs.wtn <- numeric(k)
  ns.POPs.btn <- ns.SPs.btn <- numeric(n.srvy.prs)
  
  # Loop over surveys
  for (smp.ind in 1:k) {
    # Find family data for current sample
    fam.samp <- pop.cap.hist[pop.cap.hist[, 3 + smp.ind] == 1, ]

    # Find number of POPs within current sample
    ns.POPs.wtn[smp.ind] <- sum(
      fam.samp$mum %in% fam.samp$ID, 
      fam.samp$dad %in% fam.samp$ID
    )
    
    # Find number of HSPs within current sample
    ns.HSPs.wtn[smp.ind] <- sum(
      choose(table(fam.samp$mum), 2),
      choose(table(fam.samp$dad), 2),
      -2 * choose(table(fam.samp$mum, fam.samp$dad), 2)
    )
    
    # Same mother pairs
    ns.SMPs.wtn[smp.ind] <- sum(choose(table(fam.samp$mum), 2))
  }
  
  # Set sample pair count to zero
  smp.pr.cnt <- 0
  
  # Loop over all but last survey
  for (smp.ind.1 in 1:(k - 1)) {
    # Find family data for first sample
    fam.smp.1 <- pop.cap.hist[pop.cap.hist[, 3 + smp.ind.1] == 1, ]
    
    # Loop over surveys with greater indices than first
    for (smp.ind.2 in (smp.ind.1 + 1):k) {
      # Find family data for second sample
      fam.smp.2 <- pop.cap.hist[pop.cap.hist[, 3 + smp.ind.2] == 1, ]
      
      # Increment sample pair count
      smp.pr.cnt <- smp.pr.cnt + 1
      
      # Find number of POPs between pair of samples
      ns.POPs.btn[smp.pr.cnt] <- sum(
        fam.smp.1$mum %in% fam.smp.2$ID, 
        fam.smp.1$dad %in% fam.smp.2$ID,
        fam.smp.2$mum %in% fam.smp.1$ID, 
        fam.smp.2$dad %in% fam.smp.1$ID
      )
      
      # Find number of SPs between pair of samples
      ns.SPs.btn[smp.pr.cnt] <- sum(fam.smp.1$ID %in% fam.smp.2$ID)
    }
  }  

  # Return numbers of kinpairs
  list(
    ns.APs.wtn = ns.APs.wtn,
    ns.APs.btn = ns.APs.btn,
    ns.POPs.wtn = ns.POPs.wtn, 
    ns.HSPs.wtn = ns.HSPs.wtn,
    ns.SMPs.wtn = ns.SMPs.wtn,
    ns.POPs.btn = ns.POPs.btn,
    ns.SPs.btn = ns.SPs.btn
  )
}
