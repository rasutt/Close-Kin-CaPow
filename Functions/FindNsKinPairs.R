FindNsKinPairs <- function(k, n.srvy.prs, pop.cap.hist) {
  # Create vectors for POPs within each sample, and POPs and SPs between each
  # pair of samples
  ns.POPs.wtn <- ns.HSPs.wtn <- numeric(k)
  ns.SPs.btn <- ns.POPs.btn <- ns.HSPs.btn <- numeric(n.srvy.prs)
  
  # Loop over surveys
  for (smp.ind in 1:k) {
    # Find family data for current sample
    fam.samp <- pop.cap.hist[pop.cap.hist[, 3 + smp.ind] == 1, ]

    # Find number of POPs within current sample
    ns.POPs.wtn[smp.ind] <- sum(
      fam.samp$mum %in% fam.samp$ID, 
      fam.samp$dad %in% fam.samp$ID
    )
    
    # Find number of HSPs within current sample. Selecting before calling choose
    # is 10x faster
    tbl = table(fam.samp$mum, fam.samp$dad)
    ns.HSPs.wtn[smp.ind] <- sum(
      choose(table(fam.samp$mum), 2),
      choose(table(fam.samp$dad), 2),
      -2 * choose(tbl[tbl > 1], 2)
    )
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
      
      # Find number of SPs between pair of samples
      ns.SPs.btn[smp.pr.cnt] <- sum(fam.smp.1$ID %in% fam.smp.2$ID)
      
      # Find number of POPs between pair of samples
      ns.POPs.btn[smp.pr.cnt] <- sum(
        fam.smp.1$mum %in% fam.smp.2$ID, 
        fam.smp.1$dad %in% fam.smp.2$ID,
        fam.smp.2$mum %in% fam.smp.1$ID, 
        fam.smp.2$dad %in% fam.smp.1$ID
      )
      
      # Same-mother and self-pairs between pair of samples
      max.mum = max(c(fam.smp.1$mum, fam.smp.2$mum), na.rm = T)
      mum.tab.1 = tabulate(fam.smp.1$mum, max.mum)
      mum.tab.2 = tabulate(fam.smp.2$mum, max.mum)
      ns.SMSPs.btn = mum.tab.1 %*% mum.tab.2
    
      # Same-father and self-pairs between pair of samples
      max.dad = max(c(fam.smp.1$dad, fam.smp.2$dad), na.rm = T)
      dad.tab.1 = tabulate(fam.smp.1$dad, max.dad)
      dad.tab.2 = tabulate(fam.smp.2$dad, max.dad)
      ns.SFSPs.btn = dad.tab.1 %*% dad.tab.2
      
      # Full-sibling and self-pairs
      mums.1 = fam.smp.1$mum[!is.na(fam.smp.1$mum)]
      mums.2 = fam.smp.2$mum[!is.na(fam.smp.2$mum)]
      mums = unique(c(mums.1, mums.2))
      dads.1 = fam.smp.1$dad[!is.na(fam.smp.1$dad)]
      dads.2 = fam.smp.2$dad[!is.na(fam.smp.2$dad)]
      dads = unique(c(dads.1, dads.2))
      
      # Two-way frequency tables of mothers and fathers of animals in each survey
      # year, including all possible parents in levels so that dimensions match
      tbl.1 = table(
        factor(mums.1, levels = mums), factor(dads.1, levels = dads)
      )
      tbl.2 = table(
        factor(mums.2, levels = mums), factor(dads.2, levels = dads)
      )
      
      # Full-sibling pairs, including self-pairs with known parents
      ns.FSSPs.btn = sum(tbl.1 * tbl.2)
      
      # Half-sibling pairs
      ns.HSPs.btn[smp.pr.cnt] = ns.SMSPs.btn + ns.SFSPs.btn - 2 * ns.FSSPs.btn
    }
  }  

  # Return numbers of kinpairs
  list(
    wtn = rbind(ns.POPs.wtn, ns.HSPs.wtn),
    btn = rbind(ns.SPs.btn, ns.POPs.btn, ns.HSPs.btn)
  )
}
