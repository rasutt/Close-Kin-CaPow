# Function to find numbers of kin pairs among samples within surveys, and
# between pairs of surveys
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
    fm.smp.1.i <- pop.cap.hist[pop.cap.hist[, 3 + smp.ind.1] == 1, ]
    
    # Loop over surveys with greater indices than first
    for (smp.ind.2 in (smp.ind.1 + 1):k) {
      # Find family data for second sample
      fm.smp.2.i <- pop.cap.hist[pop.cap.hist[, 3 + smp.ind.2] == 1, ]
      
      # Increment sample pair count
      smp.pr.cnt <- smp.pr.cnt + 1
      
      # Find number of SPs between pair of samples
      ns.SPs.btn[smp.pr.cnt] <- sum(fm.smp.1.i$ID %in% fm.smp.2.i$ID)
      
      # Find number of POPs between pair of samples
      ns.POPs.btn[smp.pr.cnt] <- sum(
        fm.smp.1.i$mum %in% fm.smp.2.i$ID, 
        fm.smp.1.i$dad %in% fm.smp.2.i$ID,
        fm.smp.2.i$mum %in% fm.smp.1.i$ID, 
        fm.smp.2.i$dad %in% fm.smp.1.i$ID
      )
      
      # Same-mother and self-pairs between pair of samples
      if (any(!is.na(c(fm.smp.1.i$mum, fm.smp.2.i$mum)))) {
        max.mum = max(c(fm.smp.1.i$mum, fm.smp.2.i$mum), na.rm = T)
        mum.tab.1 = tabulate(fm.smp.1.i$mum, max.mum)
        mum.tab.2 = tabulate(fm.smp.2.i$mum, max.mum)
        ns.SMSPs.btn = mum.tab.1 %*% mum.tab.2
      } else {
        ns.SMSPs.btn = 0
      }
        
    
      # Same-father and self-pairs between pair of samples
      if (any(!is.na(c(fm.smp.1.i$dad, fm.smp.2.i$dad)))) {
        max.dad = max(c(fm.smp.1.i$dad, fm.smp.2.i$dad), na.rm = T)
        dad.tab.1 = tabulate(fm.smp.1.i$dad, max.dad)
        dad.tab.2 = tabulate(fm.smp.2.i$dad, max.dad)
        ns.SFSPs.btn = dad.tab.1 %*% dad.tab.2
      } else {
        ns.SFSPs.btn = 0
      }
      
      # Full-sibling and self-pairs
      mums.1 = fm.smp.1.i$mum[!is.na(fm.smp.1.i$mum)]
      mums.2 = fm.smp.2.i$mum[!is.na(fm.smp.2.i$mum)]
      mums = unique(c(mums.1, mums.2))
      dads.1 = fm.smp.1.i$dad[!is.na(fm.smp.1.i$dad)]
      dads.2 = fm.smp.2.i$dad[!is.na(fm.smp.2.i$dad)]
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

# Function to find numbers of offset kin pairs, those among offset set of pairs
# of samples within surveys, and between pairs of surveys
FindNsOKPs <- function(k, n.srvy.prs, std, osiips, osyips) {
  # Family information for offset pairs
  offst.fm.smp.1 = std[osiips[, 1], 1:3]
  offst.fm.smp.2 = std[osiips[, 2], 1:3]
  
  # Create vectors for POPs within each sample, and POPs and SPs between each
  # pair of samples
  ns.POPs.wtn <- ns.HSPs.wtn <- numeric(k)
  ns.SPs.btn <- ns.POPs.btn <- ns.HSPs.btn <- numeric(n.srvy.prs)
  
  # Loop over surveys
  for (smp.ind in 1:k) {
    # Find family data for current sample
    inds = osyips[, 1] == smp.ind - 1 & osyips[, 2] == smp.ind - 1
    fm.smp.1.i = offst.fm.smp.1[inds, ]
    fm.smp.2.i = offst.fm.smp.2[inds, ]
    
    # Find number of POPs within current sample
    ns.POPs.wtn[smp.ind] <- sum(
      fm.smp.1.i$mum == fm.smp.2.i$ID | fm.smp.1.i$ID == fm.smp.2.i$mum |
        fm.smp.1.i$dad == fm.smp.2.i$ID | fm.smp.1.i$ID == fm.smp.2.i$dad,
      na.rm = T
    )
    
    # Find number of HSPs within current sample
    ns.HSPs.wtn[smp.ind] <- sum(
      (fm.smp.1.i$mum == fm.smp.2.i$mum | fm.smp.1.i$dad == fm.smp.2.i$dad) &
        !(fm.smp.1.i$mum == fm.smp.2.i$mum & fm.smp.1.i$dad == fm.smp.2.i$dad),
      na.rm = T 
    )
  }
  
  # Set sample pair count to zero
  smp.pr.cnt <- 0
  
  # Loop over all but last survey
  for (smp.ind.1 in 1:(k - 1)) {
    # Loop over surveys with greater indices than first
    for (smp.ind.2 in (smp.ind.1 + 1):k) {
      # Find family data for second sample
      inds = osyips[, 1] == smp.ind.1 - 1 & osyips[, 2] == smp.ind.2 - 1
      fm.smp.1.i = offst.fm.smp.1[inds, ]
      fm.smp.2.i = offst.fm.smp.2[inds, ]
      
      # Increment sample pair count
      smp.pr.cnt <- smp.pr.cnt + 1
      
      # Find number of SPs between pair of samples
      ns.SPs.btn[smp.pr.cnt] <- sum(fm.smp.1.i$ID == fm.smp.2.i$ID)
      
      # Find number of POPs between pair of samples
      ns.POPs.btn[smp.pr.cnt] <- sum(
        fm.smp.1.i$mum == fm.smp.2.i$ID | fm.smp.1.i$dad == fm.smp.2.i$ID,
        fm.smp.2.i$mum == fm.smp.1.i$ID | fm.smp.2.i$dad == fm.smp.1.i$ID,
        na.rm = T
      )

      # Half-sibling pairs
      ns.HSPs.btn[smp.pr.cnt] = sum(
        (fm.smp.1.i$mum == fm.smp.2.i$mum | 
           fm.smp.1.i$dad == fm.smp.2.i$dad) &
          !(fm.smp.1.i$mum == fm.smp.2.i$mum & 
              fm.smp.1.i$dad == fm.smp.2.i$dad),
        na.rm = T
      )
    }
  }  
  
  # Return numbers of kinpairs
  list(
    wtn = rbind(ns.POPs.wtn, ns.HSPs.wtn),
    btn = rbind(ns.SPs.btn, ns.POPs.btn, ns.HSPs.btn)
  )
}
