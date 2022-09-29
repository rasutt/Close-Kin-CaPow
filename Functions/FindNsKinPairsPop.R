# Function to find numbers of kin pairs in population within, and between pairs
# of, survey years
FindNsKinPairsPop = function(pop.cap.hist, s.yr.inds, k) {
  # Total numbers of pairs
  N.s.yrs = attributes(pop.cap.hist)$N.t.vec[s.yr.inds]
  ns.APs.wtn.pop = choose(N.s.yrs, 2)
  ns.APs.btn.pop = combn(N.s.yrs, 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2])
  
  # Which animals alive in survey years 
  alv.s.yrs = attributes(pop.cap.hist)$alv.s.yrs

  # Self-pairs between survey years
  ID = attributes(pop.cap.hist)$ID
  ns.SPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
    sum(ID[alv.s.yrs[, s.inds[1]]] %in% ID[alv.s.yrs[, s.inds[2]]])
  }))
  
  # List of frequency tables of mums and dads in each survey year
  mum = attributes(pop.cap.hist)$mum
  dad = attributes(pop.cap.hist)$dad
  max.mum = max(mum, na.rm = T)
  max.dad = max(dad, na.rm = T)
  mum.tab.lst = lapply(1:k, function(s.ind) {
    tabulate(mum[alv.s.yrs[, s.ind]], max.mum)
  })
  dad.tab.lst = lapply(1:k, function(s.ind) {
    tabulate(dad[alv.s.yrs[, s.ind]], max.dad)
  })
  
  # Parent-offspring pairs within survey years
  ns.POPs.wtn = sapply(1:k, function(s.ind) {
    sum(
      mum[alv.s.yrs[, s.ind]] %in% ID[alv.s.yrs[, s.ind]], 
      dad[alv.s.yrs[, s.ind]] %in% ID[alv.s.yrs[, s.ind]]
    )
  })
  
  # Same-mother pairs within survey years
  ns.SMPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(mum.tab.lst[[s.ind]], 2))
  })
  
  # Same-father pairs within survey years
  ns.SFPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(dad.tab.lst[[s.ind]], 2))
  })
  
  # Full-sibling pairs within survey years
  ns.FSPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(table(mum[alv.s.yrs[, s.ind]], dad[alv.s.yrs[, s.ind]]), 2))
  })
  
  # Half-sibling pairs within survey years
  ns.HSPs.wtn.pop = ns.SMPs.wtn.pop + ns.SFPs.wtn.pop - 2 * ns.FSPs.wtn.pop
  
  # Same-mother pairs between survey years
  ns.SMPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
    mum.tab.lst[[s.inds[1]]] %*% mum.tab.lst[[s.inds[2]]]
  }))
  
  # # Same-father pairs between survey years
  # ns.SFPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
  #   dad.tab.lst[[s.inds[1]]] %*% dad.tab.lst[[s.inds[2]]]
  # }))
  
  # Return as list
  list(
    wtn = rbind(
      ns.APs.wtn.pop, ns.POPs.wtn, ns.SMPs.wtn.pop, ns.SFPs.wtn.pop, 
      ns.FSPs.wtn.pop, ns.HSPs.wtn.pop
    ),
    btn = rbind(
      ns.APs.btn.pop, ns.SPs.btn.pop, ns.SMPs.btn.pop
      # , ns.SFPs.btn.pop
    )
  )
}