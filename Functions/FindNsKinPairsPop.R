# Function to find numbers of kin pairs in population within, and between pairs
# of, survey years
FindNsKinPairsPop = function(pop.cap.hist, s.yr.inds, k) {
  # Total numbers of pairs
  N.s.yrs = attributes(pop.cap.hist)$N.t.vec[s.yr.inds]
  ns.APs.wtn.pop = choose(N.s.yrs, 2)
  ns.APs.btn.pop = combn(N.s.yrs, 2, function(x) x[1] * x[2])
  
  # Self-pairs between survey years
  alv_s_yrs = attributes(pop.cap.hist)$alv.mat[, s.yr.inds]
  ID = attributes(pop.cap.hist)$ID
  ns.SPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
    sum(ID[alv_s_yrs[, s.inds[1]]] %in% ID[alv_s_yrs[, s.inds[2]]])
  }))
  
  # List of frequency tables of mums and dads in each survey year
  mum = attributes(pop.cap.hist)$mum
  dad = attributes(pop.cap.hist)$dad
  max_mum = max(mum, na.rm = T)
  max_dad = max(dad, na.rm = T)
  mum_tab_lst = lapply(1:k, function(s.ind) {
    tabulate(mum[alv_s_yrs[, s.ind]], max_mum)
  })
  dad_tab_lst = lapply(1:k, function(s.ind) {
    tabulate(dad[alv_s_yrs[, s.ind]], max_dad)
  })
  
  # Same-mother pairs within survey years
  ns.SMPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(mum_tab_lst[[s.ind]], 2))
  })
  
  # Same-father pairs within survey years
  ns.SFPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(dad_tab_lst[[s.ind]], 2))
  })
  
  # Full-sibling pairs within survey years
  ns.FSPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(table(mum[alv_s_yrs[, s.ind]], dad[alv_s_yrs[, s.ind]]), 2))
  })
  
  # Half-sibling pairs within survey years
  ns.HSPs.wtn.pop = ns.SMPs.wtn.pop + ns.SFPs.wtn.pop - 2 * ns.FSPs.wtn.pop
  
  # Same-mother pairs between survey years
  ns.SMPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
    mum_tab_lst[[s.inds[1]]] %*% mum_tab_lst[[s.inds[2]]]
  }))
  
  # Same-father pairs between survey years
  ns.SFPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
    dad_tab_lst[[s.inds[1]]] %*% dad_tab_lst[[s.inds[2]]]
  }))
  
  # Return as list
  list(
    wtn = rbind(
      ns.APs.wtn.pop, ns.SMPs.wtn.pop, ns.SFPs.wtn.pop, ns.FSPs.wtn.pop,
      ns.HSPs.wtn.pop
    ),
    btn = rbind(
      ns.APs.btn.pop, ns.SPs.btn.pop, ns.SMPs.btn.pop, ns.SFPs.btn.pop
    )
  )
}