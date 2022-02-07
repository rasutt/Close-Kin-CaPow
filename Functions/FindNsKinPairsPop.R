FindNsKinPairsPop = function(
  N.s.yrs, alv_s_yrs, ID, mum, dad, k
) {
  # Find total numbers of pairs in whole population within each survey year,
  # and between each pair of survey years
  ns.APs.wtn.pop = choose(N.s.yrs, 2)
  ns.APs.btn.pop = combn(N.s.yrs, 2, function(x) x[1] * x[2])
  
  # Self-pairs in whole population between survey years
  ns.SPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
    sum(ID[alv_s_yrs[, s.inds[1]]] %in% ID[alv_s_yrs[, s.inds[2]]])
  }))
  
  # List of frequency tables of mums and dads in each survey year
  mum_vec_lst = lapply(1:k, function(s.ind) mum[alv_s_yrs[, s.ind]])
  dad_vec_lst = lapply(1:k, function(s.ind) dad[alv_s_yrs[, s.ind]])
  max_mum = max(unlist(mum_vec_lst), na.rm = T)
  max_dad = max(unlist(dad_vec_lst), na.rm = T)
  mum_tab_lst = lapply(1:k, function(s.ind) {
    tabulate(mum_vec_lst[[s.ind]], max_mum)
  })
  dad_tab_lst = lapply(1:k, function(s.ind) {
    tabulate(dad_vec_lst[[s.ind]], max_dad)
  })
  
  # Same-mother pairs whole population in survey years
  ns.SMPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(mum_tab_lst[[s.ind]], 2))
  })
  
  # Same-father pairs whole population in survey years
  ns.SFPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(dad_tab_lst[[s.ind]], 2))
  })
  
  # Full-sibling pairs whole population in survey years
  ns.FSPs.wtn.pop = sapply(1:k, function(s.ind) {
    sum(choose(table(mum[alv_s_yrs[, s.ind]], dad[alv_s_yrs[, s.ind]]), 2))
  })
  
  # Same-mother pairs in whole population between survey years
  ns.SMPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
    mum_tab_lst[[s.inds[1]]] %*% mum_tab_lst[[s.inds[2]]]
  }))
  
  # Same-father pairs in whole population between survey years
  ns.SFPs.btn.pop = as.vector(combn(1:k, 2, function(s.inds) {
    dad_tab_lst[[s.inds[1]]] %*% dad_tab_lst[[s.inds[2]]]
  }))
  
  # Return as list
  list(
    ns.APs.wtn.pop = ns.APs.wtn.pop,
    ns.APs.btn.pop = ns.APs.btn.pop,
    ns.SPs.btn.pop = ns.SPs.btn.pop,
    ns.SMPs.wtn.pop = ns.SMPs.wtn.pop,
    ns.SFPs.wtn.pop = ns.SFPs.wtn.pop,
    ns.FSPs.wtn.pop = ns.FSPs.wtn.pop,
    ns.HSPs.wtn.pop = ns.SMPs.wtn.pop + ns.SFPs.wtn.pop - 2 * ns.FSPs.wtn.pop,
    ns.SMPs.btn.pop = ns.SMPs.btn.pop,
    ns.SFPs.btn.pop = ns.SFPs.btn.pop
  )
}