FindNsKinPairsPop = function(
  pop.cap.hist, hist.len, srvy.yrs, f.year, n.srvy.prs, k
) {
  # Find total numbers of pairs in whole population within each survey year,
  # and between each pair of survey years
  ns.alv.srvy.yrs = 
    attributes(pop.cap.hist)$N.t.vec[hist.len + srvy.yrs - f.year]
  ns.APs.wtn.pop = choose(ns.alv.srvy.yrs, 2)
  ns.APs.btn.pop = 
    combn(ns.alv.srvy.yrs, 2, function(x) x[1] * x[2])
  
  # Self-pairs in whole population between survey years
  alv_mat = attributes(pop.cap.hist)$alv.mat
  IDs = attributes(pop.cap.hist)$ID
  ns.SPs.btn.pop = numeric(n.srvy.prs)
  pr.cnt = 0
  for (s1 in 1:(k - 1)) {
    srvy.yr.1 = srvy.yrs[s1]
    IDs1 = IDs[alv_mat[, hist.len - f.year + srvy.yr.1] == 1]
    for (s2 in (s1 + 1):k) {
      pr.cnt = pr.cnt + 1
      srvy.yr.2 = srvy.yrs[s2]
      IDs2 = IDs[alv_mat[, hist.len - f.year + srvy.yr.2] == 1]
      ns.SPs.btn.pop[pr.cnt] = sum(IDs1 %in% IDs2)
    }
  }
  
  # Same-mother pairs whole population in survey years
  ns.SMPs.wtn.pop = numeric(k)
  for (s.ind in 1:k) {
    alv.s.yr = alv_mat[, hist.len - f.year + srvy.yrs[s.ind]] == 1
    ns.SMPs.wtn.pop[s.ind] = 
      sum(choose(table(attributes(pop.cap.hist)$mum[alv.s.yr]), 2))
  }
  
  # Return as list
  list(
    ns.APs.wtn.pop = ns.APs.wtn.pop,
    ns.APs.btn.pop = ns.APs.btn.pop,
    ns.SPs.btn.pop = ns.SPs.btn.pop,
    ns.SMPs.wtn.pop = ns.SMPs.wtn.pop
  )
}