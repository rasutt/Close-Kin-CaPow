# Function to find expected numbers of kin-pairs in whole population
FindEstNsKPsPop = function(
    exp.N.t, s.yr.inds, phi, lambda, alpha, srvy.yrs, k
) {
  # Intermediate results
  
  # Expected population sizes in survey years
  exp.N.srvy.yrs = exp.N.t[s.yr.inds]
  # Probability that two animals are new-born independently
  prb.nw.brn.sq = (1 - phi / lambda)^2
  # Probability that an animal is mature
  prb.mtr = (lambda / phi)^alpha
  # Lambda minus phi-squared
  lmb.m.ph.sq = lambda - phi^2
  
  # Expected numbers of kin-pairs
  
  # Approximation for expected number of pairs of animals in survey years
  exp.ns.APs.wtn = choose(exp.N.srvy.yrs, 2)
  
  # Parent-offspring pairs within survey years
  exp.ns.POPs.wtn = 2 * exp.N.srvy.yrs * phi * (lambda - phi) / lmb.m.ph.sq
  
  # Same-mother pairs within survey years
  exp.ns.SMPs.wtn = 2 * exp.N.srvy.yrs * prb.nw.brn.sq * prb.mtr * lambda *
    phi^2 / lmb.m.ph.sq^2
  
  # Same-father pairs within survey years
  exp.ns.SFPs.wtn = prb.nw.brn.sq * prb.mtr * lambda * 
    (exp.N.srvy.yrs / lmb.m.ph.sq * (2 * phi^3 / lmb.m.ph.sq + lambda) - 
       lambda^alpha / (1 - phi^2))
  
  # Full-sibling pairs within survey years
  exp.ns.FSPs.wtn = 4 * prb.nw.brn.sq * prb.mtr^2 * phi^4 / 
    (lambda - phi^3) / (1 - phi^2)
  
  # Half-sibling pairs within survey years
  exp.ns.HSPs.wtn = exp.ns.SMPs.wtn + exp.ns.SFPs.wtn - 2 * exp.ns.FSPs.wtn
  
  # All pairs within surveys and between pairs of surveys
  exp.ns.APs.btn = as.vector(combn(exp.N.srvy.yrs, 2, function(x) x[1] * x[2]))
  
  # Self-pairs between pairs of survey years
  exp.ns.SPs.btn = as.vector(combn(1:k, 2, function(s.pr.inds) {
    phi^(srvy.yrs[s.pr.inds[2]] - srvy.yrs[s.pr.inds[1]]) * 
      exp.N.srvy.yrs[s.pr.inds[1]]
  }))
  
  # Same-mother pairs between survey years
  exp.ns.POPs.btn = as.vector(combn(1:k, 2, function(s.pr.inds) {
    s.gap.lim = s.pr.inds[2] - s.pr.inds[1] - 1 - alpha
    4 * exp.N.srvy.yrs[s.pr.inds[2]] * phi * (lambda - phi) / lmb.m.ph.sq *
      (phi / lambda)^(s.pr.inds[2] - s.pr.inds[1]) +
      2 * (1 - phi / lambda) * exp.N.srvy.yrs[s.pr.inds[2]] * 
      sum((phi / lambda)^pmax(0:(s.gap.lim + alpha), s.gap.lim))
  }))
  
  # Same-mother pairs between survey years
  exp.ns.SMPs.btn = as.vector(combn(1:k, 2, function(s.pr.inds) {
    exp.ns.SMPs.wtn[s.pr.inds[1]] * 
      phi^(srvy.yrs[s.pr.inds[2]] - srvy.yrs[s.pr.inds[1]]) *
      (lmb.m.ph.sq / phi^2 * 
         (srvy.yrs[s.pr.inds[2]] - srvy.yrs[s.pr.inds[1]]) + 2)
  }))
  
  # Return as list
  list(
    wtn = cbind(
      exp.N.srvy.yrs, exp.ns.APs.wtn, exp.ns.POPs.wtn, exp.ns.SMPs.wtn, 
      exp.ns.SFPs.wtn, exp.ns.FSPs.wtn, exp.ns.HSPs.wtn
    ),
    btn = cbind(
      exp.ns.APs.btn, exp.ns.POPs.btn, exp.ns.SPs.btn, exp.ns.SMPs.btn
    )
  )
}
