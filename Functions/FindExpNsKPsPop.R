# Function to find expected numbers of kin-pairs in whole population
FindEstNsKPsPop = function(
    exp.N.t, s.yr.inds, phi, lambda, alpha, srvy.yrs, k
) {
  # Intermediate results
  
  # Expected population sizes in survey years
  exp.N.s.yrs = exp.N.t[s.yr.inds]
  # Probability that two animals are new-born independently
  prb.nw.brn.sq = (1 - phi / lambda)^2
  # Probability that an animal is mature
  prb.mtr = (lambda / phi)^alpha
  # Lambda minus phi-squared
  lmb.m.ph.sq = lambda - phi^2
  
  # Expected numbers of kin-pairs
  
  # Approximation for expected number of pairs of animals in survey years
  exp.ns.APs.wtn = choose(exp.N.s.yrs, 2)
  
  # Parent-offspring pairs within survey years
  exp.ns.POPs.wtn = 2 * exp.N.s.yrs * phi * (lambda - phi) / lmb.m.ph.sq
  
  # Same-mother pairs within survey years
  exp.ns.SMPs.wtn = 2 * exp.N.s.yrs * prb.nw.brn.sq * prb.mtr * lambda *
    phi^2 / lmb.m.ph.sq^2
  
  # Same-father pairs within survey years
  exp.ns.SFPs.wtn = prb.nw.brn.sq * prb.mtr * lambda * 
    (exp.N.s.yrs / lmb.m.ph.sq * (2 * phi^3 / lmb.m.ph.sq + lambda) - 
       lambda^alpha / (1 - phi^2))
  
  # Full-sibling pairs within survey years
  exp.ns.FSPs.wtn = 4 * prb.nw.brn.sq * prb.mtr^2 * phi^4 / 
    (lambda - phi^3) / (1 - phi^2)
  
  # Half-sibling pairs within survey years
  exp.ns.HSPs.wtn = exp.ns.SMPs.wtn + exp.ns.SFPs.wtn - 2 * exp.ns.FSPs.wtn
  
  # All pairs within surveys and between pairs of surveys
  exp.ns.APs.btn = as.vector(combn(exp.N.s.yrs, 2, function(x) x[1] * x[2]))
  
  # Self-pairs between pairs of survey years
  exp.ns.SPs.btn = as.vector(combn(1:k, 2, function(s.inds) {
    phi^(srvy.yrs[s.inds[2]] - srvy.yrs[s.inds[1]]) * 
      exp.N.s.yrs[s.inds[1]]
  }))
  
  # Parent-offspring pairs between survey years
  exp.ns.POPs.btn = as.vector(combn(1:k, 2, function(s.inds) {
    t1 = srvy.yrs[s.inds[1]]
    t2 = srvy.yrs[s.inds[2]]
    4 * exp.N.s.yrs[s.inds[2]] * phi * (lambda - phi) / lmb.m.ph.sq *
      (phi / lambda)^(t2 - t1) +
      2 * (1 - phi / lambda) * exp.N.s.yrs[s.inds[2]] * (phi / lambda)^t2 * 
      (((phi / lambda)^(-(t1 + 1)) -
          (phi / lambda)^(-(min(t1 + alpha, t2) + 1))) /
         (1 - (phi / lambda)^(-1)) + 
         ifelse(
           t1 + alpha < t2, 
           (t2 - (t1 + alpha)) * (phi / lambda)^(-(t1 + alpha)),
           0
         ))
  }))
  
  # Birth rate among mature females
  beta = 2 * (1 - phi/lambda) * (lambda/phi)^alpha
  
  # Same-mother pairs between survey years
  exp.ns.SMPs.btn = as.vector(combn(1:k, 2, function(s.inds) {
    # Gap between survey years
    srvy.gap = srvy.yrs[s.inds[2]] - srvy.yrs[s.inds[1]]
    
    # Same-mother pairs
    2 * exp.ns.SMPs.wtn[s.inds[1]] * phi^srvy.gap +
      2 * srvy.gap * exp.N.s.yrs[s.inds[2]] * beta * (1 - phi/lambda) * 
      lambda / lmb.m.ph.sq * (phi/lambda)^srvy.gap
  }))

  # Return as list
  list(
    wtn = cbind(
      N.s.yrs = exp.N.s.yrs, APs = exp.ns.APs.wtn, POPs = exp.ns.POPs.wtn, 
      SMPs = exp.ns.SMPs.wtn, SFPs = exp.ns.SFPs.wtn, FSPs = exp.ns.FSPs.wtn, 
      HSPs = exp.ns.HSPs.wtn
    ),
    btn = cbind(
      APs = exp.ns.APs.btn, POPs = exp.ns.POPs.btn, SPs = exp.ns.SPs.btn, 
      SMPs = exp.ns.SMPs.btn
    )
  )
}
