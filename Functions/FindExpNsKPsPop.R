# Function to find expected numbers of kin-pairs in whole population
FindEstNsKPsPop = function(
    exp.N.t, s.yr.inds, phi, lambda, alpha, srvy.yrs, k
) {
  # Intermediate results
  
  # Expected population sizes in survey years
  exp.N.s.yrs = exp.N.t[s.yr.inds]
  # Probability not new-born (phi over lambda)
  p.o.l = phi / lambda
  # Probability that two animals are new-born independently
  prb.nw.brn.sq = (1 - p.o.l)^2
  # Reciprocal of probability that an animal is mature
  rcl.prb.mtr = (lambda / phi)^alpha
  # Lambda minus phi-squared
  lmb.m.ph.sq = lambda - phi^2
  # Birth rate among mature females
  beta = 2 * (1 - p.o.l) * (lambda/phi)^alpha
  
  # Expected numbers of kin-pairs
  
  # Approximation for expected number of pairs of animals in survey years
  exp.ns.APs.wtn = choose(exp.N.s.yrs, 2)
  
  # Parent-offspring pairs within survey years
  exp.ns.POPs.wtn = 2 * exp.N.s.yrs * phi * (lambda - phi) / lmb.m.ph.sq
  
  # Same-mother pairs within survey years
  exp.ns.SMPs.wtn = exp.N.s.yrs * beta * (lambda - phi) * phi^2 / lmb.m.ph.sq^2
  
  # Same-father pairs within survey years
  exp.ns.SFPs.diff.b.yrs.wtn = phi * exp.ns.SMPs.wtn
  exp.ns.SFPs.same.b.yr.wtn = beta^2 * phi^(alpha + 1) / 4 * 
    (exp.N.s.yrs / (lambda^(alpha - 1) * lmb.m.ph.sq) - 1 / (1 - phi^2))
  exp.ns.SFPs.wtn = exp.ns.SFPs.diff.b.yrs.wtn + exp.ns.SFPs.same.b.yr.wtn
  
  # Full-sibling pairs within survey years
  exp.ns.FSPs.wtn = 2 * beta * rcl.prb.mtr * (lambda - phi) * phi^4 / 
    (lambda * (lambda - phi^3) * (1 - phi^2))
  
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
      (p.o.l)^(t2 - t1) +
      2 * (1 - p.o.l) * exp.N.s.yrs[s.inds[2]] * (p.o.l)^t2 * 
      (((p.o.l)^(-(t1 + 1)) -
          (p.o.l)^(-(min(t1 + alpha, t2) + 1))) /
         (1 - (p.o.l)^(-1)) + 
         ifelse(
           t1 + alpha < t2, 
           (t2 - (t1 + alpha)) * (p.o.l)^(-(t1 + alpha)),
           0
         ))
  }))
  
  # Same-mother pairs between survey years
  exp.ns.SMPs.btn = as.vector(combn(1:k, 2, function(s.inds) {
    # Gap between survey years
    srvy.gap = srvy.yrs[s.inds[2]] - srvy.yrs[s.inds[1]]
    
    # Same-mother pairs
    2 * exp.ns.SMPs.wtn[s.inds[1]] * phi^srvy.gap +
      srvy.gap * exp.N.s.yrs[s.inds[2]] * beta * (1 - p.o.l) * 
      lambda / lmb.m.ph.sq * (p.o.l)^srvy.gap
  }))
  
  # Same-mother pairs between surveys, ages known (five and zero)
  exp.ns.SMPs.kwn.age.btn = as.vector(combn(1:k, 2, function(s.inds) {
    exp.N.s.yrs[s.inds[2]] * (1 - p.o.l) * beta * phi^5 *
      (p.o.l)^(srvy.yrs[s.inds[2]] - srvy.yrs[s.inds[1]] + 5)
  }))
  
  # Same-father pairs between survey years
  exp.ns.SFPs.btn = phi * exp.ns.SMPs.btn + 
    as.vector(combn(1:k, 2, function(s.inds) {
      2 * phi^(srvy.yrs[s.inds[2]] - srvy.yrs[s.inds[1]]) * 
        exp.ns.SFPs.same.b.yr.wtn[s.inds[1]]
    }))
  
  # Full-sibling pairs between survey years
  exp.ns.FSPs.btn = 
    as.vector(combn(1:k, 2, function(s.inds) {
      t1 = srvy.yrs[s.inds[1]]
      t2 = srvy.yrs[s.inds[2]]
      2 * exp.ns.FSPs.wtn[s.inds[1]] * phi^(t2 - t1) +
      2 * beta * (1 - p.o.l) * rcl.prb.mtr * phi^(t1 + t2 + 1) * 
        (p.o.l^(t1 + 1) - p.o.l^(t2 + 1)) / 
        ((phi^3 / lambda)^t1 * (1 - phi^3 / lambda) * (1 - p.o.l))
    }))
  
  # Half-sibling pairs between surveys
  exp.ns.HSPs.btn = exp.ns.SMPs.btn + exp.ns.SFPs.btn - 2 * exp.ns.FSPs.btn
  
  # Return as list
  list(
    wtn = cbind(
      N.s.yrs = exp.N.s.yrs, APs = exp.ns.APs.wtn, POPs = exp.ns.POPs.wtn, 
      SMPs = exp.ns.SMPs.wtn, SFPs = exp.ns.SFPs.wtn, FSPs = exp.ns.FSPs.wtn, 
      HSPs = exp.ns.HSPs.wtn
    ),
    btn = cbind(
      APs = exp.ns.APs.btn, POPs = exp.ns.POPs.btn, SPs = exp.ns.SPs.btn, 
      SMPs = exp.ns.SMPs.btn, SMPs.kwn.age = exp.ns.SMPs.kwn.age.btn,
      SFPs = exp.ns.SFPs.btn, FSPs = exp.ns.FSPs.btn, HSPs = exp.ns.HSPs.btn
    )
  )
}
