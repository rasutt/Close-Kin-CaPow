# Function to find kin-pair probabilities
FindKPProbs = function(
  exp.N.srvy.yrs, exp.ns.APs.wtn, phi, lambda, alpha, srvy.yrs, k
) {
  # Intermediate results
  
  # Probability that an animal is new-born
  prb.nw.brn = 1 - phi / lambda
  # Probability that an animal is mature
  prb.mtr = (lambda / phi)^alpha
  # Lambda minus phi-squared
  lmb.m.ph.sq = lambda - phi^2
  
  # Kin-pair probabilities
  
  # Same-mother pairs within survey years
  prbs.SMPs.wtn = 4 / (exp.N.srvy.yrs - 1) * prb.nw.brn^2 * prb.mtr * lambda *
    phi^2 / lmb.m.ph.sq^2
  
  # Same-father pairs within survey years
  prbs.SFPs.wtn = prb.nw.brn^2 * prb.mtr * lambda * 
    ((2 * phi^3 / lmb.m.ph.sq^2 + lambda / lmb.m.ph.sq) / 
       (exp.N.srvy.yrs - 1) - lambda^alpha / (1 - phi^2) / exp.ns.APs.wtn)
  
  # Full-sibling pairs within survey years
  prbs.FSPs.wtn = 4 / exp.ns.APs.wtn * prb.nw.brn^2 * prb.mtr^2 * phi^4 / 
    (lambda - phi^3) / (1 - phi^2)
  
  # Half-sibling pairs within survey years
  prbs.HSPs.wtn = prbs.SMPs.wtn + prbs.SFPs.wtn - 2 * prbs.FSPs.wtn
  
  # Self-pairs between pairs of survey years
  prbs.SPs.btn = as.vector(combn(1:k, 2, function(s.pr.inds) {
    phi^(srvy.yrs[s.pr.inds[2]] - srvy.yrs[s.pr.inds[1]]) / 
      exp.N.srvy.yrs[s.pr.inds[2]]
  }))
  
  # Same-mother pairs between survey years
  prbs.SMPs.btn = as.vector(combn(1:k, 2, function(s.pr.inds) {
    prbs.SMPs.wtn[s.pr.inds[1]] * 
      phi^(srvy.yrs[s.pr.inds[2]] - srvy.yrs[s.pr.inds[1]]) *
      (lmb.m.ph.sq / phi^2 * 
         (srvy.yrs[s.pr.inds[2]] - srvy.yrs[s.pr.inds[1]]) + 2)
  }))
  
  # Return as list
  list(
    wtn = rbind(prbs.SMPs.wtn, prbs.SFPs.wtn, prbs.FSPs.wtn, prbs.HSPs.wtn),
    btn = rbind(prbs.SPs.btn, prbs.SMPs.btn)
  )
}