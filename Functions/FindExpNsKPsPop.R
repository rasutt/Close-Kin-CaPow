# Function to find expected numbers of kin-pairs in whole population
FindEstNsKPsPop = function(
    exp.N.t, s.yr.inds, phi, rho, lambda, alpha, srvy.yrs, k
) {
  ## Intermediate results
  
  # Expected population sizes in survey years
  exp.N.s.yrs = exp.N.t[s.yr.inds]
  # Lambda minus phi-squared
  lmb.m.ph.sq = lambda - phi^2
  # Probability not new-born (phi over lambda)
  p.o.l = phi / lambda
  # Reciprocal of probability that an animal is mature
  rcl.prb.mtr = (lambda / phi)^alpha
  # Birth rate among mature females
  beta = 2 * (1 - p.o.l) * rcl.prb.mtr
  
  ### Expected numbers of kin-pairs
  
  ## Within surveys
  
  # Approximation for expected number of pairs of animals in survey years
  exp.ns.APs.wtn = choose(exp.N.s.yrs, 2)
  
  # Parent-offspring pairs within survey years
  exp.ns.POPs.wtn = 2 * exp.N.s.yrs * phi * rho / lmb.m.ph.sq
  
  # Same-mother pairs within survey years
  exp.ns.SMPs.wtn = exp.N.s.yrs * beta * rho * phi^2 / lmb.m.ph.sq^2
  
  # Same-father pairs within survey years, split into same and different birth
  # years as used separately later
  exp.ns.SFPs.diff.b.yrs.wtn = phi * exp.ns.SMPs.wtn
  exp.ns.SFPs.same.b.yr.wtn = beta^2 * phi^(alpha + 1) / 4 * 
    (exp.N.s.yrs / (lambda^(alpha - 1) * lmb.m.ph.sq) - 1 / (1 - phi^2))
  exp.ns.SFPs.wtn = exp.ns.SFPs.diff.b.yrs.wtn + exp.ns.SFPs.same.b.yr.wtn
  
  # Full-sibling pairs within survey years, constant over time but gets repeated
  # by cbind when returned
  exp.ns.FSPs.wtn = 2 * beta * rcl.prb.mtr * rho * phi^4 / 
    (lambda * (lambda - phi^3) * (1 - phi^2))
  
  # Half-sibling pairs within survey years
  exp.ns.HSPs.wtn = exp.ns.SMPs.wtn + exp.ns.SFPs.wtn - 2 * exp.ns.FSPs.wtn
  
  ## Between surveys
  
  # Survey-pair indices
  s.pr.inds = combn(k, 2)
  s.inds.1 = s.pr.inds[1, ]
  s.inds.2 = s.pr.inds[2, ]
  exp.N.s.yrs.1 = exp.N.s.yrs[s.inds.1]
  exp.N.s.yrs.2 = exp.N.s.yrs[s.inds.2]
  s.yrs.1 = srvy.yrs[s.inds.1]
  s.yrs.2 = srvy.yrs[s.inds.2]
  s.gaps = s.yrs.2 - s.yrs.1
  p.t.s.gs = phi^s.gaps
  
  # All pairs between pairs of surveys
  exp.ns.APs.btn = exp.N.s.yrs.1 * exp.N.s.yrs.2

  # Self-pairs between pairs of survey years
  exp.ns.SPs.btn = p.t.s.gs * exp.N.s.yrs.1
  
  # Parent-offspring pairs between survey years
  s.yrs.1.p.a = s.yrs.1 + alpha
  exp.ns.POPs.btn = 
    4 * exp.N.s.yrs.2 * phi * rho / lmb.m.ph.sq * p.o.l^s.gaps +
    2 * (1 - p.o.l) * exp.N.s.yrs.2 * p.o.l^s.yrs.2 *
    ((p.o.l^(-(s.yrs.1 + 1)) - p.o.l^(-(pmin(s.yrs.1.p.a, s.yrs.2) + 1))) /
       (1 - p.o.l^(-1)) +
       ifelse(
         s.yrs.1.p.a < s.yrs.2,
         (s.yrs.2 - s.yrs.1.p.a) * p.o.l^(-s.yrs.1.p.a),
         0
       ))

  # Same-mother pairs between survey years
  exp.ns.SMPs.btn = 2 * exp.ns.SMPs.wtn[s.inds.1] * p.t.s.gs +
    s.gaps * exp.N.s.yrs.2 * beta * (1 - p.o.l) *
    lambda / lmb.m.ph.sq * p.o.l^s.gaps
  
  # Same-mother pairs between surveys, ages known (five and zero)
  exp.ns.SMPs.kwn.age.btn = exp.N.s.yrs.2 * (1 - p.o.l) * beta * 
    phi^5 * p.o.l^(s.gaps + 5)

  # Same-father pairs between survey years
  exp.ns.SFPs.btn = phi * exp.ns.SMPs.btn + 
    2 * p.t.s.gs * exp.ns.SFPs.same.b.yr.wtn[s.inds.1]

  # Full-sibling pairs between survey years, note the predicted number within
  # surveys is constant
  exp.ns.FSPs.btn = 2 * exp.ns.FSPs.wtn * p.t.s.gs +
    2 * beta * (1 - p.o.l) * rcl.prb.mtr * phi^(s.yrs.1 + s.yrs.2 + 1) * 
    (p.o.l^(s.yrs.1 + 1) - p.o.l^(s.yrs.2 + 1)) / 
    ((phi^3 / lambda)^s.yrs.1 * (1 - phi^3 / lambda) * (1 - p.o.l))

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
