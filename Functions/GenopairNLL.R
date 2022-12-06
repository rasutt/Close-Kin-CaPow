# Find genopair model negative log likelihood
GenopairNLLR <- function(
    params, k, f.year, srvy.yrs, alpha, gpps, smp.yr.ind.prs
) {
  # Unpack parameters
  rho <- params[1]
  phi <- params[2]
  N.fin <- params[3]
  lambda <- rho + phi
  
  # Rho must be strictly greater than zero for the closed form series below
  # to hold, otherwise nlminb gives warnings
  if (rho == 0) return(Inf)
  
  # Lambda minus phi-squared
  lmb.m.ph.sq = lambda - phi^2
  # Probability not new-born (phi over lambda)
  p.o.l = phi / lambda
  # Reciprocal of probability not new-born (phi over lambda)
  l.o.p = lambda / phi
  
  # Find the expected number alive at the survey year
  exp.N.s.yrs <- N.fin / lambda^(f.year - srvy.yrs)

  # All pairs in survey years
  exp.ns.APs.wtn = choose(exp.N.s.yrs, 2)
  
  # Parent-offspring pairs within survey years
  exp.ns.POPs.wtn = exp.N.s.yrs * rho * (1 + phi) / lmb.m.ph.sq
  
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
  exp.ns.POPs.btn = 2 * exp.ns.POPs.wtn[s.inds.1] * p.t.s.gs +
    2 * exp.N.s.yrs.2 * (1 - p.o.l) * p.o.l^s.yrs.2 *
    ((l.o.p^(s.yrs.1 + 1) - l.o.p^(pmin(s.yrs.1.p.a, s.yrs.2) + 1)) /
       (1 - l.o.p) +
       ifelse(
         s.yrs.1.p.a < s.yrs.2,
         (s.yrs.2 - s.yrs.1.p.a) * l.o.p^s.yrs.1.p.a,
         0
       ))

  # Create matrices for probabilities of POPs and SPs over survey-pairs
  prb.POPs.mat <- prb.SPs.mat <- matrix(0, k, k)
  diag(prb.POPs.mat) = exp.ns.POPs.wtn / exp.ns.APs.wtn
  prb.POPs.mat[upper.tri(prb.POPs.mat)] = exp.ns.POPs.btn / exp.ns.APs.btn
  prb.SPs.mat[upper.tri(prb.SPs.mat)] = exp.ns.SPs.btn / exp.ns.APs.btn
  prb.UPs.mat = 1 - prb.POPs.mat - prb.SPs.mat
  
  # Find negative log-likelihood from genopair probabilities given kinships,
  # n_pairs x n_kinships (UP, POP, SP), survey-year index-pairs, n_pairs x 2,
  # and kinship probabilities given parameters, k x k for each kinship
  nll = -sum(log(
    gpps[, 1] * prb.UPs.mat[smp.yr.ind.prs] +
      gpps[, 2] * prb.POPs.mat[smp.yr.ind.prs] +
      gpps[, 3] * prb.SPs.mat[smp.yr.ind.prs]
  ))
  print(nll, digits = 20)
  
  # Return negative log likelihood
  nll
}