# Function to estimate numbers of same-mother/father pairs in the population
# including animals born in each year in the population history
FindExpNsKPsT = function(
  exp.N.fin, phi, lambda, alpha, hist.len, exp.N.t, n.yrs.chk.t
) {
  # Intermediate results
  
  # Probability that two animals are new-born independently
  prb.nw.brn.sq = (1 - phi / lambda)^2
  # Probability that an animal is mature
  prb.mtr = (lambda / phi)^alpha
  # Phi-squared over lambda
  ph.sq.ovr.lmd = phi^2 / lambda
  # Each year in history except first and last
  t.vec = (hist.len - n.yrs.chk.t):(hist.len - 1)
  # Final year minus each prior year except the first
  f.yr.m.t.vec = hist.len - t.vec
  
  # Estimates
  
  # Same-mother pairs
  exp.ns.SMPs.f.yr = 
    2 * exp.N.fin * prb.nw.brn.sq * prb.mtr * ph.sq.ovr.lmd^f.yr.m.t.vec
  
  # Same-father pairs with different ages (one always age zero)
  exp.ns.SFPs.f.yr = exp.ns.SMPs.f.yr * phi
  
  # Same-father pairs with same ages (for all ages up to length of population
  # history)
  exp.ns.SFPs.t.t = phi^(2 * f.yr.m.t.vec) * 
    (exp.N.fin / lambda^(alpha + f.yr.m.t.vec) - 1) *
    (lambda^2 / phi)^alpha * lambda * prb.nw.brn.sq
  
  # Same-mother pairs between each year and final year, born in two years
  # preceding each year
  exp.ns.SMPs.tm2.tm1.t.f.yr = 4 * exp.N.t[t.vec] * 
    prb.nw.brn.sq * prb.mtr * ph.sq.ovr.lmd^2 * phi^f.yr.m.t.vec
  
  # Same-mother pairs between each year and final year, born in years
  # preceding both
  exp.ns.SMPs.tm1.f.yrm1.t.f.yr = 2 * exp.N.t[t.vec] * 
    prb.nw.brn.sq * prb.mtr * ph.sq.ovr.lmd * phi^f.yr.m.t.vec
  
  # Same-mother pairs between each year and final year, born year before each
  # year, and all years in between each year and final year
  exp.ns.SMPs.tm1.btwn.t.f.yr.t.f.yr = 2 * exp.N.t[t.vec] * 
    prb.nw.brn.sq * prb.mtr * ph.sq.ovr.lmd * phi^f.yr.m.t.vec * f.yr.m.t.vec
  
  # Same-mother pairs between each year and final year, born in first year of
  # history, and all years in between first of history and each year
  exp.ns.SMPs.fst.yr.btwn.t.f.yr = 4 * exp.N.t[t.vec] * prb.nw.brn.sq * 
    prb.mtr * ph.sq.ovr.lmd^(t.vec - 1) * phi^f.yr.m.t.vec * (t.vec - 1)

  # Combine and return
  cbind(
    exp.ns.SMPs.f.yr, exp.ns.SFPs.f.yr, exp.ns.SFPs.t.t, 
    exp.ns.SMPs.tm2.tm1.t.f.yr, exp.ns.SMPs.tm1.f.yrm1.t.f.yr,
    exp.ns.SMPs.tm1.btwn.t.f.yr.t.f.yr, exp.ns.SMPs.fst.yr.btwn.t.f.yr
  )
}