# Function to find numbers of same-mother/father pairs in the population
# including animals born in each year in the population history
FindNsKPsT <- function(pop.cap.hist, hist.len, n.yrs.chk.t) {
  # Parents of animals born in final year
  brn.f.yr = attributes(pop.cap.hist)$f.age == 0
  mums.of.brn.f.yr = attributes(pop.cap.hist)$mum[brn.f.yr]
  dads.of.brn.f.yr = attributes(pop.cap.hist)$dad[brn.f.yr]
  
  # Animals alive in final year
  alv.f.yr = attributes(pop.cap.hist)$alive == 1
  
  # Vectors for numbers of kin-pairs
  ns.SMPs.fnl.b1.fnl = ns.SFPs.fnl.b1.fnl = ns.SFPs.fnl.b = integer(n.yrs.chk.t)
  
  # Loop over check-years from earliest to latest
  for (t in 1:n.yrs.chk.t) {
    # Parents of animals born in current year
    brn.yr.t = attributes(pop.cap.hist)$f.age == t & alv.f.yr
    mums.of.brn.yr.t = attributes(pop.cap.hist)$mum[brn.yr.t]
    dads.of.brn.yr.t = attributes(pop.cap.hist)$dad[brn.yr.t]
    chk.yr.ind = n.yrs.chk.t + 1 - t
    
    # Pairs in the final year, with one born in the current year, and one born
    # in the final year
    
    # Same-mother pairs (max one per mum)
    ns.SMPs.fnl.b1.fnl[chk.yr.ind] = sum(mums.of.brn.yr.t %in% mums.of.brn.f.yr)
    
    # Same-father pairs (many possible per dad)
    max.dad.id = max(dads.of.brn.yr.t, dads.of.brn.f.yr)
    ns.SFPs.fnl.b1.fnl[chk.yr.ind] = 
      tabulate(dads.of.brn.yr.t, max.dad.id) %*% 
      tabulate(dads.of.brn.f.yr, max.dad.id)

    # Same-father pairs in the final year, both born in the current year (many
    # possible per dad)
    ns.SFPs.fnl.b[chk.yr.ind] = sum(choose(tabulate(dads.of.brn.yr.t), 2))
    
    # # Same-mother pairs between each year in population history, and final year.
    # 
    # # Both born before first year - birth years are the two years before the
    # # first year.
    # born.tm2.alv.t = attributes(pop.cap.hist)$f.age == t + 2 & 
    #   attributes(pop.cap.hist)$alv.mat[, hist.len - t] == 1
    # born.tm1.alv.f.yr = attributes(pop.cap.hist)$f.age == t + 1 & alv.f.yr
    # born.tm2.alv.f.yr = attributes(pop.cap.hist)$f.age == t + 2 & alv.f.yr
    # born.tm1.alv.t = attributes(pop.cap.hist)$f.age == t + 1 & 
    #   attributes(pop.cap.hist)$alv.mat[, hist.len - t] == 1
    # ns.kps.t.mat[4, hist.len - 1 - t] = 
    #   sum(
    #     attributes(pop.cap.hist)$mum[born.tm2.alv.t] %in% 
    #       attributes(pop.cap.hist)$mum[born.tm1.alv.f.yr],
    #     attributes(pop.cap.hist)$mum[born.tm2.alv.f.yr] %in% 
    #       attributes(pop.cap.hist)$mum[born.tm1.alv.t]
    #   )
    # 
    # # One born after first year - birth years are years before first and second
    # # years.
    # born.f.yrm1.alv.f.yr = attributes(pop.cap.hist)$f.age == 1 & alv.f.yr
    # ns.kps.t.mat[5, hist.len - 1 - t] = 
    #   sum(
    #     attributes(pop.cap.hist)$mum[born.tm1.alv.t] %in% 
    #       attributes(pop.cap.hist)$mum[born.f.yrm1.alv.f.yr]
    #   )
    # 
    # # One born after first year - birth years are years before first years, and
    # # all years between first and second years.  Doesn't match expression...
    # born.btwn.t.f.yr.alv.f.yr = 
    #   attributes(pop.cap.hist)$f.age < t & alv.f.yr
    # ns.kps.t.mat[6, hist.len - 1 - t] = 
    #   sum(
    #     attributes(pop.cap.hist)$mum[born.tm1.alv.t] %in% 
    #       attributes(pop.cap.hist)$mum[born.btwn.t.f.yr.alv.f.yr]
    #   )
    # 
    # # Both born before first year - birth years are first year in history, and
    # # all years between it and first years.  Just always gives zero...
    # born.fst.yr.alv.t = attributes(pop.cap.hist)$f.age == hist.len - 1 & 
    #   attributes(pop.cap.hist)$alv.mat[, hist.len - t] == 1
    # born.aftr.fst.yr.bfr.t.alv.f.yr = 
    #   attributes(pop.cap.hist)$f.age < hist.len - 1 &
    #   attributes(pop.cap.hist)$f.age > t & alv.f.yr
    # born.fst.yr.alv.f.yr = 
    #   attributes(pop.cap.hist)$f.age == hist.len - 1 & alv.f.yr
    # born.aftr.fst.yr.bfr.t.alv.t = 
    #   attributes(pop.cap.hist)$f.age < hist.len - 1 &
    #   attributes(pop.cap.hist)$f.age > t & 
    #   attributes(pop.cap.hist)$alv.mat[, hist.len - t] == 1
    # ns.kps.t.mat[7, hist.len - 1 - t] = 
    #   sum(
    #     attributes(pop.cap.hist)$mum[born.fst.yr.alv.t] %in% 
    #       attributes(pop.cap.hist)$mum[born.aftr.fst.yr.bfr.t.alv.f.yr],
    #     attributes(pop.cap.hist)$mum[born.fst.yr.alv.f.yr] %in% 
    #       attributes(pop.cap.hist)$mum[born.aftr.fst.yr.bfr.t.alv.t]
    #   )
  }
  
  cbind(ns.SMPs.fnl.b1.fnl, ns.SFPs.fnl.b1.fnl, ns.SFPs.fnl.b)
}