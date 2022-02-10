# Function to find numbers of same-mother/father pairs in the population
# including animals born in each year in the population history
FindSMFPsT <- function(pop.cap.hist, hist.len) {
  # Parents of animals born in final year
  brn.f.yr = attributes(pop.cap.hist)$f.age == 0
  mums.of.brn.f.yr = attributes(pop.cap.hist)$mum[brn.f.yr]
  dads.of.brn.f.yr = attributes(pop.cap.hist)$dad[brn.f.yr]
  
  # Animals alive in final year
  alv.f.yr = attributes(pop.cap.hist)$alv.mat[, hist.len] == 1
  
  # Matrix for numbers of kin-pairs
  SMFPs.t.mat = matrix(NA, 7, hist.len - 2)
  
  # Loop over years between second and last in population history
  for (t in (hist.len - 2):1) {
    # Parents of animals born in current year
    brn.yr.t = attributes(pop.cap.hist)$f.age == t & alv.f.yr
    mums.of.brn.yr.t = attributes(pop.cap.hist)$mum[brn.yr.t]
    dads.of.brn.yr.t = attributes(pop.cap.hist)$dad[brn.yr.t]
    
    # Same-mother pairs between animals born in current and final years (max one
    # per mum)
    # ns.SMPs.t.f.yr[hist.len - 1 - t] = 
    SMFPs.t.mat[1, hist.len - 1 - t] = 
      sum(mums.of.brn.yr.t %in% mums.of.brn.f.yr)
    
    # Same-father pairs between animals born in current and final years (many
    # possible per dad)
    max.dad.id = max(dads.of.brn.yr.t, dads.of.brn.f.yr)
    # ns.SFPs.t.f.yr[hist.len - 1 - t] = 
    SMFPs.t.mat[2, hist.len - 1 - t] = 
      tabulate(dads.of.brn.yr.t, max.dad.id) %*% 
      tabulate(dads.of.brn.f.yr, max.dad.id)
    
    # Same-father pairs born in the current year (many possible per dad)
    # ns.SFPs.t.t[hist.len - 1 - t] = sum(choose(tabulate(dads.of.brn.yr.t), 2))
    SMFPs.t.mat[3, hist.len - 1 - t] = sum(choose(tabulate(dads.of.brn.yr.t), 2))
      
    # Same-mother pairs between each year in population history, and final year.  

    # Both born before first year - birth years are the two years before the
    # first year.
    born.tm2.alv.t = attributes(pop.cap.hist)$f.age == t + 2 & 
      attributes(pop.cap.hist)$alv.mat[, hist.len - t] == 1
    born.tm1.alv.f.yr = attributes(pop.cap.hist)$f.age == t + 1 & alv.f.yr
    born.tm2.alv.f.yr = attributes(pop.cap.hist)$f.age == t + 2 & alv.f.yr
    born.tm1.alv.t = attributes(pop.cap.hist)$f.age == t + 1 & 
      attributes(pop.cap.hist)$alv.mat[, hist.len - t] == 1
    # ns.SMPs.tm2.tm1.t.f.yr[hist.len - 1 - t] =
    SMFPs.t.mat[4, hist.len - 1 - t] = 
      sum(
        attributes(pop.cap.hist)$mum[born.tm2.alv.t] %in% 
          attributes(pop.cap.hist)$mum[born.tm1.alv.f.yr],
        attributes(pop.cap.hist)$mum[born.tm2.alv.f.yr] %in% 
          attributes(pop.cap.hist)$mum[born.tm1.alv.t]
      )
    
    # One born after first year - birth years are years before first and second
    # years.
    born.f.yrm1.alv.f.yr = attributes(pop.cap.hist)$f.age == 1 & alv.f.yr
    # ns.SMPs.tm1.f.yrm1.t.f.yr[hist.len - 1 - t] =
    SMFPs.t.mat[5, hist.len - 1 - t] = 
      sum(
        attributes(pop.cap.hist)$mum[born.tm1.alv.t] %in% 
          attributes(pop.cap.hist)$mum[born.f.yrm1.alv.f.yr]
      )
    
    # One born after first year - birth years are years before first years, and
    # all years between first and second years.  Doesn't match expression...
    born.btwn.t.f.yr.alv.f.yr = 
      attributes(pop.cap.hist)$f.age < t & alv.f.yr
    # ns.SMPs.tm1.btwn.t.f.yr.t.f.yr[hist.len - 1 - t] =
    SMFPs.t.mat[6, hist.len - 1 - t] = 
      sum(
        attributes(pop.cap.hist)$mum[born.tm1.alv.t] %in% 
          attributes(pop.cap.hist)$mum[born.btwn.t.f.yr.alv.f.yr]
      )
    
    # Both born before first year - birth years are first year in history, and
    # all years between it and first years.  Just always gives zero...
    born.fst.yr.alv.t = attributes(pop.cap.hist)$f.age == hist.len - 1 & 
      attributes(pop.cap.hist)$alv.mat[, hist.len - t] == 1
    born.aftr.fst.yr.bfr.t.alv.f.yr = 
      attributes(pop.cap.hist)$f.age < hist.len - 1 &
      attributes(pop.cap.hist)$f.age > t & alv.f.yr
    born.fst.yr.alv.f.yr = 
      attributes(pop.cap.hist)$f.age == hist.len - 1 & alv.f.yr
    born.aftr.fst.yr.bfr.t.alv.t = 
      attributes(pop.cap.hist)$f.age < hist.len - 1 &
      attributes(pop.cap.hist)$f.age > t & 
      attributes(pop.cap.hist)$alv.mat[, hist.len - t] == 1
    # ns.SMPs.fst.yr.btwn.t.f.yr[hist.len - 1 - t] =
    SMFPs.t.mat[7, hist.len - 1 - t] = 
      sum(
        attributes(pop.cap.hist)$mum[born.fst.yr.alv.t] %in% 
          attributes(pop.cap.hist)$mum[born.aftr.fst.yr.bfr.t.alv.f.yr],
        attributes(pop.cap.hist)$mum[born.fst.yr.alv.f.yr] %in% 
          attributes(pop.cap.hist)$mum[born.aftr.fst.yr.bfr.t.alv.t]
      )
  }
  
  SMFPs.t.mat
}