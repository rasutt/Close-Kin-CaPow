# Separate functions to find numbers of kin-pairs

# Find self-pairs between survey years
find.SPs = function(pop.atts, k) {
  as.vector(combn(1:k, 2, function(s.inds) {
    sum(pop.atts$ID[pop.atts$alv.s.yrs[, s.inds[1]]] %in% 
          pop.atts$ID[pop.atts$alv.s.yrs[, s.inds[2]]])
  }))
}

# Find proportions with unknown parents in survey-years
find.pns.UPs.wtn = function(pop.atts, k) {
  sapply(1:k, function(s.ind){
    mean(is.na(pop.atts$mum[pop.atts$alv.s.yrs[, s.ind]]))
  })
}

# Find proportions with unknown parents in pairs of survey-years
find.pns.UPs.btn = function(pop.atts, k) {
  as.vector(combn(1:k, 2, function(s.inds) {
    mean(c(
      is.na(pop.atts$mum[pop.atts$alv.s.yrs[, s.inds[1]]]),
      is.na(pop.atts$mum[pop.atts$alv.s.yrs[, s.inds[2]]])
    ))
  }))
}

# Find parent-offspring pairs in survey-years
find.POPs.wtn = function(pop.atts, k) {
  alv.s.yrs = pop.atts$alv.s.yrs
  sapply(1:k, function(s.ind) {
    sum(
      pop.atts$mum[alv.s.yrs[, s.ind]] %in% pop.atts$ID[alv.s.yrs[, s.ind]], 
      pop.atts$dad[alv.s.yrs[, s.ind]] %in% pop.atts$ID[alv.s.yrs[, s.ind]]
    )
  })
}

# Find parent-offspring pairs between pairs of survey-years
find.POPs.btn = function(pop.atts, k) {
  alv.s.yrs = pop.atts$alv.s.yrs
  ID = pop.atts$ID
  as.vector(combn(1:k, 2, function(s.inds) {
    sum(
      pop.atts$mum[alv.s.yrs[, s.inds[1]]] %in% ID[alv.s.yrs[, s.inds[2]]], 
      pop.atts$dad[alv.s.yrs[, s.inds[1]]] %in% ID[alv.s.yrs[, s.inds[2]]],
      pop.atts$mum[alv.s.yrs[, s.inds[2]]] %in% ID[alv.s.yrs[, s.inds[1]]], 
      pop.atts$dad[alv.s.yrs[, s.inds[2]]] %in% ID[alv.s.yrs[, s.inds[1]]]
    )
  }))
}

# Find same-mother pairs with ages unknown, in survey-years
find.SMPs.wtn = function(pop.atts, k) {
  # List of frequency tables of mums in each survey year
  max.mum = max(pop.atts$mum, na.rm = T)
  mum.tab.lst = lapply(1:k, function(s.ind) {
    tabulate(pop.atts$mum[pop.atts$alv.s.yrs[, s.ind]], max.mum)
  })
  
  # Same-mother pairs within survey years
  sapply(1:k, function(s.ind) sum(choose(mum.tab.lst[[s.ind]], 2)))
}

# Find same-mother pairs with ages unknown, and self-pairs, between survey-years
find.SMSPs.btn = function(pop.atts, k) {
  # List of frequency tables of mums in each survey year
  max.mum = max(pop.atts$mum, na.rm = T)
  mum.tab.lst = lapply(1:k, function(s.ind) {
    tabulate(pop.atts$mum[pop.atts$alv.s.yrs[, s.ind]], max.mum)
  })
  
  # Same-mother and self-pairs between survey years
  as.vector(combn(1:k, 2, function(s.inds) {
    mum.tab.lst[[s.inds[1]]] %*% mum.tab.lst[[s.inds[2]]]
  }))
}

# Find same-mother pairs with ages known (five and zero), between survey-years
find.SMPs.age.known.btn = function(pop.atts, k, fnl.year, srvy.yrs) {
  alv.s.yrs = pop.atts$alv.s.yrs
  new.born = sapply(1:k, function(s.ind) {
    pop.atts$f.age == fnl.year - srvy.yrs[s.ind]
  })
  fv.yrs.old = sapply(1:k, function(s.ind) {
    pop.atts$f.age == fnl.year - srvy.yrs[s.ind] + 5
  })
  
  # List of frequency tables of mums of animals born in each survey year
  max.mum = max(pop.atts$mum, na.rm = T)
  mum.tab.lst.nb = lapply(1:k, function(s.ind) {
    tabulate(pop.atts$mum[alv.s.yrs[, s.ind] & new.born[, s.ind]], max.mum)
  })
  mum.tab.lst.fyo = lapply(1:k, function(s.ind) {
    tabulate(pop.atts$mum[alv.s.yrs[, s.ind] & fv.yrs.old[, s.ind]], max.mum)
  })
  
  # Same-mother pairs between survey years - No self-pairs as different birth
  # years
  as.vector(combn(1:k, 2, function(s.inds) {
    mum.tab.lst.fyo[[s.inds[1]]] %*% mum.tab.lst.nb[[s.inds[2]]]
  }))
}

# Find same-father pairs with ages unknown, in survey-years
find.SFPs.wtn = function(pop.atts, k) {
  # List of frequency tables of dads in each survey year
  max.dad = max(pop.atts$dad, na.rm = T)
  dad.tab.lst = lapply(1:k, function(s.ind) {
    tabulate(pop.atts$dad[pop.atts$alv.s.yrs[, s.ind]], max.dad)
  })
  
  # Same-father pairs within survey years
  sapply(1:k, function(s.ind) sum(choose(dad.tab.lst[[s.ind]], 2)))
}

# Find same-father pairs with ages unknown, and self-pairs, between survey-years
find.SFSPs.btn = function(pop.atts, k) {
  # List of frequency tables of dads in each survey year
  max.dad = max(pop.atts$dad, na.rm = T)
  dad.tab.lst = lapply(1:k, function(s.ind) {
    tabulate(pop.atts$dad[pop.atts$alv.s.yrs[, s.ind]], max.dad)
  })
  
  # Same-father and self-pairs between survey years
  as.vector(combn(1:k, 2, function(s.inds) {
    dad.tab.lst[[s.inds[1]]] %*% dad.tab.lst[[s.inds[2]]]
  }))
}

# Find full-sibling pairs with ages unknown, in survey-years
find.FSPs.wtn = function(pop.atts, k) {
  # Full-sibling pairs within survey years
  sapply(1:k, function(s.ind) {
    tbl = table(
      pop.atts$mum[pop.atts$alv.s.yrs[, s.ind]], 
      pop.atts$dad[pop.atts$alv.s.yrs[, s.ind]]
    )
    
    # Selecting before calling choose is 10x faster
    sum(choose(tbl[tbl > 1], 2))
  })
}

# Self-pairs with known parents, to subtract from sibling pairs found through
# parents
find.SPs.prnts.kwn = function(pop.atts, k) {
  # Indices of animals alive each survey, with known parents
  alv.KwnPs.lst = lapply(1:k, function(s.ind) {
    which(pop.atts$alv.s.yrs[, s.ind] & !is.na(pop.atts$mum))
  })
  as.vector(combn(1:k, 2, function(s.inds) {
    sum(alv.KwnPs.lst[[s.inds[1]]] %in% alv.KwnPs.lst[[s.inds[2]]])
  }))
}

# Find full-sibling pairs and self-pairs, with ages unknown, between survey-years
find.FSSPs.btn = function(pop.atts, k) {
  # Indices of animals alive each survey, with known parents
  alv.KwnPs.lst = lapply(1:k, function(s.ind) {
    which(pop.atts$alv.s.yrs[, s.ind] & !is.na(pop.atts$mum))
  })
  
  # List of frequency tables of mums and dads in each survey year
  mum.lst = lapply(alv.KwnPs.lst, function(alv.KwnPs) pop.atts$mum[alv.KwnPs])
  dad.lst = lapply(alv.KwnPs.lst, function(alv.KwnPs) pop.atts$dad[alv.KwnPs])
  mums = unique(unlist(mum.lst))
  dads = unique(unlist(dad.lst))
  
  # Two-way frequency tables of mothers and fathers of animals in each survey
  # year, including all possible parents in levels so that dimensions match
  tbl.lst = lapply(1:k, function(s.ind) {
    table(
      factor(mum.lst[[s.ind]], levels = mums), 
      factor(dad.lst[[s.ind]], levels = dads)
    )
  })
  
  # Full-sibling pairs, including self-pairs with known parents
  as.vector(combn(1:k, 2, function(s.inds) {
    sum(tbl.lst[[s.inds[1]]] * tbl.lst[[s.inds[2]]])
  }))
}

# Same-mother pairs with ages known, in final year
find.SMPs.age.knwn = function(pop.atts, n.yrs.chk.t) {
  sapply(n.yrs.chk.t:1, function(t) {
    # Same-mother pairs in the final year, with one born t years before the
    # final year, and one born in the final year (max one per mum)
    sum(
      pop.atts$mum[pop.atts$f.age == t & pop.atts$alive == 1] %in% 
        pop.atts$mum[pop.atts$f.age == 0]
    )
  })
}

# Same-father pairs with ages known, in final year
find.SFPs.age.knwn = function(pop.atts, n.yrs.chk.t) {
  # Fathers of animals born in final year
  dads.of.brn.f.yr = pop.atts$dad[pop.atts$f.age == 0]
  
  # Same-father pairs in the final year, with one born t years before the final
  # year, and one born in the final year (many possible per dad)
  sapply(n.yrs.chk.t:1, function(t) {
    dads.of.brn.yr.t = pop.atts$dad[pop.atts$f.age == t & pop.atts$alive == 1]
    max.dad.id = max(dads.of.brn.yr.t, dads.of.brn.f.yr)
    tabulate(dads.of.brn.yr.t, max.dad.id) %*% 
      tabulate(dads.of.brn.f.yr, max.dad.id)
  })
}

# Same-father pairs with same ages known, in final year
find.SFPs.same.age = function(pop.atts, n.yrs.chk.t) {
  # Same-father pairs in the final year, both born in the current year
  # (many possible per dad)
  sapply(n.yrs.chk.t:1, function(t) {
    sum(choose(
      tabulate(pop.atts$dad[pop.atts$f.age == t & pop.atts$alive == 1]), 2
    ))
  })
  
}

