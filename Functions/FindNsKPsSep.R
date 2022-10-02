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

# Find same-mother pairs with ages unknown, between survey-years
find.SMPs.btn = function(pop.atts, k) {
  # List of frequency tables of mums in each survey year
  max.mum = max(pop.atts$mum, na.rm = T)
  mum.tab.lst = lapply(1:k, function(s.ind) {
    tabulate(pop.atts$mum[pop.atts$alv.s.yrs[, s.ind]], max.mum)
  })
  
  # Same-mother pairs between survey years
  as.vector(combn(1:k, 2, function(s.inds) {
    mum.tab.lst[[s.inds[1]]] %*% mum.tab.lst[[s.inds[2]]]
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

# Find same-father pairs with ages unknown, between survey-years
find.SFPs.btn = function(pop.atts, k) {
  # List of frequency tables of dads in each survey year
  max.dad = max(pop.atts$dad, na.rm = T)
  dad.tab.lst = lapply(1:k, function(s.ind) {
    tabulate(pop.atts$dad[pop.atts$alv.s.yrs[, s.ind]], max.dad)
  })
  
  # Same-father pairs between survey years
  as.vector(combn(1:k, 2, function(s.inds) {
    dad.tab.lst[[s.inds[1]]] %*% dad.tab.lst[[s.inds[2]]]
  }))
}

# Find full-sibling pairs with ages unknown, in survey-years
find.FSPs.wtn = function(pop.atts, k) {
  # Full-sibling pairs within survey years
  sapply(1:k, function(s.ind) {
    sum(choose(table(
        pop.atts$mum[pop.atts$alv.s.yrs[, s.ind]], 
        pop.atts$dad[pop.atts$alv.s.yrs[, s.ind]]
      ), 2))
  })
}

