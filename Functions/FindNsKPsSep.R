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