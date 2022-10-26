## Predictions and proportional errors

# Estimated numbers of kin-pairs for whole population
est.ns.kps.pop = reactive({
  FindEstNsKPsPop(
    exp.N.t(), s.yr.inds(), phi(), rho(), lambda(), alpha(), srvy.yrs(), k()
  )
})

# Temporally estimated numbers of kin-pairs
est.ns.kps.t = reactive({
  FindEstNsKPsT(
    exp.N.fin(), phi(), lambda(), alpha(), hist.len(), exp.N.t(), n.yrs.chk.t()
  )
})

# Function to find estimate errors as proportions of estimates
find.errs = function(vals, ests, samp = F) {
  # For sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) ests = rep(ests, each = n.sims())
  vals / ests - 1
}

# Function to find list of estimate errors
l.fnd.errs = function(ns, preds) {
  lapply(1:length(ns), function(i) find.errs(ns[[i]], preds[[i]]))
}

# Function to combine predictions within and between survey-years
comb.preds = function(id) {
  reactive(list(est.ns.kps.pop()$wtn[, id], est.ns.kps.pop()$btn[, id]))
}

# Predictions, can't be done with lists as want to keep reactives separate
ns.SPs.preds = reactive(rep(list(est.ns.kps.pop()$btn[, "SPs"]), 2))
preds.lst = reactive({
  lst = lapply(rglr.kp.ids, function(kp.id) {
    list(est.ns.kps.pop()$wtn[, kp.id], est.ns.kps.pop()$btn[, kp.id])
  })
  names(lst) = rglr.kp.ids
  lst
})
# ns.APs.preds = comb.preds("APs")
ns.APs.preds = reactive(preds.lst()[["APs"]])
ns.POPs.preds = comb.preds("POPs")
ns.SMPs.preds = comb.preds("SMPs")
ns.SFPs.preds = comb.preds("SFPs")
ns.FSPs.preds = comb.preds("FSPs")
ns.HSPs.preds = comb.preds("HSPs")
ns.SibPs.preds = reactive(c(ns.FSPs.preds, ns.HSPs.preds))

# Errors, can't be done with lists as want to keep reactives separate
N.errs = reactive(list(find.errs(N.s.yrs(), est.ns.kps.pop()$wtn[, 1])))
phi.errs = reactive({
  list(matrix(
    find.errs(avg.phi.obs(), phi()), n.sims(), 
    dimnames = list(NULL, "All times")
  ))
})
ns.APs.errs = reactive(l.fnd.errs(ns.APs(), preds.lst()[["APs"]]))
ns.SPs.errs = reactive(l.fnd.errs(ns.SPs(), ns.SPs.preds()))
ns.POPs.errs = reactive(l.fnd.errs(ns.POPs(), ns.POPs.preds()))
ns.SMPs.errs = reactive(l.fnd.errs(ns.SMPs(), ns.SMPs.preds()))
ns.SFPs.errs = reactive(l.fnd.errs(ns.SFPs(), ns.SFPs.preds()))
ns.FSPs.errs = reactive(l.fnd.errs(ns.FSPs(), ns.FSPs.preds()))
ns.HSPs.errs = reactive(l.fnd.errs(ns.HSPs(), ns.HSPs.preds()))
