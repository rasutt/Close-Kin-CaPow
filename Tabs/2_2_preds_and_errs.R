## Predictions and proportional errors

# Predicted numbers of kin-pairs for whole population
pred.ns.kps.pop = reactive({
  FindPredNsKPsPop(
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
  reactive(list(pred.ns.kps.pop()$wtn[, id], pred.ns.kps.pop()$btn[, id]))
}

# Predictions
ns.SPs.preds = reactive(rep(list(pred.ns.kps.pop()$btn[, "SPs"]), 2))
preds.lst = reactive({
  lst = lapply(rglr.kp.ids, function(kp.id) {
    list(pred.ns.kps.pop()$wtn[, kp.id], pred.ns.kps.pop()$btn[, kp.id])
  })
  names(lst) = rglr.kp.ids
  lst
})
ns.SibPs.preds = reactive(c(preds.lst()[["FSPs"]], preds.lst()[["HSPs"]]))

# Errors, can't be done with lists as want to keep reactives separate
N.errs = reactive(list(find.errs(N.s.yrs(), pred.ns.kps.pop()$wtn[, 1])))
phi.errs = reactive({
  list(matrix(
    find.errs(avg.phi.obs(), phi()), n.sims(), 
    dimnames = list(NULL, "All times")
  ))
})
ns.APs.errs = reactive(l.fnd.errs(ns.APs(), preds.lst()[["APs"]]))
ns.SPs.errs = reactive(l.fnd.errs(ns.SPs(), ns.SPs.preds()))
ns.POPs.errs = reactive(l.fnd.errs(ns.POPs(), preds.lst()[["POPs"]]))
ns.SMPs.errs = reactive(l.fnd.errs(ns.SMPs(), preds.lst()[["SMPs"]]))
ns.SFPs.errs = reactive(l.fnd.errs(ns.SFPs(), preds.lst()[["SFPs"]]))
ns.FSPs.errs = reactive(l.fnd.errs(ns.FSPs(), preds.lst()[["FSPs"]]))
ns.HSPs.errs = reactive(l.fnd.errs(ns.HSPs(), preds.lst()[["HSPs"]]))
