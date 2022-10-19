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

# Population sizes ----
ns.wtn.errs = reactive(list(find.errs(N.s.yrs(), est.ns.kps.pop()$wtn[, 1])))
# Survival rates ----
phi.errs = reactive({
  list(matrix(
    find.errs(avg.phi.obs(), phi()), n.sims(), 
    dimnames = list(NULL, "All times")
  ))
})
# All-pairs ----
ns.APs.preds = reactive({
  list(est.ns.kps.pop()$wtn[, "APs"], est.ns.kps.pop()$btn[, "APs"])
})
ns.APs.errs = reactive(l.fnd.errs(ns.APs(), ns.APs.preds()))
# Self-pairs ----
ns.SPs.preds = reactive(rep(list(est.ns.kps.pop()$btn[, "SPs"]), 2))
ns.SPs.errs = reactive(l.fnd.errs(ns.SPs(), ns.SPs.preds()))
# Parent-offspring pairs ----
ns.POPs.preds = reactive(list(
  est.ns.kps.pop()$wtn[, "POPs"], est.ns.kps.pop()$btn[, "POPs"]
))
ns.POPs.errs = reactive(l.fnd.errs(ns.POPs(), ns.POPs.preds()))
# Same-mother pairs ----
ns.SMPs.preds = reactive(list(
  est.ns.kps.t()[, "SMPs.kwn.age"],  est.ns.kps.pop()$wtn[, "SMPs"], 
  est.ns.kps.pop()$btn[, "SMPs"], est.ns.kps.pop()$btn[, "SMPs.kwn.age"]
))
ns.SMPs.errs = reactive(l.fnd.errs(ns.SMPs.t(), ns.SMPs.preds()))
# Same-father pairs ----
ns.SFPs.preds = reactive(list(
  est.ns.kps.t()[, "SFPs.kwn.age"],  est.ns.kps.t()[, "SFPs.sm.age"], 
  est.ns.kps.pop()$wtn[, "SFPs"], est.ns.kps.pop()$btn[, "SFPs"]
))
ns.SFPs.errs = reactive(l.fnd.errs(ns.SFPs.t(), ns.SFPs.preds()))
# Sibling-pairs ----
ns.FSPs.preds = reactive(list(
  est.ns.kps.pop()$wtn[, "FSPs"], est.ns.kps.pop()$btn[, "FSPs"]
))
ns.HSPs.preds = reactive(list(
  est.ns.kps.pop()$wtn[, "HSPs"], est.ns.kps.pop()$btn[, "HSPs"]
))
ns.SibPs.preds = reactive(c(ns.FSPs.preds, ns.HSPs.preds))
ns.FSPs.errs = reactive(l.fnd.errs(ns.FSPs(), ns.FSPs.preds()))
ns.HSPs.errs = reactive(l.fnd.errs(ns.HSPs(), ns.HSPs.preds()))
