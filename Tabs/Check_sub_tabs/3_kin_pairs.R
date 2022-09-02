# Outputs for kin-pairs sub-tab of checks tab

# Show percentage of animals captured for which the parents are unknown
output$percUnknPrnts = renderText({
  paste0(
    "Percentage of captured animals with unknown parents: ", 
    round(mean(checks.lst()$prpn.prnts.unkn.vec) * 100, 1), "%"
  )
})

# Estimated numbers of kin-pairs for whole population
est.ns.kps.pop.lst = reactive({
  FindEstNsKPsPop(
    exp.N.t(), s.yr.inds(), phi(), lambda(), alpha(), srvy.yrs(), k()
  )
})

### Kin-pair estimator biases (tables of average percentage differences)

# Function to find estimate errors as proportions of estimates
find.est.errs = function(vals, ests, samp = F) {
  # For sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) ests = rep(ests, each = n.sims())
  vals / ests - 1
}

# Errors in estimates for numbers of kin-pairs in whole population within
# surveys
ns.kps.pop.wtn.est.errs = reactive({
  find.est.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn)
})
ns.kps.pop.btn.est.errs = reactive({
  find.est.errs(checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn)
})
ns.kps.prb.wtn.est.errs = reactive({
  # Remove population sizes and total numbers of pairs then divide by the latter
  find.est.errs(
    checks.lst()$ns.kps.pop.wtn.arr[, , -1:-2] / 
      array(
        rep(checks.lst()$ns.kps.pop.wtn.arr[, , 2], n.kp.tps.prb.wtn), 
        c(n.sims(), k(), n.kp.tps.prb.wtn)
      ), 
    est.ns.kps.pop.lst()$wtn[, -1:-2] / est.ns.kps.pop.lst()$wtn[, 2]
  )
})
ns.kps.prb.btn.est.errs = reactive({
  # Remove population sizes and total numbers of pairs then divide by the latter
  find.est.errs(
    checks.lst()$ns.kps.pop.btn.arr[, , -1] / 
      array(
        rep(checks.lst()$ns.kps.pop.btn.arr[, , 1], n.kp.tps.prb.btn), 
        c(n.sims(), n.srvy.prs(), n.kp.tps.prb.btn)
      ), 
    est.ns.kps.pop.lst()$btn[, -1] / est.ns.kps.pop.lst()$btn[, 1]
  )
})
ns.kps.cap.wtn.est.errs = reactive({
  find.est.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn, T)
})
ns.kps.cap.btn.est.errs = reactive({
  find.est.errs(checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn, T)
})

## Function to find them from true and expected values
find.est.bias = function(est.errs, type_names) {
  df = data.frame(matrix(perc(colMeans(est.errs, dims = 2)), nrow = 1))
  names(df) = type_names
  df
}
est.bias = function(vals, exp.vals, type_names, samp = F) {
  # For the sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) exp.vals = rep(exp.vals, each = n.sims())
  df = data.frame(
    matrix(perc(colMeans(vals / exp.vals, dims = 2) - 1), nrow = 1)
  )
  names(df) = type_names
  df
}

## Numbers in whole population

# Within surveys
output$biasNsKPsPopWtn = renderTable({
  find.est.bias(ns.kps.pop.wtn.est.errs(), kp.tps.pop.wtn)
})

# Between surveys
output$biasNsKPsPopBtn = renderTable({
  find.est.bias(ns.kps.pop.btn.est.errs(), kp.tps.pop.btn)
})

## Probabilities (numbers divided by total numbers of pairs)

# Within surveys
output$biasProbsKPsWtn = renderTable({
  find.est.bias(ns.kps.prb.wtn.est.errs(), kp.tps.prb.wtn)
})

# Between surveys
output$biasProbsKPsBtn = renderTable({
  find.est.bias(ns.kps.prb.btn.est.errs(), kp.tps.prb.btn)
})

## Numbers among sampled animals

# Within surveys
output$biasNsKPsCapWtn = renderTable({
  find.est.bias(ns.kps.cap.wtn.est.errs(), kp.tps.cap.wtn)
})

# Between surveys
output$biasNsKPsCapBtn = renderTable({
  find.est.bias(ns.kps.cap.btn.est.errs(), kp.tps.cap.btn)
})
