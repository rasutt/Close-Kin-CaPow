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
  if (!samp) ests = rep(ests, each = n_sims())
  vals / ests - 1
}

# Errors in estimates for numbers of kin-pairs in whole population within
# surveys
ns.kps.pop.wtn.est.errs = reactive({
  find.est.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn)
})

## Function to find them from true and expected values
est.bias = function(vals, exp.vals, type_names, samp = F) {
  # For the sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) exp.vals = rep(exp.vals, each = n_sims())
  df = data.frame(
    matrix(perc(colMeans(vals / exp.vals, dims = 2) - 1), nrow = 1)
  )
  names(df) = type_names
  df
}

## Numbers in whole population

# Within surveys
output$biasNsKPsPopWtn = renderTable({
  est.bias(
    checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn, kp.tps.pop.wtn
  )
})

# Between surveys
output$biasNsKPsPopBtn = renderTable({
  est.bias(
    checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn, kp.tps.pop.btn
  )
})

## Probabilities (numbers divided by total numbers of pairs)

# Within surveys
output$biasProbsKPsWtn = renderTable({
  # Remove population sizes and total numbers of pairs then divide by the latter
  est.bias(
    checks.lst()$ns.kps.pop.wtn.arr[, , -1:-2] / 
      array(
        rep(checks.lst()$ns.kps.pop.wtn.arr[, , 2], n.kp.tps.prb.wtn), 
        c(n_sims(), k(), n.kp.tps.prb.wtn)
      ), 
    est.ns.kps.pop.lst()$wtn[, -1:-2] / est.ns.kps.pop.lst()$wtn[, 2], 
    kp.tps.prbs.wtn
  )
})

# Between surveys
output$biasProbsKPsBtn = renderTable({
  est.bias(
    # Remove total numbers of pairs then divide by them
    checks.lst()$ns.kps.pop.btn.arr[, , -1] / 
      array(
        rep(checks.lst()$ns.kps.pop.btn.arr[, , 1], n.kp.tps.prb.btn), 
        c(n_sims(), n.srvy.prs(), n.kp.tps.prb.btn)
      ), 
    est.ns.kps.pop.lst()$btn[, -1] / est.ns.kps.pop.lst()$btn[, 1], 
    kp.tps.prbs.btn
  )
})

## Numbers among sampled animals

# Within surveys
output$biasNsKPsCapWtn = renderTable({
  est.bias(
    checks.lst()$ns.kps.cap.wtn.arr,
    checks.lst()$exp.ns.kps.cap.wtn.arr, kp.tps.cap.wtn,
    cap = T
  )
})

# Between surveys
output$biasNsKPsCapBtn = renderTable({
  est.bias(
    checks.lst()$ns.kps.cap.btn.arr,
    checks.lst()$exp.ns.kps.cap.btn.arr, kp.tps.cap.btn,
    cap = T
  )
})
