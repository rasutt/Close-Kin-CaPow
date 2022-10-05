# Outputs for biases sub-tab of checks tab

# Show percentage of animals sampled for which the parents are unknown
output$percUnknPrnts = renderTable({
  df = data.frame(matrix(c(
    perc(mean(pns.UPs()$pns.UPs.wtn)),
    perc(mean(pns.UPs()$pns.UPs.btn))
  ), nrow = 1))
  names(df) = c("Survey-years", "Survey-pairs")
  df
})

# Combine/find errors for estimates in sets
ns.kps.pop.wtn.errs = reactive({
  arr = array(
    c(
      ns.wtn.errs(), ns.APs.errs()[[1]], ns.POPs.errs()[[1]], 
      ns.SMPs.errs()[[2]], ns.SFPs.wtn.errs(), ns.FSPs.wtn.errs(), 
      ns.HSPs.wtn.errs() 
    ), 
    dim = c(n.sims(), k(), n.kp.tps.pop.wtn),
    dimnames = list(NULL, Survey = srvy.yrs(), kp.type = kp.tps.pop.wtn)
  )
})
ns.kps.pop.btn.errs = reactive({
  arr = array(
    c(ns.APs.errs()[[2]], ns.SPs.errs(), ns.POPs.errs()[[2]], 
      ns.SMPs.errs()[[3]]),
    dim = c(n.sims(), n.srvy.prs(), n.kp.tps.pop.btn),
    dimnames = list(NULL, Survey_pair = srvy.prs(), kp.type = kp.tps.pop.btn)
  )
})
ns.kps.t.errs = reactive({
  arr = array(
    c(
      ns.SMPs.errs()[[1]], ns.SFPs.age.knwn.errs(), ns.SFPs.same.age.errs()
    ),
    dim = c(n.sims(), n.yrs.chk.t(), n.kp.tps.t),
    dimnames = list(NULL, Year = yrs.chk.t(), kp.type = kp.tps.t)
  )
})
ns.kps.prb.wtn.errs = reactive({
  # Remove population sizes and total numbers of pairs then divide by the latter
  find.errs(
    checks.lst()$ns.kps.pop.wtn.arr[, , -1:-2] / 
      array(
        rep(checks.lst()$ns.kps.pop.wtn.arr[, , 2], n.kp.tps.prb.wtn), 
        c(n.sims(), k(), n.kp.tps.prb.wtn)
      ), 
    est.ns.kps.pop.lst()$wtn[, -1:-2] / est.ns.kps.pop.lst()$wtn[, 2]
  )
})
ns.kps.prb.btn.errs = reactive({
  # Remove population sizes and total numbers of pairs then divide by the latter
  # (drop = F retains array dimensions when only one survey pair)
  find.errs(
    checks.lst()$ns.kps.pop.btn.arr[, , -1, drop = F] / 
      array(
        rep(checks.lst()$ns.kps.pop.btn.arr[, , 1], n.kp.tps.prb.btn), 
        c(n.sims(), n.srvy.prs(), n.kp.tps.prb.btn)
      ), 
    est.ns.kps.pop.lst()$btn[, -1] / est.ns.kps.pop.lst()$btn[, 1]
  )
})
# ns.kps.cap.wtn.errs = reactive({
#   find.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn, T)
# })
# ns.kps.cap.btn.errs = reactive({
#   find.errs(checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn, T)
# })

### Estimator biases (tables of average percentage differences)

## Numbers in whole population

# Temporal estimates
output$tempEstBiasNote = renderText(
  "These outputs seem to be higher variance, but 1000 studies seems to be 
  enough."
)
output$biasNsKPsTemp = renderTable(find.bias(ns.kps.t.errs()))

# Within surveys
output$biasNsKPsPopWtn = renderTable(find.bias(ns.kps.pop.wtn.errs()))

# Between surveys
output$biasNsKPsPopBtn = renderTable(find.bias(ns.kps.pop.btn.errs()))

## Probabilities (numbers divided by total numbers of pairs)

# Within surveys
output$biasProbsKPsWtn = renderTable(find.bias(ns.kps.prb.wtn.errs()))

# Between surveys
output$biasProbsKPsBtn = renderTable(find.bias(ns.kps.prb.btn.errs()))

## Numbers among sampled animals

# Within surveys
output$biasNsKPsCapWtn = renderTable(find.bias(ns.kps.cap.wtn.errs()))

# Between surveys
output$biasNsKPsCapBtn = renderTable(find.bias(ns.kps.cap.btn.errs()))
