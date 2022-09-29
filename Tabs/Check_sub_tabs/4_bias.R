# Outputs for kin-pairs sub-tab of checks tab

### Kin-pair estimator biases (tables of average percentage differences)

## Numbers in whole population

# Temporal estimates
output$tempEstBiasNote = renderText(
  "These outputs seem to be higher variance, but 1000 studies seems to be enough."
)
output$biasNsKPsTemp = renderTable({
  find.est.bias(ns.kps.t.est.errs(), kp.tps.t)
})

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

# Show percentage of animals sampled for which the parents are unknown
output$percUnknPrnts = renderText({
  paste0(
    "Percentage of sampled animals with unknown parents: ", 
    round(mean(checks.lst()$prpn.prnts.unkn.vec) * 100, 1), "%"
  )
})

# Within surveys
output$biasNsKPsCapWtn = renderTable({
  find.est.bias(ns.kps.cap.wtn.est.errs(), kp.tps.cap.wtn)
})

# Between surveys
output$biasNsKPsCapBtn = renderTable({
  find.est.bias(ns.kps.cap.btn.est.errs(), kp.tps.cap.btn)
})
