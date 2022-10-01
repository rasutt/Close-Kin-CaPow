# Outputs for biases sub-tab of checks tab

# Show percentage of animals sampled for which the parents are unknown
output$percUnknPrnts = renderTable({
  df = data.frame(matrix(c(
    perc(mean(prpn.unkn.prnts()$prpn.unkn.prnts.wtn)),
    perc(mean(prpn.unkn.prnts()$prpn.unkn.prnts.btn))
  ), nrow = 1))
  names(df) = c("Surveys", "Survey-pairs")
  df
})

### Estimator biases (tables of average percentage differences)

## Numbers in whole population

# Temporal estimates
output$tempEstBiasNote = renderText(
  "These outputs seem to be higher variance, but 1000 studies seems to be enough."
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
