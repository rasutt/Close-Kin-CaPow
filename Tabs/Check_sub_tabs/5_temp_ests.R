# Outputs for temporal estimates sub-tab of checks tab

## Temporal estimates vs observed averages (mainly for debugging)

# Numbers of same-mother/father pairs in the population including animals born
# in each year in the population history
output$nsKPsTemp = renderTable({
  # Find average values simulated
  mean.kps.t = t(apply(checks.lst()$kps.t.arr, 2:3, mean))
  
  # Find estimated values
  exp.kps.t = FindExpNsKPsT(
    exp.N.fin(), phi(), lambda(), alpha(), hist.len(), exp.N.t()
  )
  
  # Combine for output
  df = rbind(mean.kps.t, exp.kps.t)
  rownames(df) = paste0(
    rep(c("Avg", "Exp"), each = n.kp.t.tps),
    rep(kp.t.tps, 2)
  )
  colnames(df) = sim.yrs()[-c(1, hist.len())]
  df[rep(seq(1:n.kp.t.tps), each = 2) + c(0, n.kp.t.tps), ]
}, rownames = T, digits = 1)

