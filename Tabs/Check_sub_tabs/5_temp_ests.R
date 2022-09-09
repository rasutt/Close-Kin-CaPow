# Outputs for temporal estimates sub-tab of checks tab

## Temporal estimates vs observed averages (mainly for debugging)

output$tempEstNote = renderText(
  "These outputs are for checking intermediate results that are used in the derivation of estimators of numbers of close-kin pairs, and which depend on time in some way."
)

# Numbers of same-mother/father pairs in the population including animals born
# in each year in the population history
output$nsKPsTemp = renderTable({
  # Find average values simulated for each year and kin-pair type, and combine
  # with years and estimates
  df = data.frame(
    yrs.chk.t(), apply(checks.lst()$kps.t.arr, 2:3, mean), est.ns.kps.t()
  )
  colnames(df) = 
    c("Year", paste0(rep(c("Avg", "Exp"), each = n.kp.tps.t), rep(kp.tps.t, 2)))

  # Reorder for output with averages next to estimates
  df[n.yrs.chk.t():1, 
     c(1, rep(seq(1:n.kp.tps.t), each = 2) + c(0, n.kp.tps.t) + 1)]
}, digits = 1)

