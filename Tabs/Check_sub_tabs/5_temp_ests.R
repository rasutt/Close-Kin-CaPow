# Outputs for temporal estimates sub-tab of checks tab

## Temporal estimates vs observed averages (mainly for debugging)

output$tempEstNote = renderText(
  "These outputs are for checking intermediate results that are used in the derivation of estimators of numbers of close-kin pairs, and which depend on time in some way."
)

# Numbers of same-mother/father pairs in the population including animals born
# in each year in the population history
output$nsKPsTemp = renderTable({
  # Find average values simulated for each each and each kin-pair type
  mean.kps.t = apply(checks.lst()$kps.t.arr, 2:3, mean)
  
  # Combine and reorder for output with averages next to estimates
  # df = data.frame(sim.yrs()[-c(1, hist.len())], mean.kps.t, exp.kps.t)
  df = data.frame(
    sim.yrs()[(hist.len() - n.yrs.chk.t):(hist.len() - 1)], mean.kps.t, 
    est.ns.kps.t()
  )
  colnames(df) = 
    c("Year", paste0(rep(c("Avg", "Exp"), each = n.kp.tps.t), rep(kp.tps.t, 2)))
  df[n.yrs.chk.t:1,
  # df[(hist.len() - 2):2, 
     c(1, rep(seq(1:n.kp.tps.t), each = 2) + c(0, n.kp.tps.t) + 1)]
}, digits = 1)

