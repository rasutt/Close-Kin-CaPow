# Outputs for temporal estimates sub-tab of checks tab

# Table of temporal estimates vs observed averages (mainly for debugging)
output$tempEstNote = renderText(
  "These outputs are for checking intermediate results that are used in the derivation of estimators of numbers of close-kin pairs, and which depend on time in some way.  SMP{t=fnl,b1,b2=fnl} represents the number of same-mother pairs in the final year, with one born in the year indicated, and one born in the final year."
)

# Numbers of same-mother/father pairs in the population including animals born
# in each year in the population history
output$nsKPsTemp = renderTable({
  # Find average values simulated for each year and kin-pair type, and combine
  # with estimates
  df = data.frame(apply(checks.lst()$kps.t.arr, 2:3, mean), est.ns.kps.t())
  colnames(df) = paste0(rep(c("Avg", "Est"), each = n.kp.tps.t), kp.tps.t)

  # Reorder for output with last year at top and averages next to estimates
  df[n.yrs.chk.t():1, rep(seq(1:n.kp.tps.t), each = 2) + c(0, n.kp.tps.t)]
}, digits = 1, rownames = T)
