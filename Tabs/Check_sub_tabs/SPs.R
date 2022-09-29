# Outputs for self-pairs sub-tab in checks tab

# Bias tables
output$biasSPsPop = renderTable(find.bias.srvy(ns.SPs.errs()))

# Boxplots
output$nsSPsPop = renderPlot(nsKPsPlot(ns.SPs.errs(), kp.tps.pop.btn[2]))
