# Outputs for all-pairs sub-tab in checks tab

# Bias tables
output$biasAPsPopWtn = renderTable(find.bias.srvy(ns.APs.wtn.errs()))
output$biasAPsPopBtn = renderTable(find.bias.srvy(ns.APs.btn.errs()))

# Boxplots
output$nsAPsWtnPop = renderPlot(nsKPsPlot(ns.APs.wtn.errs(), kp.tps.pop.wtn[2]))
output$nsAPsBtnPop = renderPlot(nsKPsPlot(ns.APs.btn.errs(), kp.tps.pop.btn[1]))
