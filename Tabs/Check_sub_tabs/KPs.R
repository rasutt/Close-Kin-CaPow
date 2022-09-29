# Outputs for kin-pairs sub-tabs in checks tab

## All-pairs

# Bias tables
output$biasAPsPopWtn = renderTable(find.bias.srvy(ns.APs.wtn.errs()))
output$biasAPsPopBtn = renderTable(find.bias.srvy(ns.APs.btn.errs()))

# Box plots
output$nsAPsWtnPop = renderPlot(nsKPsPlot(ns.APs.wtn.errs(), kp.tps.pop.wtn[2]))
output$nsAPsBtnPop = renderPlot(nsKPsPlot(ns.APs.btn.errs(), kp.tps.pop.btn[1]))

## Self-pairs

# Bias tables
output$biasSPsPop = renderTable(find.bias.srvy(ns.SPs.errs()))

# Box plots
output$nsSPsPop = renderPlot(nsKPsPlot(ns.SPs.errs(), kp.tps.pop.btn[2]))

## Parent-offspring pairs

## Same-mother pairs

# Bias tables
output$biasSMPsPopWtn = renderTable(find.bias.srvy(ns.SMPs.wtn.errs()))
output$biasSMPsPopBtn = renderTable(find.bias.srvy(ns.SMPs.btn.errs()))

# Box plots
output$nsSMPsWtnPop = renderPlot(nsKPsPlot(ns.SMPs.wtn.errs(), kp.tps.pop.wtn[3]))
output$nsSMPsBtnPop = renderPlot(nsKPsPlot(ns.SMPs.btn.errs(), kp.tps.pop.btn[3]))

## Same-father pairs

# Bias tables
output$biasSFPsPopWtn = renderTable(find.bias.srvy(ns.SFPs.wtn.errs()))
# output$biasSFPsPopBtn = renderTable(find.bias.srvy(ns.SFPs.btn.errs()))

# Box plots
output$nsSFPsWtnPop = renderPlot(nsKPsPlot(ns.SFPs.wtn.errs(), kp.tps.pop.wtn[4]))
# output$nsSFPsBtnPop = renderPlot(nsKPsPlot(ns.SFPs.btn.errs(), kp.tps.pop.btn[3]))
