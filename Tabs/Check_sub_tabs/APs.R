# Outputs for all-pairs sub-tab in checks tab

# Bias tables
output$biasAPsPopWtn = renderTable({
  find.est.bias.srvy(ns.kps.pop.wtn.est.errs()[, , 2])
})
output$biasAPsPopBtn = renderTable({
  find.est.bias.srvy(ns.kps.pop.btn.est.errs()[, , 1])
})

# Boxplots
output$nsAPsWtnPop = renderPlot(
  nsKPsPlot(ns.kps.pop.wtn.est.errs()[, , 2], kp.tps.pop.wtn[2])
)
output$nsAPsBtnPop = renderPlot(
  nsKPsPlot(ns.kps.pop.btn.est.errs()[, , 1], kp.tps.pop.btn[1])
)
