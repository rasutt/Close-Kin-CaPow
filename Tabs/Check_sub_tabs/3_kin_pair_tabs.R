# Outputs for kin-pairs sub-tabs in checks tab

# Function to find estimate errors as proportions of estimates
find.errs = function(vals, ests, samp = F) {
  # For sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) ests = rep(ests, each = n.sims())
  vals / ests - 1
}

# Function to plot simulated versus expected numbers of kin-pairs for one type
# of kin-pair
nsKPsPlot = function(errs, kp.type) {
  boxplot(
    errs, main = kp.type, xlab = names(dimnames(errs))[2],
    ylab = "Proportional errors", show.names = T
  )
  abline(h = 0, col = 'red')
  abline(h = mean(errs), col = 'blue')
  legend(
    "topleft", col = c(2, 4), lty = 1,
    legend = c("Estimated error (zero)", "Average error"),
  )
}

# Find errors for population sizes and survival rates
ns.wtn.errs = reactive({
  find.errs(N.s.yrs(), est.ns.kps.pop.lst()$wtn[, 1])
})
phi.errs = reactive({
  matrix(
    find.errs(avg.phi.obs(), phi()), n.sims(), 
    dimnames = list(NULL, "All times")
  )
})
output$biasPhi = renderTable(find.bias.srvy(phi.errs()))
output$errsPhi = renderPlot(nsKPsPlot(phi.errs(), "Phi"))

## Find errors and output bias tables and box plots for numbers of kin-pairs
## separately

# Changing to list approach
# All-pairs
ns.APs = reactive(list(ns.APs.wtn.pop(), ns.APs.btn.pop()))
ns.APs.preds = reactive(list(
  est.ns.kps.pop.lst()$wtn[, 2], est.ns.kps.pop.lst()$btn[, 1]
))
ns.APs.errs = reactive({
  lapply(1:length(ns.APs()), function(i) {
    find.errs(ns.APs()[[i]], ns.APs.preds()[[i]])
  })
})
ns.APs.bs.tbls = c("biasAPsPopWtn", "biasAPsPopBtn")
ns.APs.err.plts = c("nsAPsWtnPop", "nsAPsBtnPop")
lapply(1:length(ns.APs.bs.tbls), function(i) {
  output[[ns.APs.bs.tbls[i]]] = renderTable(find.bias.srvy(ns.APs.errs()[[i]]))
  output[[ns.APs.err.plts[i]]] = 
    renderPlot(nsKPsPlot(ns.APs.errs()[[i]], "All-pairs"))
})

# Parent-offspring pairs
ns.POPs.preds = reactive(list(
  est.ns.kps.pop.lst()$wtn[, 3], est.ns.kps.pop.lst()$btn[, 3]
))
ns.POPs.errs = reactive({
  lapply(1:length(ns.POPs()), function(i) {
    find.errs(ns.POPs()[[i]], ns.POPs.preds()[[i]])
  })
})
ns.POPs.bs.tbls = c("biasPOPsPopWtn", "biasPOPsPopBtn")
ns.POPs.err.plts = c("nsPOPsWtnPop", "nsPOPsBtnPop")
lapply(1:length(ns.POPs.bs.tbls), function(i) {
  output[[ns.POPs.bs.tbls[i]]] = renderTable(find.bias.srvy(ns.POPs.errs()[[i]]))
  output[[ns.POPs.err.plts[i]]] = 
    renderPlot(nsKPsPlot(ns.POPs.errs()[[i]], "Parent-offspring pairs"))
})

# Self-pairs
ns.SPs.errs = reactive(find.errs(ns.SPs(), est.ns.kps.pop.lst()$btn[, 2]))
output[["biasSPsPop"]] = renderTable(find.bias.srvy(ns.SPs.errs()))
output[["nsSPsPop"]] = renderPlot(nsKPsPlot(ns.SPs.errs(), "Self-pairs"))

# Same-mother pairs
ns.SMPs.preds = reactive(list(
  est.ns.kps.t()[, 1],  est.ns.kps.pop.lst()$wtn[, 4], 
  est.ns.kps.pop.lst()$btn[, 4]
))
ns.SMPs.errs = reactive({
  lapply(1:length(ns.SMPs()), function(i) {
    find.errs(ns.SMPs()[[i]], ns.SMPs.preds()[[i]])
  })
})
ns.SMPs.bs.tbls = c("biasSMPsAgeKnwn", "biasSMPsPopWtn", "biasSMPsPopBtn")
ns.SMPs.err.plts = c("nsSMPsAgeKnwn", "nsSMPsWtnPop", "nsSMPsBtnPop")
lapply(1:length(ns.SMPs.bs.tbls), function(i) {
  output[[ns.SMPs.bs.tbls[i]]] = renderTable(find.bias.srvy(ns.SMPs.errs()[[i]]))
  output[[ns.SMPs.err.plts[i]]] = 
    renderPlot(nsKPsPlot(ns.SMPs.errs()[[i]], "Same-mother pairs"))
})

# Parent-offspring pairs
unknPrntsServer("POPs.tab", pns.UPs)

# Same-mother pairs
unknPrntsServer("SMPs.tab", pns.UPs)
# ns.SMPs.age.knwn.errs = reactive({
#   find.errs(ns.SMPs()[["ns.SMPs.age.knwn"]], est.ns.kps.t()[, 1])
# })
# ns.SMPs.wtn.errs = reactive({
#   find.errs(ns.SMPs()[["ns.SMPs.wtn"]], est.ns.kps.pop.lst()$wtn[, 4])
# })
# ns.SMPs.btn.errs = reactive({
#   find.errs(ns.SMPs()[["ns.SMPs.btn"]], est.ns.kps.pop.lst()$btn[, 4])
# })
# output$biasSMPsAgeKnwn = renderTable(find.bias.srvy(ns.SMPs.age.knwn.errs()))
# output$biasSMPsPopWtn = renderTable(find.bias.srvy(ns.SMPs.wtn.errs()))
# output$biasSMPsPopBtn = renderTable(find.bias.srvy(ns.SMPs.btn.errs()))
# output$nsSMPsAgeKnwn = 
#   renderPlot(nsKPsPlot(ns.SMPs.age.knwn.errs(), kp.tps.t[1]))
# output$nsSMPsWtnPop = 
#   renderPlot(nsKPsPlot(ns.SMPs.wtn.errs(), kp.tps.pop.wtn[4]))
# output$nsSMPsBtnPop = 
#   renderPlot(nsKPsPlot(ns.SMPs.btn.errs(), kp.tps.pop.btn[4]))

# Same-father pairs
unknPrntsServer("SFPs.tab", pns.UPs)
ns.SFPs.age.knwn.errs = reactive({
  find.errs(ns.SFPs()[["ns.SFPs.age.knwn"]], est.ns.kps.t()[, 2])
})
ns.SFPs.same.age.errs = reactive({
  find.errs(ns.SFPs()[["ns.SFPs.same.age"]], est.ns.kps.t()[, 3])
})
ns.SFPs.wtn.errs = reactive({
  find.errs(ns.SFPs()[["ns.SFPs.wtn"]], est.ns.kps.pop.lst()$wtn[, 5])
})
# ns.SFPs.btn.errs = reactive({
#   find.errs(ns.SFPs()[["ns.SFPs.btn"]], est.ns.kps.pop.lst()$btn[, 6])
# })
output$biasSFPsAgeKnwn = renderTable(find.bias.srvy(ns.SFPs.age.knwn.errs()))
output$biasSFPsSameAge = renderTable(find.bias.srvy(ns.SFPs.same.age.errs()))
output$biasSFPsPopWtn = renderTable(find.bias.srvy(ns.SFPs.wtn.errs()))
# output$biasSFPsPopBtn = renderTable(find.bias.srvy(ns.SFPs.btn.errs()))
output$nsSFPsAgeKnwn = 
  renderPlot(nsKPsPlot(ns.SFPs.age.knwn.errs(), kp.tps.t[2]))
output$nsSFPsSameAge = 
  renderPlot(nsKPsPlot(ns.SFPs.same.age.errs(), kp.tps.t[3]))
output$nsSFPsWtnPop = 
  renderPlot(nsKPsPlot(ns.SFPs.wtn.errs(), kp.tps.pop.wtn[5]))
# output$nsSFPsBtnPop =
#   renderPlot(nsKPsPlot(ns.SFPs.btn.errs(), kp.tps.pop.btn[6]))

# Sibling-pairs
unknPrntsServer("SibPs.tab", pns.UPs)
ns.FSPs.wtn.errs = reactive({
  find.errs(ns.SibPs()[["ns.FSPs.wtn"]], est.ns.kps.pop.lst()$wtn[, 6])
})
ns.HSPs.wtn.errs = reactive({
  find.errs(ns.SibPs()[["ns.HSPs.wtn"]], est.ns.kps.pop.lst()$wtn[, 7])
})
# ns.HSPs.wtn.errs = reactive({
#   find.errs(ns.HSPs()[["ns.HSPs.wtn"]], est.ns.kps.pop.lst()$wtn[, 5])
# })
# ns.FSPs.btn.errs = reactive({
#   find.errs(ns.FSPs()[["ns.FSPs.btn"]], est.ns.kps.pop.lst()$btn[, 6])
# })
output$biasFSPsPopWtn = renderTable({
  find.bias.srvy(ns.FSPs.wtn.errs())
})
output$biasHSPsPopWtn = renderTable(find.bias.srvy(ns.HSPs.wtn.errs()))
# output$biasHSPsPopBtn = renderTable(find.bias.srvy(ns.HSPs.btn.errs()))
# output$biasFSPsPopBtn = renderTable(find.bias.srvy(ns.FSPs.btn.errs()))
output$nsFSPsWtnPop = 
  renderPlot(nsKPsPlot(ns.FSPs.wtn.errs(), kp.tps.pop.wtn[6]))
output$nsHSPsWtnPop = 
  renderPlot(nsKPsPlot(ns.HSPs.wtn.errs(), kp.tps.pop.wtn[7]))
# output$nsHSPsBtnPop =
#   renderPlot(nsKPsPlot(ns.HSPs.btn.errs(), kp.tps.pop.btn[6]))
# output$nsFSPsBtnPop =
#   renderPlot(nsKPsPlot(ns.FSPs.btn.errs(), kp.tps.pop.btn[6]))
