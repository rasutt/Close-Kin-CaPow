# Outputs for kin-pairs sub-tabs in checks tab

# Function to find estimate errors as proportions of estimates
find.errs = function(vals, ests, samp = F) {
  # For sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) ests = rep(ests, each = n.sims())
  vals / ests - 1
}

# Find errors for estimates in sets
ns.kps.pop.wtn.errs = reactive({
  find.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn)
})
ns.kps.pop.btn.errs = reactive({
  find.errs(checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn)
})
ns.kps.t.errs = reactive({
  find.errs(checks.lst()$kps.t.arr, est.ns.kps.t())
})
ns.kps.prb.wtn.errs = reactive({
  # Remove population sizes and total numbers of pairs then divide by the latter
  find.errs(
    checks.lst()$ns.kps.pop.wtn.arr[, , -1:-2] / 
      array(
        rep(checks.lst()$ns.kps.pop.wtn.arr[, , 2], n.kp.tps.prb.wtn), 
        c(n.sims(), k(), n.kp.tps.prb.wtn)
      ), 
    est.ns.kps.pop.lst()$wtn[, -1:-2] / est.ns.kps.pop.lst()$wtn[, 2]
  )
})
ns.kps.prb.btn.errs = reactive({
  # Remove population sizes and total numbers of pairs then divide by the latter
  # (drop = F retains array dimensions when only one survey pair)
  find.errs(
    checks.lst()$ns.kps.pop.btn.arr[, , -1, drop = F] / 
      array(
        rep(checks.lst()$ns.kps.pop.btn.arr[, , 1], n.kp.tps.prb.btn), 
        c(n.sims(), n.srvy.prs(), n.kp.tps.prb.btn)
      ), 
    est.ns.kps.pop.lst()$btn[, -1] / est.ns.kps.pop.lst()$btn[, 1]
  )
})
# ns.kps.cap.wtn.errs = reactive({
#   find.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn, T)
# })
# ns.kps.cap.btn.errs = reactive({
#   find.errs(checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn, T)
# })

# Function to find biases over all surveys from array of proportional errors for
# multiple estimators
find.bias = function(errs) {
  df = data.frame(matrix(perc(colMeans(errs, dims = 2)), nrow = 1))
  names(df) = dimnames(errs)[["kp.type"]]
  df
}

# Function to find biases in each survey from matrix of proportional errors for
# a single estimator
find.bias.srvy = function(errs) {
  df = data.frame(matrix(perc(colMeans(errs)), nrow = 1))
  names(df) = colnames(errs)
  df
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

# Find errors for population sizes
ns.wtn.errs = reactive(find.errs(N.s.yrs(), est.ns.kps.pop.lst()$wtn[, 1]))

## Find errors and output bias tables and box plots for numbers of kin-pairs separately

# All-pairs
ns.APs.wtn.errs = reactive({
  find.errs(ns.APs.wtn.pop(), est.ns.kps.pop.lst()$wtn[, 2])
})
ns.APs.btn.errs = reactive({
  find.errs(ns.APs.btn.pop(), est.ns.kps.pop.lst()$btn[, 1])
})
output$biasAPsPopWtn = renderTable(find.bias.srvy(ns.APs.wtn.errs()))
output$biasAPsPopBtn = renderTable(find.bias.srvy(ns.APs.btn.errs()))
output$nsAPsWtnPop = renderPlot(nsKPsPlot(ns.APs.wtn.errs(), kp.tps.pop.wtn[2]))
output$nsAPsBtnPop = renderPlot(nsKPsPlot(ns.APs.btn.errs(), kp.tps.pop.btn[1]))

# Self-pairs
ns.SPs.errs = reactive(find.errs(ns.SPs(), est.ns.kps.pop.lst()$btn[, 2]))
output$biasSPsPop = renderTable(find.bias.srvy(ns.SPs.errs()))
output$nsSPsPop = renderPlot(nsKPsPlot(ns.SPs.errs(), kp.tps.pop.btn[2]))

# Parent-offspring pairs
ns.POPs.wtn.errs = reactive({
  find.errs(ns.POPs()[["ns.POPs.wtn"]], est.ns.kps.pop.lst()$wtn[, 3])
})
ns.POPs.btn.errs = reactive({
  find.errs(ns.POPs()[["ns.POPs.btn"]], est.ns.kps.pop.lst()$btn[, 3])
})
output$biasPOPsPopWtn = renderTable(find.bias.srvy(ns.POPs.wtn.errs()))
output$biasPOPsPopBtn = renderTable(find.bias.srvy(ns.POPs.btn.errs()))
output$nsPOPsWtnPop =
  renderPlot(nsKPsPlot(ns.POPs.wtn.errs(), kp.tps.pop.wtn[3]))
output$nsPOPsBtnPop =
  renderPlot(nsKPsPlot(ns.POPs.btn.errs(), kp.tps.pop.btn[3]))

# Same-mother pairs
ns.SMPs.age.knwn.errs = reactive({
  find.errs(ns.SMPs()[["ns.SMPs.age.knwn"]], est.ns.kps.t()[, 1])
})
ns.SMPs.wtn.errs = reactive({
  find.errs(ns.SMPs()[["ns.SMPs.wtn"]], est.ns.kps.pop.lst()$wtn[, 4])
})
ns.SMPs.btn.errs = reactive({
  find.errs(ns.SMPs()[["ns.SMPs.btn"]], est.ns.kps.pop.lst()$btn[, 4])
})
output$biasSMPsAgeKnwn = renderTable(find.bias.srvy(ns.SMPs.age.knwn.errs()))
output$biasSMPsPopWtn = renderTable(find.bias.srvy(ns.SMPs.wtn.errs()))
output$biasSMPsPopBtn = renderTable(find.bias.srvy(ns.SMPs.btn.errs()))
output$nsSMPsAgeKnwn = 
  renderPlot(nsKPsPlot(ns.SMPs.age.knwn.errs(), kp.tps.t[1]))
output$nsSMPsWtnPop = 
  renderPlot(nsKPsPlot(ns.SMPs.wtn.errs(), kp.tps.pop.wtn[4]))
output$nsSMPsBtnPop = 
  renderPlot(nsKPsPlot(ns.SMPs.btn.errs(), kp.tps.pop.btn[4]))

# Same-father pairs
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
output$biasFSPsPopWtn = renderTable(find.bias.srvy(ns.FSPs.wtn.errs()))
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
