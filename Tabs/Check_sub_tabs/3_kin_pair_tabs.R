# Outputs for kin-pairs sub-tabs in checks tab

# Function to find estimate errors as proportions of estimates
find.errs = function(vals, ests, samp = F) {
  # For sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) ests = rep(ests, each = n.sims())
  vals / ests - 1
}

# Function to find list of estimate errors
l.fnd.errs = function(ns, preds) {
  lapply(1:length(ns), function(i) find.errs(ns[[i]], preds[[i]]))
}

# Functions to display estimate errors as bias tables and box plots
show.errs = function(errs, kp.tp, kp.vrtns, kp.tp.nm) {
  tbls = paste0("bs", kp.tp, kp.vrtns)
  plts = paste0("errs", kp.tp, kp.vrtns)

  lapply(1:length(tbls), function(i) {
    output[[tbls[i]]] = renderTable(find.bias.srvy(errs()[[i]]))
    output[[plts[i]]] = renderPlot(nsKPsPlot(errs()[[i]], kp.tp.nm))
  })
}

# Function to plot simulated versus expected numbers of kin-pairs for one type
# of kin-pair ----
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

# Unknown-parents module-servers ----
unknPrntsServer("POPs.tab", pns.UPs)
unknPrntsServer("SMPs.tab", pns.UPs)
unknPrntsServer("SFPs.tab", pns.UPs)
unknPrntsServer("SibPs.tab", pns.UPs)

## ----
## Find errors and output bias tables and box plots for estimators separately

# Population sizes ----
ns.wtn.errs = reactive(list(find.errs(N.s.yrs(), est.ns.kps.pop()$wtn[, 1])))
show.errs(ns.wtn.errs, "Ns", "WtnPop", "Population sizes")

# Survival rates ----
phi.errs = reactive({
  list(matrix(
    find.errs(avg.phi.obs(), phi()), n.sims(), 
    dimnames = list(NULL, "All times")
  ))
})
show.errs(phi.errs, "Phi", "", "Survival rates")

# All-pairs ----
ns.APs.preds = reactive({
  list(est.ns.kps.pop()$wtn[, "APs"], est.ns.kps.pop()$btn[, "APs"])
})
ns.APs.errs = reactive(l.fnd.errs(ns.APs(), ns.APs.preds()))
show.errs(ns.APs.errs, "APs", c("WtnPop", "BtnPop"), "All-pairs")

# Self-pairs ----
ns.SPs.errs = reactive({
  list(find.errs(ns.SPs(), est.ns.kps.pop()$btn[, "SPs"]))
})
show.errs(ns.SPs.errs, "SPs", "Pop", "Self-pairs")

# Parent-offspring pairs ----
ns.POPs.preds = reactive(list(
  est.ns.kps.pop()$wtn[, "POPs"], est.ns.kps.pop()$btn[, "POPs"]
))
ns.POPs.errs = reactive(l.fnd.errs(ns.POPs(), ns.POPs.preds()))
show.errs(
  ns.POPs.errs, "POPs", c("WtnPop", "BtnPop"), "Parent-offspring pairs"
)

# Same-mother pairs ----
ns.SMPs.preds = reactive(list(
  est.ns.kps.t()[, "SMPs.kwn.age"],  est.ns.kps.pop()$wtn[, "SMPs"], 
  est.ns.kps.pop()$btn[, "SMPs"], est.ns.kps.pop()$btn[, "SMPs.kwn.age"]
))
ns.SMPs.errs = reactive(l.fnd.errs(ns.SMPs.t(), ns.SMPs.preds()))
show.errs(
  ns.SMPs.errs, "SMPs", c("AgeKnwn", "WtnPop", "BtnPop", "BtnAgeKnwnPop"), 
  "Same-mother pairs"
)

# Same-father pairs
ns.SFPs.preds = reactive(list(
  est.ns.kps.t()[, "SFPs.kwn.age"],  est.ns.kps.t()[, "SFPs.sm.age"], 
  est.ns.kps.pop()$wtn[, "SFPs"], est.ns.kps.pop()$btn[, "SFPs"]
))
ns.SFPs.errs = reactive(l.fnd.errs(ns.SFPs.t(), ns.SFPs.preds()))
show.errs(
  ns.SFPs.errs, "SFPs", c("AgeKnwn", "SameAge", "WtnPop", "BtnPop"), 
  "Same-father pairs"
)

# Sibling-pairs
ns.SibPs.preds = reactive(list(
  est.ns.kps.pop()$wtn[, "FSPs"], est.ns.kps.pop()$btn[, "FSPs"], 
  est.ns.kps.pop()$wtn[, "HSPs"], est.ns.kps.pop()$btn[, "HSPs"]
))
ns.FSPs.errs = reactive(l.fnd.errs(ns.SibPs()[1:2], ns.SibPs.preds()[1:2]))
ns.HSPs.errs = reactive(l.fnd.errs(ns.SibPs()[3:4], ns.SibPs.preds()[3:4]))
show.errs(ns.FSPs.errs, "FSPs", c("WtnPop", "BtnPop"), "Full-sibling pairs")
show.errs(ns.HSPs.errs, "HSPs", c("WtnPop", "BtnPop"), "Half-sibling pairs")

