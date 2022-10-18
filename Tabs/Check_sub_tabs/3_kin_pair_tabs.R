# Outputs for kin-pairs sub-tabs in checks tab


# Function to display estimate errors as bias tables and box plots
show.errs = function(errs, kp.tp, kp.vrtns, kp.tp.nm) {
  tbls = paste0("bs", kp.tp, kp.vrtns)
  plts = paste0("errs", kp.tp, kp.vrtns)

  lapply(1:length(tbls), function(i) {
    output[[tbls[i]]] = renderTable(find.bias.srvy(errs()[[i]]))
    output[[plts[i]]] = renderPlot(plot.errs(errs()[[i]], kp.tp.nm))
  })
}

# Unknown-parents module-servers ----
unknPrntsServer("POPs.tab", pns.UPs)
unknPrntsServer("SMPs.tab", pns.UPs)
unknPrntsServer("SFPs.tab", pns.UPs)
unknPrntsServer("SibPs.tab", pns.UPs)

## Output bias tables and box plots for prediction errors ----

# Population sizes ----
show.errs(ns.wtn.errs, "Ns", "WtnPop", "Population sizes")

# Survival rates ----
show.errs(phi.errs, "Phi", "", "Survival rates")

# All-pairs ----
# output[["nsAPsWtnPop"]] = renderPlot(
#   plot.vals(ns.APs()[[1]], ns.APs.preds()[[1]], "All-pairs")
# )
# output[["nsAPsBtnPop"]] = renderPlot(
#   plot.vals(ns.APs()[[2]], ns.APs.preds()[[2]], "All-pairs")
# )
# show.errs(ns.APs.errs, "APs", c("WtnPop", "BtnPop"), "All-pairs")
VPE.srvr(
  "APs", ns.APs, ns.APs.preds, ns.APs.errs, c("WtnPop", "BtnPop"), "All-pairs"
)

# Self-pairs ----
show.errs(ns.SPs.errs, "SPs", "Pop", "Self-pairs")

# Parent-offspring pairs ----
show.errs(
  ns.POPs.errs, "POPs", c("WtnPop", "BtnPop"), "Parent-offspring pairs"
)

# Same-mother pairs ----
show.errs(
  ns.SMPs.errs, "SMPs", c("AgeKnwn", "WtnPop", "BtnPop", "BtnAgeKnwnPop"), 
  "Same-mother pairs"
)

# Same-father pairs ----
show.errs(
  ns.SFPs.errs, "SFPs", c("AgeKnwn", "SameAge", "WtnPop", "BtnPop"), 
  "Same-father pairs"
)

# Sibling-pairs ----
show.errs(ns.FSPs.errs, "FSPs", c("WtnPop", "BtnPop"), "Full-sibling pairs")
show.errs(ns.HSPs.errs, "HSPs", c("WtnPop", "BtnPop"), "Half-sibling pairs")

