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
VPE.srvr(
  "APs", ns.APs, ns.APs.preds, ns.APs.errs, c("WtnPop", "BtnPop"), "All-pairs"
)
# Self-pairs ----
VPE.srvr(
  "SPs", ns.SPs, ns.SPs.preds, ns.SPs.errs, c("BtnPop", "PrntsKwn"), "Self-pairs"
)
# Parent-offspring pairs ----
VPE.srvr(
  "POPs", ns.POPs, ns.POPs.preds, ns.POPs.errs, c("WtnPop", "BtnPop"), 
  "Parent-offspring pairs"
)
# Same-mother pairs ----
VPE.srvr(
  "SMPs", ns.SMPs.t, ns.SMPs.preds, ns.SMPs.errs, 
  c("AgeKnwn", "WtnPop", "BtnPop", "BtnAgeKnwnPop"), "Same-mother pairs"
)
# Same-father pairs ----
VPE.srvr(
  "SFPs", ns.SFPs.t, ns.SFPs.preds, ns.SFPs.errs, 
  c("AgeKnwn", "SameAge", "WtnPop", "BtnPop"), "Same-father pairs"
)
# Sibling-pairs ----
VPE.srvr(
  "FSPs", ns.FSPs, ns.FSPs.preds, ns.FSPs.errs, 
  c("WtnPop", "BtnPop"), "Full-sibling pairs"
)
VPE.srvr(
  "HSPs", ns.HSPs, ns.HSPs.preds, ns.HSPs.errs, 
  c("WtnPop", "BtnPop"), "Half-sibling pairs"
)

