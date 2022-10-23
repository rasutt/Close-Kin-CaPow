# Outputs for kin-pairs sub-tabs in checks tab

## Output bias tables and box plots for values and prediction errors, can't be
## done with lists as want to keep reactives separate

VPE.srvr(
  "N", reactive(list(N.s.yrs())), reactive(list(est.ns.kps.pop()$wtn[, 1])), 
  N.errs, "WtnPop", "Population sizes"
)
VPE.srvr(
  "phi", 
  reactive(list(matrix(
    avg.phi.obs(), ncol = 1, dimnames = list(NULL, "All times")
  ))), 
  phi, phi.errs, "All", "Survival rates"
)
VPE.srvr(
  "APs", ns.APs, ns.APs.preds, ns.APs.errs, c("WtnPop", "BtnPop"), "All pairs"
)
VPE.srvr(
  "SPs", ns.SPs, ns.SPs.preds, ns.SPs.errs, c("BtnPop", "PrntsKwn"), 
  "Self-pairs"
)
VPE.srvr(
  "POPs", ns.POPs, ns.POPs.preds, ns.POPs.errs, c("WtnPop", "BtnPop"), 
  "Parent-offspring pairs"
)
VPE.srvr(
  "SMPs", ns.SMPs.t, ns.SMPs.preds, ns.SMPs.errs, 
  c("AgeKnwn", "WtnPop", "BtnPop", "BtnAgeKnwnPop"), "Same-mother pairs"
)
VPE.srvr(
  "SFPs", ns.SFPs.t, ns.SFPs.preds, ns.SFPs.errs, 
  c("AgeKnwn", "SameAge", "WtnPop", "BtnPop"), "Same-father pairs"
)
VPE.srvr(
  "FSPs", ns.FSPs, ns.FSPs.preds, ns.FSPs.errs, 
  c("WtnPop", "BtnPop"), "Full-sibling pairs"
)
VPE.srvr(
  "HSPs", ns.HSPs, ns.HSPs.preds, ns.HSPs.errs, 
  c("WtnPop", "BtnPop"), "Half-sibling pairs"
)

