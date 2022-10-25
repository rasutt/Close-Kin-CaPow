# Outputs for kin-pairs sub-tabs in checks tab

## Output bias tables and box plots for values and prediction errors, can't be
## done with lists as want to keep reactives separate

VPE.srvr(
  "N", reactive(list(N.s.yrs())), reactive(list(est.ns.kps.pop()$wtn[, 1])), 
  N.errs, "In survey-years", "Population sizes"
)
# VPE.srvr(
#   "phi", 
#   reactive(list(matrix(
#     avg.phi.obs(), ncol = 1, dimnames = list(NULL, "All times")
#   ))), 
#   phi, phi.errs, "All", "Survival rates"
# )
VPE.srvr(
  "APs", ns.APs, ns.APs.preds, ns.APs.errs, wtn_btn_headings, "All pairs"
)
VPE.srvr(
  "SPs", ns.SPs, ns.SPs.preds, ns.SPs.errs, 
  c("All self-pairs", "Self-pairs with known parents"), "Self-pairs"
)
VPE.srvr(
  "POPs", ns.POPs, ns.POPs.preds, ns.POPs.errs, wtn_btn_headings, 
  "Parent-offspring pairs"
)
VPE.srvr(
  "SMPs", ns.SMPs.t[2:3], ns.SMPs.preds[2:3], ns.SMPs.errs[2:3], 
  wtn_btn_headings,
  # c("AgeKnwn", "WtnPop", "BtnPop", "BtnAgeKnwnPop"), 
  "Same-mother pairs"
)
VPE.srvr(
  "SFPs", ns.SFPs.t[2:3], ns.SFPs.preds[2:3], ns.SFPs.errs[2:3], 
  wtn_btn_headings,
  # c("AgeKnwn", "SameAge", "WtnPop", "BtnPop"), 
  "Same-father pairs"
)
VPE.srvr(
  "FSPs", ns.FSPs, ns.FSPs.preds, ns.FSPs.errs, wtn_btn_headings, 
  "Full-sibling pairs"
)
VPE.srvr(
  "HSPs", ns.HSPs, ns.HSPs.preds, ns.HSPs.errs, wtn_btn_headings, 
  "Half-sibling pairs"
)

