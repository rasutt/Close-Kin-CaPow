# Outputs for kin-pairs sub-tabs in checks tab

# Bias tables and box plots for values and proportional errors, can't be done
# with lists as want to keep reactives separate
VPE.srvr(
  "N", reactive(list(N.s.yrs())), 
  reactive(list(est.ns.kps.pop()$wtn[, 1])), 
  N.errs, "In survey-years"
)
# VPE.srvr(
#   "phi", 
#   reactive(list(matrix(
#     avg.phi.obs(), ncol = 1, dimnames = list(NULL, "All times")
#   ))), 
#   phi, phi.errs, "All", "Survival rates"
# )
VPE.srvr("APs", ns.APs, ns.APs.preds, ns.APs.errs)
# VPE.srvr("APs", ns.APs, preds, ns.APs.errs)
VPE.srvr(
  "SPs", ns.SPs, ns.SPs.preds, ns.SPs.errs, 
  c("All self-pairs", "Self-pairs with known parents")
)
VPE.srvr("POPs", ns.POPs, ns.POPs.preds, ns.POPs.errs)
VPE.srvr("SMPs", ns.SMPs, ns.SMPs.preds, ns.SMPs.errs)
VPE.srvr("SFPs", ns.SFPs, ns.SFPs.preds, ns.SFPs.errs)
VPE.srvr("FSPs", ns.FSPs, ns.FSPs.preds, ns.FSPs.errs)
VPE.srvr("HSPs", ns.HSPs, ns.HSPs.preds, ns.HSPs.errs)

