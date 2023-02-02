# Outputs for kin-pairs sub-tabs in checks tab

# Bias tables and box plots for values and proportional errors, can't be done
# with lists as want to keep reactives separate
VPE.srvr(
  "N", reactive(list(N.s.yrs())), 
  reactive(list(pred.ns.kps.pop()$wtn[, 1])), 
  N.errs, "In survey-years"
)
VPE.srvr.rglr("APs", ns.APs, preds.lst, ns.APs.errs)
VPE.srvr(
  "SPs", ns.SPs, ns.SPs.preds, ns.SPs.errs, 
  c("All self-pairs", "Self-pairs with known parents")
)
VPE.srvr.rglr("POPs", ns.POPs, preds.lst, ns.POPs.errs)
VPE.srvr.rglr("SMPs", ns.SMPs, preds.lst, ns.SMPs.errs)
VPE.srvr.rglr("SFPs", ns.SFPs, preds.lst, ns.SFPs.errs)
VPE.srvr.rglr("FSPs", ns.FSPs, preds.lst, ns.FSPs.errs)
VPE.srvr.rglr("HSPs", ns.HSPs, preds.lst, ns.HSPs.errs)

