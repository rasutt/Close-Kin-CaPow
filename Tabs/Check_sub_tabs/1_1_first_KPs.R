# Outputs for first study kin-pairs sub-tab of checks tab

# Get attributes of first study
FS.atts = reactive(attributes(frst.std()))

# Create general optimizer starting-values and bounds, NAs filled in below
ck.start = reactive({
  vec = c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
  names(vec) = c("rho", "phi", "N.final")
  vec
})

# True kinship likelihood TMB objective function - True kinships 
tk.obj = reactive({
  # Get numbers of animals captured in each survey
  ns.caps = FS.atts()$ns.caps
  
  # Find numbers of kin pairs
  ns.kps.lst = FindNsKinPairs(k(), n.srvy.prs(), frst.std())
  
  # Create TMB function
  MakeTMBObj(
    ck.start(), "true kinship",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), knshp_st_bool = all.knshps.bln,
    ns_SPs_btn = ns.kps.lst$btn[1, ], ns_POPs_wtn = ns.kps.lst$wtn[1, ], 
    ns_POPs_btn = ns.kps.lst$btn[2, ], ns_HSPs_wtn = ns.kps.lst$wtn[2, ],
    ns_HSPs_btn = ns.kps.lst$btn[3, ], ns_caps = ns.caps
  )
})

# Kinpair probabilities for true parameter values from TMB objective function
frst.KP.prbs.TVs.lst = reactive({
  tk.obj()$report(c(rho(), phi(), exp.N.fin()))
})

# Function to format kinpair probabilities for display
frmt.kp.prbs = function(kp.prbs.lst) {
  mat = t(
    sapply(
      kp.prbs.lst[c("prbs_SPs", "prbs_POPs", "prbs_HSPs")], 
      function(kp.prbs.mat) {
        # Pull out and combine values within and between surveys
        c(diag(kp.prbs.mat), kp.prbs.mat[upper.tri(kp.prbs.mat)])
      }
    )
  )
  mat = rbind(c(kp.prbs.lst[["exp_N_s_yrs"]], rep(NA, n.srvy.prs())), mat)
  df = data.frame(format(mat, digits = 3, scientific = T))
  colnames(df) = c(srvy.yrs(), srvy.prs())
  rownames(df) = c("Population size", knshp.chcs)
  df
}

# Kinpair probabilities for first study given true parameter values, from TMB
output$firstKPPrbsTMB = renderTable({
  frmt.kp.prbs(frst.KP.prbs.TVs.lst())
}, rownames = T)

# Kinpair probabilities for first study given true parameter values, from R
output$firstKPPrbsR = renderTable({
  wtn = pred.ns.kps.pop()$wtn
  btn = pred.ns.kps.pop()$btn

  prbs.wtn = cbind(wtn[, "N.s.yrs"], 0, wtn[, c("POPs", "HSPs")] / wtn[, "APs"])
  prbs.btn = cbind(NA, btn[, c("SPs", "POPs", "HSPs")] / btn[, "APs"])
  
  df = data.frame(format(t(rbind(prbs.wtn, prbs.btn)), dig = 3, sci = T))
  
  colnames(df) = c(srvy.yrs(), srvy.prs())
  rownames(df) = c("Population size", knshp.chcs)
  
  df
}, rownames = T)

# Function to format table of integers
frmt.tbl = function(data, rw.nms, cl.nms) {
  mode(data) = "integer"
  df = data.frame(matrix(data, ncol = length(cl.nms)), row.names = rw.nms)
  names(df) = cl.nms
  df
}

### Numbers of kin-pairs in whole population

## Within surveys

# Population sizes
FS.Nts.wtn = reactive(FS.atts()$N.t.vec[s.yr.inds()])

# Show table
output$firstNsKPsWtnPop = renderTable({
  SMPs = find.SMPs.wtn(FS.atts(), k())
  SFPs = find.SFPs.wtn(FS.atts(), k())
  FSPs = find.FSPs.wtn(FS.atts(), k())
  frmt.tbl(
    rbind(
      FS.Nts.wtn(), choose(FS.Nts.wtn(), 2), 
      find.POPs.wtn(FS.atts(), k()), 
      SMPs, SFPs, FSPs, SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps.pop.wtn, srvy.yrs()
  )
}, rownames = T)

## Between surveys

# Show table
output$firstNsKPsBtnPop = renderTable({
  SPs.prnts.kwn = find.SPs.prnts.kwn(FS.atts(), k())
  SMPs = find.SMSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  SFPs = find.SFSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  FSPs = find.FSSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  frmt.tbl(
    rbind(
      # Total numbers of pairs
      combn(FS.Nts.wtn(), 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2]), 
      
      # Self-pairs with parents unknown and known
      find.SPs(FS.atts(), k()), SPs.prnts.kwn, 
      
      # Other kin-pairs
      find.POPs.btn(FS.atts(), k()),
      SMPs, SFPs, FSPs, SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps.pop.btn, srvy.prs()
  )
}, rownames = T)

## Estimated numbers of kin-pairs in whole population

# Within surveys
output$firstEstNsKPsWtnPop = renderTable({
  frmt.tbl(t(pred.ns.kps.pop()$wtn), kp.tps.pop.wtn, srvy.yrs())
}, rownames = T)

# Between surveys - repeats predictions for self-pairs as haven't implemented
# for parents unknown
output$firstEstNsKPsBtnPop = renderTable({
  preds = t(pred.ns.kps.pop()$btn)
  frmt.tbl(preds[c(1:2, 2, 3:4, 6:8), ], kp.tps.pop.btn, srvy.prs())
}, rownames = T)

### Numbers of kin-pairs among sampled individuals

## Within surveys

# Find numbers of kin pairs
frst.ns.kps.lst = reactive(FindNsKinPairs(k(), n.srvy.prs(), frst.std()))

# Show table
output$firstNsKPsWtnSmp = renderTable({
  frmt.tbl(
    rbind(
      FS.atts()$ns.caps, 
      choose(FS.atts()$ns.caps, 2), frst.ns.kps.lst()$wtn
    ), 
    c(
      "Total number sampled", "All pairs", "Parent-offspring pairs",
      "Half-sibling pairs"
    ), 
    srvy.yrs()
  )
}, rownames = T)

## Between surveys

# Show table
output$firstNsKPsBtnSmp = renderTable({
  frmt.tbl(
    rbind(
      combn(
        FS.atts()$ns.caps, 2, function(ns.caps.pr) ns.caps.pr[1] * ns.caps.pr[2]
      ),
      frst.ns.kps.lst()$btn
    ), 
    c(
      "All pairs", "Self-pairs (all)", "Parent-offspring pairs",
      "Half-sibling pairs"
    ), 
    srvy.prs()
  )
}, rownames = T)

## Estimated numbers of kin-pairs among sampled individuals

# # Within surveys
# output$firstEstNsKPsWtnSmp = renderTable({
#   frmt.tbl(t(pred.ns.kps.pop()$wtn), kp.tps.pop.wtn, srvy.yrs())
# }, rownames = T)
# 
# # Between surveys - repeats predictions for self-pairs as haven't implemented
# # for parents unknown
# output$firstEstNsKPsBtnSmp = renderTable({
#   preds = t(pred.ns.kps.pop()$btn)
#   frmt.tbl(preds[c(1:2, 2, 3:4, 6:8), ], kp.tps.pop.btn, srvy.prs())
# }, rownames = T)

