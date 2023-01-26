# Outputs for first study kin-pairs sub-tab of checks tab

# Get first study from list
fst.std = reactive(sim.lst()$hists.lst[[1]])

# Get attributes of first study
FS.atts = reactive(attributes(fst.std()))

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
  ns.kps.lst = FindNsKinPairs(k(), n.srvy.prs(), fst.std())
  
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
frst.KP.prbs.lst = reactive({
  print(str(tk.obj()))
  tk.obj()$report(tk.obj()$par)
})

# Kinpair probabilities for first study given true parameter values, from TMB
frmt.kp.prbs = function(kp.prbs) {
  mat = matrix(format(kp.prbs, digits = 3, scientific = T), nrow = k())
  df = data.frame(matrix(c(diag(mat), mat[upper.tri(mat)]), nrow = 1))
  colnames(df) = c(paste(srvy.yrs(), srvy.yrs(), sep = "-"), srvy.prs())
  df
}
output$firstKPPrbsTMB = renderTable({
  df = rbind(
    frmt.kp.prbs(frst.KP.prbs.lst()[["prbs_SPs"]]),
    frmt.kp.prbs(frst.KP.prbs.lst()[["prbs_POPs"]]),
    frmt.kp.prbs(frst.KP.prbs.lst()[["prbs_HSPs"]])
  )
  rownames(df) = knshp.chcs
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
frst.ns.kps.lst = reactive(FindNsKinPairs(k(), n.srvy.prs(), fst.std()))

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

