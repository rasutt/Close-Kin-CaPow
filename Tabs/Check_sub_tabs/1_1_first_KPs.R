# Outputs for first study kin-pairs sub-tab of checks tab

# Get attributes of first study
FS.atts = reactive(attributes(frst.std()))

# Function to format table of integers
FrmtTbl = function(data, rw.nms, cl.nms) {
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
  SPs.prnts.kwn = find.SPs.prnts.kwn(FS.atts(), k())
  SMPs = c(
    find.SMPs.wtn(FS.atts(), k()), 
    find.SMSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  )
  SFPs = c(
    find.SFPs.wtn(FS.atts(), k()),
    find.SFSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  )
  FSPs = c(
    find.FSPs.wtn(FS.atts(), k()),
    find.FSSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  )
  
  FrmtTbl(
    rbind(
      # Population sizes
      c(FS.Nts.wtn(), rep(NA, n.srvy.prs())),
      
      # Total numbers of pairs
      c(
        choose(FS.Nts.wtn(), 2),
        combn(FS.Nts.wtn(), 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2])
      ), 
      
      # Self-pairs with parents unknown and known
      c(rep(NA, k()), find.SPs(FS.atts(), k())), 
      c(rep(NA, k()), SPs.prnts.kwn), 
      
      # Other kin-pairs
      c(find.POPs.wtn(FS.atts(), k()), find.POPs.btn(FS.atts(), k())),
      SMPs, SFPs, FSPs, SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps, c(srvy.yrs(), srvy.prs())
  )
}, rownames = T)

## Between surveys

# Show table
output$firstNsKPsBtnPop = renderTable({
  SPs.prnts.kwn = find.SPs.prnts.kwn(FS.atts(), k())
  SMPs = find.SMSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  SFPs = find.SFSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  FSPs = find.FSSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  FrmtTbl(
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
  FrmtTbl(t(pred.ns.kps.pop()$wtn), kp.tps.pop.wtn, srvy.yrs())
}, rownames = T)

# Between surveys - repeats predictions for self-pairs as haven't implemented
# for parents unknown
output$firstEstNsKPsBtnPop = renderTable({
  preds = t(pred.ns.kps.pop()$btn)
  FrmtTbl(preds[c(1:2, 2, 3:4, 6:8), ], kp.tps.pop.btn, srvy.prs())
}, rownames = T)

### Numbers of kin-pairs among sampled individuals

## Within surveys

# Find numbers of kin pairs
frst.ns.kps.lst = reactive(FindNsKinPairs(k(), n.srvy.prs(), frst.std()))

# Show table
output$firstNsKPsWtnSmp = renderTable({
  FrmtTbl(
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
  FrmtTbl(
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

