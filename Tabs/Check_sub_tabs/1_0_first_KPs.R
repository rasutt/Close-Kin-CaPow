# Outputs for first study kin-pairs sub-tab of checks tab

# Get first study from list
fst.std = reactive(sim.lst()$hists.lst[[1]])

# Get attributes of first study
FS.atts = reactive(attributes(fst.std()))

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

# Total numbers of pairs
FS.APs.wtn = reactive(choose(FS.atts()$N.t.vec[s.yr.inds()], 2))

# Numbers of parent-offspring pairs
FS.POPs.wtn = reactive(find.POPs.wtn(FS.atts(), k()))

# Show table
output$firstNsKPsWtnPop = renderTable({
  SMPs = find.SMPs.wtn(FS.atts(), k())
  SFPs = find.SFPs.wtn(FS.atts(), k())
  FSPs = find.FSPs.wtn(FS.atts(), k())
  frmt.tbl(
    rbind(
      FS.Nts.wtn(), FS.APs.wtn(), FS.POPs.wtn(), SMPs, SFPs, FSPs, 
      SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps.pop.wtn, srvy.yrs()
  )
}, rownames = T)

## Between surveys

# Total numbers of pairs
FS.APs.btn = reactive({
  combn(FS.Nts.wtn(), 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2])
})

# Numbers of self-pairs
FS.SPs.btn = reactive(find.SPs(FS.atts(), k()))

# Numbers of parent-offspring pairs
FS.POPs.btn = reactive(find.POPs.btn(FS.atts(), k()))

# Show table
output$firstNsKPsBtnPop = renderTable({
  SPs.prnts.kwn = find.SPs.prnts.kwn(FS.atts(), k())
  SMPs = find.SMSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  SFPs = find.SFSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  FSPs = find.FSSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  frmt.tbl(
    rbind(
      FS.APs.btn(), FS.SPs.btn(), SPs.prnts.kwn, FS.POPs.btn(),
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

