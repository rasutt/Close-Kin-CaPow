# Outputs for first study sub-tab of checks tab

# Get first study from list
fst.std = reactive(sim.lst()$hists.lst[[1]])

### First study simulated

# Function to format table of integers
frmt.tbl = function(data, rw.nms, cl.nms) {
  mode(data) = "integer"
  df = data.frame(matrix(data, ncol = length(cl.nms)), row.names = rw.nms)
  names(df) = cl.nms
  df
}

## Numbers of kin-pairs in whole population
# Within surveys
output$firstNsKPsWtnPop = renderTable({
  pop.atts = attributes(fst.std())
  SMPs = find.SMPs.wtn(pop.atts, k())
  SFPs = find.SFPs.wtn(pop.atts, k())
  FSPs = find.FSPs.wtn(pop.atts, k())
  frmt.tbl(
    rbind(
      pop.atts$N.t.vec[s.yr.inds()],
      choose(pop.atts$N.t.vec[s.yr.inds()], 2),
      find.POPs.wtn(pop.atts, k()),
      SMPs, SFPs, FSPs, SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps.pop.wtn, srvy.yrs()
  )
}, rownames = T)

# Between surveys
output$firstNsKPsBtnPop = renderTable({
  pop.atts = attributes(fst.std())
  SPs.prnts.kwn = find.SPs.prnts.kwn(pop.atts, k())
  SMPs = find.SMSPs.btn(pop.atts, k()) - SPs.prnts.kwn
  SFPs = find.SFSPs.btn(pop.atts, k()) - SPs.prnts.kwn
  FSPs = find.FSSPs.btn(pop.atts, k()) - SPs.prnts.kwn
  frmt.tbl(
    rbind(
      combn(
        pop.atts$N.t.vec[s.yr.inds()], 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2]
      ),
      find.SPs(pop.atts, k()),  SPs.prnts.kwn, find.POPs.btn(pop.atts, k()),
      SMPs, SFPs, FSPs, SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps.pop.btn, srvy.prs()
  )
}, rownames = T)

## Estimated numbers of kin-pairs in whole population
# Within surveys
output$firstEstNsKPsWtnPop = renderTable({
  frmt.tbl(t(est.ns.kps.pop()$wtn), kp.tps.pop.wtn, srvy.yrs())
}, rownames = T)

# Between surveys
output$firstEstNsKPsBtnPop = renderTable({
  frmt.tbl(t(est.ns.kps.pop()$btn), kp.tps.pop.btn, srvy.prs())
}, rownames = T)

