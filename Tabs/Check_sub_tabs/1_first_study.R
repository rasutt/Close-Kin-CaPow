# Outputs for first study sub-tab of checks tab

### First study simulated

# Function to format table of integers
frmt.tbl = function(data, rw.nms, cl.nms) {
  mode(data) = "integer"
  df = data.frame(matrix(data, ncol = length(cl.nms)), row.names = rw.nms)
  names(df) = cl.nms
  df
}

## First sample-histories
output$firstSampHists = renderTable(head(data.frame(sim.lst()$hists.lst[[1]])))

## Numbers of kin-pairs in whole population
# Within surveys
output$firstNsKPsPopWtn = renderTable({
  pop.atts = attributes(sim.lst()$hists.lst[[1]])
  smps = find.SMPs.wtn(pop.atts, k())
  sfps = find.SFPs.wtn(pop.atts, k())
  fsps = find.FSPs.wtn(pop.atts, k())
  frmt.tbl(
    matrix(
      c(
        pop.atts$N.t.vec[s.yr.inds()],
        choose(pop.atts$N.t.vec[s.yr.inds()], 2),
        find.POPs.wtn(pop.atts, k()),
        smps, sfps, fsps,
        smps + sfps - 2 * fsps
      ), ncol = k(), byrow = T
    ), 
    kp.tps.pop.wtn, srvy.yrs()
  )
}, rownames = T)

# Between surveys
output$firstNsKPsPopBtn = renderTable({
  pop.atts = attributes(sim.lst()$hists.lst[[1]])
  frmt.tbl(
    matrix(
      c(
        combn(
          pop.atts$N.t.vec[s.yr.inds()], 2, 
          function(N.s.pr) N.s.pr[1] * N.s.pr[2]
        ),
        find.SPs(pop.atts, k()),
        find.POPs.btn(pop.atts, k()),
        find.SMPs.btn(pop.atts, k())
      ), ncol = n.srvy.prs(), byrow = T
    ), 
    kp.tps.pop.btn, srvy.prs()
  )
}, rownames = T)

## Estimated numbers of kin-pairs in whole population
# Within surveys
output$firstEstNsKPsPopWtn = renderTable({
  frmt.tbl(t(est.ns.kps.pop.lst()$wtn), kp.tps.pop.wtn, srvy.yrs())
}, rownames = T)

# Between surveys
output$firstEstNsKPsPopBtn = renderTable({
  frmt.tbl(t(est.ns.kps.pop.lst()$btn), kp.tps.pop.btn, srvy.prs())
}, rownames = T)
