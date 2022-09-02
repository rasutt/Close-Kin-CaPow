# Outputs for first study sub-tab of checks tab

### First study simulated

# Function to format table of integers
form.tab = function(data, rw.nms, cl.nms) {
  mode(data) = "integer"
  df = data.frame(data, row.names = rw.nms)
  names(df) = cl.nms
  df
}

## First life-histories
output$firstLifeHists = renderTable({
  form.tab(head(attributes(sim.lst()$hists.lst[[1]])$alv.mat), NULL, sim.yrs())
})

## First sample-histories
output$firstSampHists = renderTable(head(data.frame(sim.lst()$hists.lst[[1]])))

## Numbers of kin-pairs in whole population
# Within surveys
output$firstNsKPsPopWtn = renderTable({
  form.tab(
    t(checks.lst()$ns.kps.pop.wtn.arr[1, , ]), kp.tps.pop.wtn, srvy.yrs()
  )
}, rownames = T)

# Between surveys
output$firstNsKPsPopBtn = renderTable({
  form.tab(
    t(checks.lst()$ns.kps.pop.btn.arr[1, , ]), kp.tps.pop.btn,
    apply(combn(srvy.yrs(), 2), 2, paste, collapse = "-")
  )
}, rownames = T)

## Estimated numbers of kin-pairs in whole population
# Within surveys
output$firstEstNsKPsPopWtn = renderTable({
  form.tab(
    t(est.ns.kps.pop.lst()$wtn), kp.tps.pop.wtn, srvy.yrs()
  )
}, rownames = T)

# Between surveys
output$firstEstNsKPsPopBtn = renderTable({
  form.tab(
    t(est.ns.kps.pop.lst()$btn), kp.tps.pop.btn, 
    apply(combn(srvy.yrs(), 2), 2, paste, collapse = "-")
  )
}, rownames = T)
