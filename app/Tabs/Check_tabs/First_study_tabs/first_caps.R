# Get first study from list
frst.std = reactive(sim.lst()$hists.lst[[1]])

# Summarise data for POPAN model
frst.pop.sum = reactive({
  pop.sum = as.matrix(FindPopSum(k(), frst.std(), n.cap.hists()))
  mode(pop.sum) = "integer"
  data.frame(pop.sum)
})

# First sample-histories
output$firstSampHists = renderTable(head(data.frame(frst.std())))

# First sample-history sufficient statistics
output$firstSmpHstStats = renderTable({
  df = t(frst.pop.sum())
  rownames(df) = c(
    "Number sampled for first time", "Total number sampled", 
    "Number not sampled, but sampled earlier and later",
    "Number sampled, or sampled earlier, and sampled later",
    "Number sampled for last time"
  )
  df
}, rownames = T)

# First number of sample-histories
output$firstNSmpHsts = renderText({
  paste("Total number of sample histories:", n.cap.hists())
})

