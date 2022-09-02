# Save/load tab

# Save checks
output$downloadData <- downloadHandler(
  filename = "ckc_saved_objs.Rdata",
  content = function(file) {
    saved_objs = list(
      phi = phi(),
      rho = rho(),
      lambda = lambda(),
      base.yr = base.yr(),
      exp.N.base = exp.N.base(),
      srvy.yrs = srvy.yrs(),
      p = p(),
      clvng.ints = clvng.ints(),
      clvng.p = clvng.p(),
      tmp.emgn = tmp.emgn(),
      alpha = alpha(),
      hist.len = hist.len(),
      n.sims = n.sims(),
      srvy.gaps = srvy.gaps(),
      k = k(),
      n.srvy.prs = n.srvy.prs(),
      f.year = f.year(),
      s.yr.inds = s.yr.inds(),
      exp.N.t = exp.N.t(),
      exp.N.fin = exp.N.fin(),
      exp.Ns = exp.Ns(),
      par.names = par.names(),
      par.vals = par.vals(),
      sim.opts= sim.opts(),
      sim.lst = sim.lst(),
      checks.lst = checks.lst()
    )
    save(saved_objs, file = file)
  }
)

# Load checks
observeEvent(input$file, {
  load(input$file$datapath)
  checks.lst(checks)
})

# Show number of files uploaded
up.msg = reactiveVal("No files uploaded")
observeEvent(input$file, {
  up.msg(paste0("File uploaded!"))
})
output$upMsg = renderText(up.msg())
