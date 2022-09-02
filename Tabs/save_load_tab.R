# Save/load tab

# Save checks
output$downloadData <- downloadHandler(
  filename = "ckc_saved_objs.Rdata",
  content = function(file) {
    saved_objs = list(
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
