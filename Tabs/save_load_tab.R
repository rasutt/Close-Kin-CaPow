# Save/load tab

# Save checks
output$downloadData <- downloadHandler(
  filename = "checks.Rdata",
  content = function(file) {
    checks = checks.lst()
    save(checks, file = file)
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
