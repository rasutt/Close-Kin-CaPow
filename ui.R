library(shiny)

# Define UI for app
ui <- fluidPage(
  # App title
  titlePanel("Close Kin CaPow"),
  
  p("This is a little web app in Shiny to illustrate close-kin capture-recapture
    parameter estimation via simulation. It simulates datasets from capture-
    recapture studies of a population of animals over time.  It takes the 
    simulated population growth rate as input, and shows the population sizes 
    over time, the negative log-likelihood for the first sample, and the 
    expected versus observed numbers of kin-pairs and the MLEs
    over all samples."),

  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      sliderInput(
        inputId = "lambda", 
        label = "Population growth rate (lambda):",
        min_lambda, 
        max_lambda, 
        value = 1.03, 
        step = step_lambda
      ),
      sliderInput(
        inputId = "n_sims", 
        label = "Number of populations:",
        10, 
        1000, 
        value = 10, 
        step = 10
      ),
      helpText("~3 seconds per 10 populations") # shinapps.io
      # helpText("~1 second per 10 populations") # TMB locally
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      h4("Head of first sample"),
      plotOutput(outputId = "popPlot"),
      tableOutput(outputId = "dataHead"),
      plotOutput(outputId = "nPOPsPlot"),
      textOutput(outputId = "nUnknPrnts"),
      plotOutput(outputId = "NLLPlot"),
      h4("Parameter estimates for first few studies"),
      tableOutput(outputId = "firstEsts"),
      h4("Parameter estimates for all studies"),
      plotOutput(outputId = "modComp")
      # plotOutput(outputId = "scatterPlot"),
      # plotOutput(outputId = "MLEplot"),
      # plotOutput(outputId = "CIPlot"),
      # textOutput(outputId = "ciCov")
    )
  )
)
