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
        inputId = "rho", 
        label = HTML("Birth rate (&rho;):"),
        min_rho, 
        max_rho, 
        value = 0.08, 
        step = step_rho
      ),
      sliderInput(
        inputId = "pi", 
        label = HTML("Mortality rate (&pi;):"),
        min_pi, 
        max_pi, 
        value = 0.05, 
        step = step_pi
      ),
      # checkboxGroupInput(
      #   inputId = "srvy.yrs",
      #   label = "Survey years:",
      #   choices = as.character(1995:2020),
      #   selected = as.character(c(1995:1998, 2006:2009, 2020)),
      #   inline = T
      # ),
      textInput(
        inputId = "srvy.yrs",
        label = "Survey years (ascending):",
        value = paste(c(1995:1998, 2006:2009, 2020), collapse = ", ")
      ),
      sliderInput(
        inputId = "n_sims", 
        label = "Number of populations:",
        10, 
        1000, 
        value = 40, 
        step = 10
      ),
      actionButton(
        inputId = "simulate",
        label = "Simulate studies"
      ),
      # helpText("~3 seconds per 10 populations"), # shinapps.io
      helpText("~1 second per 10 populations"), # TMB locally
      checkboxGroupInput(
        inputId = "models",
        label = "Fit models:",
        choices = c("POPAN", "Close kin"),
        selected = "Close kin",
        inline = T
      )
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      h4("Implied parameters"),
      textOutput(outputId = "phi"),
      textOutput(outputId = "lambda"),
      textOutput(outputId = "k"),
      textOutput(outputId = "f.year"),
      plotOutput(outputId = "popPlot"),
      h4("Head of first dataset"),
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
