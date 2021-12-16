library(shiny)

# Define UI for app
ui <- fluidPage(
  # App title and description
  titlePanel("Close-Kin CaPow!"),
  p("This is a web app to perform close-kin capture-recapture study power 
  analysis via simulation. It simulates multiple populations of animals over 
  time, and datasets from capture-recapture studies of those populations, checks
  that the simulations have the expected characteristics, fits close-kin and/or 
  conventional capture-recapture models to the datasets, and analyses the 
  resulting estimates."),

  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      sliderInput(
        inputId = "rho", label = HTML("Birth rate (&rho;):"),
        min_rho, max_rho, value = 0.08, step = step_rho
      ),
      sliderInput(
        inputId = "phi", label = HTML("Survival rate (&phi;):"),
        min_phi, max_phi, value = 0.95, step = step_phi
      ),
      sliderInput(
        inputId = "hist.len", label = "Length of simulations:",
        10, 100, value = 80, step = 1
      ),
      textInput(
        inputId = "srvy.yrs", label = "Survey years (ascending):",
        value = "1995:1998, 2006:2009, 2020"
      ),
      sliderInput(
        inputId = "n_sims", label = "Number of populations:",
        10, 1000, value = 10, step = 10
      ),
      actionButton(
        inputId = "simulate", label = "Simulate studies"
      ),
      helpText("~1 second per 10 populations"), # TMB locally
      checkboxGroupInput(
        inputId = "models", label = "Fit models:",
        choices = c("POPAN", "Close kin"), selected = "Close kin", inline = T
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
      textOutput(outputId = "percUnknPrnts"),
      plotOutput(outputId = "NLLPlot"),
      h4("Parameter estimates for first few studies"),
      tableOutput(outputId = "firstEsts"),
      h4("Parameter estimates for all studies"),
      plotOutput(outputId = "modComp"),
      # plotOutput(outputId = "CIPlot")
      # textOutput(outputId = "ciCov")
    )
  )
)
