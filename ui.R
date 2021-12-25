library(shiny)

# Define UI for app
ui <- fluidPage(
  # App title and description
  includeMarkdown("README.md"),
  
  sidebarLayout(
    # Inputs ----
    sidebarPanel(
      sliderInput(
        inputId = "phi", label = HTML("Survival rate (&phi;):"),
        min_phi, max_phi, value = 0.95, step = step_phi
      ),
      sliderInput(
        inputId = "rho", label = HTML("Birth rate (&rho;):"),
        min_rho, max_rho, value = 0.08, step = step_rho
      ),
      sliderInput(
        inputId = "exp.N.base", 
        label = HTML("Expected population size in base year:"),
        0, 3000, value = 1500, step = 50
      ),
      sliderInput(
        inputId = "base.yr", 
        label = HTML("Base year for expected population size:"),
        1990, 2030, value = 2009, step = 1
      ),
      textInput(
        inputId = "srvy.yrs", label = "Survey years:",
        value = "1995:1998, 2006:2009, 2020"
      ),
      sliderInput(
        inputId = "p", label = "Capture probability:",
        0, 0.2, value = 0.1, step = 0.01
      ),
      sliderInput(
        inputId = "hist.len", label = "Length of population histories:",
        0, 100, value = 80, step = 10
      ),
      sliderInput(
        inputId = "n_sims", label = "Number of studies to simulate:",
        0, 1000, value = 10, step = 10
      ),
      actionButton(
        inputId = "simulate", label = "Simulate studies"
      ),
      helpText("~1 second per 10 populations"),
      checkboxInput(
        inputId = "popan", label = "Popan model", value = T
      ),
      checkboxInput(
        inputId = "close_kin", label = "Close-kin model", value = T
      )
    ),
    # ----
    mainPanel(
      tabsetPanel(
        id = "tab",
        selected = "sim_tab",
        # Sim tab ----
        tabPanel(
          title = "Simulate studies",
          value = "sim_tab",
          h4("Values simulated"),
          tableOutput(outputId = "simParVals"),
          # plotOutput(outputId = "simExpPop"),
          h4("Values to simulate next"),
          tableOutput(outputId = "selParVals"),
          plotOutput(outputId = "selExpPop")
        ),
        # ----
        # Check tab ----
        tabPanel(
          title = "Check simulation",
          value = "check_tab",
          h4("Values simulated"),
          tableOutput(outputId = "checkParVals"),
          h4("Simulation parameters"),
          tableOutput(outputId = "checkSimVals"),
          plotOutput(outputId = "popPlot"),
          plotOutput(outputId = "nPOPsPlot"),
          textOutput(outputId = "percUnknPrnts"),
          h4("Head of first dataset"),
          tableOutput(outputId = "dataHead")
        ), 
        # ----
        # Model tab ----
        tabPanel(
          title = "Analyze model performance",
          value = "model_tab",
          h4("Values simulated"),
          tableOutput(outputId = "modParVals"),
          h4("Simulation parameters"),
          tableOutput(outputId = "modSimVals"),
          h4("Model fitting success rates"),
          tableOutput(outputId = "modStats"),
          h4("Estimates"),
          plotOutput(outputId = "modComp"),
          h4("95% confidence interval coverage"),
          tableOutput(outputId = "CICov"),
          h4("Results for first study"),
          tableOutput(outputId = "firstResults"),
          h4("95% confidence intervals for lambda"),
          plotOutput(outputId = "CIPlot")
        )
        # ----
      )
    )
  )
)