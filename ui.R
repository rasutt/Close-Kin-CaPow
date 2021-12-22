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
        50, 3050, value = 1500, step = 50
      ),
      sliderInput(
        inputId = "base.yr", 
        label = HTML("Base year for expected population size:"),
        2002, 2022, value = 2009, step = 1
      ),
      textInput(
        inputId = "srvy.yrs", label = "Survey years:",
        value = "1995:1998, 2006:2009, 2020"
      ),
      sliderInput(
        inputId = "hist.len", label = "Length of population histories:",
        40, 100, value = 80, step = 10
      ),
      sliderInput(
        inputId = "n_sims", label = "Number of studies to simulate:",
        10, 1010, value = 10, step = 10
      ),
      actionButton(
        inputId = "simulate", label = "Simulate studies"
      ),
      helpText("~1 second per 10 populations"),
      checkboxGroupInput(
        inputId = "models", label = "Fit models:",
        choices = mod_choices, selected = "Close-kin", inline = T
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
          h4("Implied parameters"),
          textOutput(outputId = "lambda"),
          textOutput(outputId = "exp.Ns"),
          plotOutput(outputId = "expPopPlot")
        ),
        # ----
        # Check tab ----
        tabPanel(
          title = "Check simulation",
          value = "check_tab",
          plotOutput(outputId = "popPlot"),
          h4("Head of first dataset"),
          tableOutput(outputId = "dataHead"),
          plotOutput(outputId = "nPOPsPlot"),
          textOutput(outputId = "percUnknPrnts")
        ), 
        # ----
        # Model tab ----
        tabPanel(
          title = "Analyze model performance",
          value = "model_tab",
          h4("Parameter estimates from all models for all studies"),
          plotOutput(outputId = "modComp"),
          h4("Model fitting success rates"),
          tableOutput(outputId = "modStats"),
          h4(
            "95% confidence interval coverage 
            (log-normal for population parameters)"
          ),
          tableOutput(outputId = "CICov"),
          h4("95% confidence intervals for lambda"),
          plotOutput(outputId = "CIPlot"),
          plotOutput(outputId = "NLLPlot"),
          h4("Parameter estimates from close-kin model for first study"),
          tableOutput(outputId = "firstEsts")
        )
        # ----
      )
    )
  )
)