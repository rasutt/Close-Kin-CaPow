library(shiny)

# Define UI for app
ui <- fluidPage(
  # App title and description ----
  includeMarkdown("README.md"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
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
    
    # Main panel for displaying outputs ----
    mainPanel(
      h3("Preparing next simulation"),
      h4("Implied parameters"),
      textOutput(outputId = "lambda"),
      textOutput(outputId = "k"),
      textOutput(outputId = "f.year"),
      h3("Checking last simulation"),
      plotOutput(outputId = "popPlot"),
      h4("Head of first dataset"),
      tableOutput(outputId = "dataHead"),
      plotOutput(outputId = "nPOPsPlot"),
      textOutput(outputId = "percUnknPrnts"),
      h3("Analysing model performance"),
      textOutput(outputId = "cnvgRate"),
      textOutput(outputId = "sesRate"),
      plotOutput(outputId = "NLLPlot"),
      h4("Parameter estimates from first model for first few studies"),
      tableOutput(outputId = "firstEsts"),
      h4("Parameter estimates from all models for all studies"),
      plotOutput(outputId = "modComp"),
      plotOutput(outputId = "CIPlot"),
      textOutput(outputId = "CICov")
    )
    # ----
  )
)
