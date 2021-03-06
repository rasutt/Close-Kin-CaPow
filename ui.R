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
        # value = "1995:1998, 2006:2009, 2020"
        value = "2000, 2010, 2020"
      ),
      sliderInput(
        inputId = "p", label = "Base level capture probability:",
        0, 0.2, value = 0.1, step = 0.01
      ),
      checkboxInput(
        inputId = "clvng.ints", 
        label = "Females breed in order of time since last breeding", 
        value = F
      ),
      sliderInput(
        inputId = "clvng.p", 
        label = "Additional capture probability when calving:",
        0, 0.5, value = 0, step = 0.05
      ),
      sliderInput(
        inputId = "tmp.emgn", 
        label = "Probability that males absent:",
        0, 1, value = 0, step = 0.05
      ),
      sliderInput(
        inputId = "alpha", label = "Age of sexual maturity:",
        1, 15, value = 8, step = 1
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
      helpText("~3 seconds per 100 populations"),
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
        selected = "check_tab",
        # Sim tab ----
        tabPanel(
          title = "Simulate studies",
          value = "sim_tab",
          h4("Last simulation"),
          tableOutput(outputId = "simParVals"),
          tableOutput(outputId = "lastSimVals"),
          plotOutput(outputId = "simExpPop"),
          h4("Next simulation"),
          tableOutput(outputId = "selParVals"),
          tableOutput(outputId = "nextSimVals"),
          plotOutput(outputId = "selExpPop")
        ),
        # ----
        # Check tab ----
        tabPanel(
          title = "Check simulation",
          value = "check_tab",
          h4("Values simulated"),
          tableOutput(outputId = "checkParVals"),
          tableOutput(outputId = "checkSimVals"),
          
          h4("Values observed"),
          plotOutput(outputId = "popPlot"),
          tableOutput(outputId = "slfPrPrb"),
          textOutput(outputId = "percUnknPrnts"),
          
          h2("Kin-pair estimator bias"),
          h3("Numbers among population"),
          h4("Within survey years"),
          tableOutput(outputId = "nsKPsPopWtn"),
          h4("Between survey years"),
          tableOutput(outputId = "nsKPsPopBtn"),
          
          h3("Probabilies"),
          tableOutput(outputId = "nsKPsProb"),
          
          h3("Numbers among sampled animals"),
          h4("Within survey years"),
          tableOutput(outputId = "nsKPsCapWtn"),
          h4("Between survey years"),
          tableOutput(outputId = "nsKPsCapBtn"),
          
          h4("Same-mother/father pairs in whole population including animals 
          born in each year in the population history"),
          tableOutput(outputId = "nsKPsT"),
          
          h4("Kin-pair estimate errors - Whole population"),
          # plotOutput(outputId = "nsAPsWtnPop"),
          # plotOutput(outputId = "nsAPsBtnPop"),
          # plotOutput(outputId = "nsSPsBtnPop"),
          # plotOutput(outputId = "nsSMPsWtnPop"),
          # plotOutput(outputId = "nsSFPsWtnPop"),
          # plotOutput(outputId = "nsSMPsBtnPop"),
          h4("Kin-pair probability errors"),
          # plotOutput(outputId = "probSPsBtnPop"),
          # plotOutput(outputId = "probSMPsWtnPop"),
          h4("Kin-pair estimate errors - Captured animals"),
          # plotOutput(outputId = "nsSPsBtn"),
          # plotOutput(outputId = "nsPOPsWtn"),
          # plotOutput(outputId = "nsPOPsBtn"),
          # plotOutput(outputId = "nsSMPsWtn"),
          # plotOutput(outputId = "nsHSPsWtn"),
          h4("Head of first life-histories table"),
          tableOutput(outputId = "alive"),
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