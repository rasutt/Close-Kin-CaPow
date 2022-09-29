library(shiny)

# Define UI for app
ui <- fluidPage(
  # App title and description
  includeMarkdown("README.md"),
  
  navbarPage(
    "Close-kin CaPow!",
    id = "nav.tab",
    selected = "check.tab",
    # Sim tab ----
    tabPanel(
      title = "Simulate studies",
      value = "sim.tab",
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
            inputId = "n.sims", label = "Number of studies to simulate:",
            0, 1000, value = 10, step = 10
          ),
          actionButton(
            inputId = "simulate", label = "Simulate studies"
          ),
          # helpText("~3 seconds per 100 populations"),
          checkboxInput(
            inputId = "popan", label = "Popan model", value = T
          ),
          checkboxInput(
            inputId = "close.kin", label = "Close-kin model", value = T
          )
        ),
        # Outputs ---- 
        mainPanel(
          h2("Last simulation"),
          h3("Parameter values"),
          tableOutput(outputId = "lastParVals"),
          h3("Simulation options"),
          tableOutput(outputId = "lastSimOpts"),
          plotOutput(outputId = "lastExpPop"),
          
          h2("Next simulation"),
          h3("Parameter values"),
          tableOutput(outputId = "nextParVals"),
          h3("Simulation options"),
          tableOutput(outputId = "nextSimOpts"),
          plotOutput(outputId = "nextExpPop")
        ),
        # ----
      )
    ),
    # Check tab ----
    tabPanel(
      title = "Check simulation",
      value = "check.tab",
      tabsetPanel(
        id = "check.sub.tabs",
        selected = "populations",
        # First study ----
        tabPanel(
          title = "First study",
          value = "first.study",
          h2("First study simulated"),
          # h3("First life-histories"),
          # tableOutput(outputId = "firstLifeHists"),
          
          h3("First sample-histories"),
          tableOutput(outputId = "firstSampHists"),
          
          h3("Numbers of kin-pairs in whole population"),
          h4("Within surveys"),
          tableOutput(outputId = "firstNsKPsPopWtn"),
          h4("Between surveys"),
          tableOutput(outputId = "firstNsKPsPopBtn"),
          
          h3("Estimated numbers of kin-pairs in whole population"),
          h4("Within surveys"),
          tableOutput(outputId = "firstEstNsKPsPopWtn"),
          h4("Between surveys"),
          tableOutput(outputId = "firstEstNsKPsPopBtn")
        ),
        # Populations ----
        tabPanel(
          title = "Population sizes",
          value = "populations",
          h2("Population sizes over time"),
          plotOutput(outputId = "checkExpPop"),
          h2("Population sizes in survey-years"),
          h3("Biases"),
          tableOutput(outputId = "biasNsPopWtn"),
          h3("Error distributions"),
          plotOutput(outputId = "nsWtnPop"),
        ),
        # All pairs ----
        tabPanel(
          title = "All pairs",
          value = "all.pairs",
          h2("Total numbers of pairs of individuals"),
          h3("Biases"),
          tableOutput(outputId = "biasAPsPopWtn"),
          tableOutput(outputId = "biasAPsPopBtn"),
          h3("Error distributions"),
          plotOutput(outputId = "nsAPsWtnPop"),
          plotOutput(outputId = "nsAPsBtnPop")
        ),
        # Parent-offspring pairs ----
        tabPanel(
          title = "Parent-offspring pairs",
          value = "POPs",
          h2("Kin-pair estimator biases"),
          h3("Numbers in whole population"),
          tableOutput(outputId = "biasNsPOPsPopWtn")
        ),
        # Temporal estimates ----
        tabPanel(
          title = "Temporal estimates",
          value = "temp.ests",
          h2("Temporal estimates vs observed averages"),
          textOutput("tempEstNote"),
          h4("Same-mother/father pairs in whole population including animals
          born in each year in the population history"),
          tableOutput(outputId = "nsKPsTemp")
        ),
        # Kin-pairs ----
        tabPanel(
          title = "Bias",
          value = "bias",
          h2("Estimator biases"),
          h3("Numbers in whole population"),
          h4("Temporal estimates"),
          textOutput("tempEstBiasNote"),
          tableOutput(outputId = "biasNsKPsTemp"),
          h4("Within surveys"),
          tableOutput(outputId = "biasNsKPsPopWtn"),
          h4("Between surveys"),
          tableOutput(outputId = "biasNsKPsPopBtn"),
          
          h3("Probabilities"),
          h4("Within surveys"),
          tableOutput(outputId = "biasProbsKPsWtn"),
          h4("Between surveys"),
          tableOutput(outputId = "biasProbsKPsBtn"),
          
          h3("Numbers among sampled animals"),
          textOutput(outputId = "percUnknPrnts"),
          h4("Within surveys"),
          # tableOutput(outputId = "biasNsKPsCapWtn"),
          h4("Between surveys")
          # tableOutput(outputId = "biasNsKPsCapBtn")
        ),
        # Error distributions ----
        tabPanel(
          title = "Error distributions",
          value = "err.dists",
          h2("Kin-pair estimator error distributions"),
          h3("Numbers in whole population"),
          h4("Temporal estimates"),
          plotOutput(outputId = "nsSMPsFnlB1Fnl"),
          plotOutput(outputId = "nsSFPsFnlB1Fnl"),
          plotOutput(outputId = "nsSFPsFnlB"),
          h4("Within surveys"),
          # plotOutput(outputId = "nsSMPsWtnPop"),
          # plotOutput(outputId = "nsSFPsWtnPop"),
          h4("Between surveys"),
          # plotOutput(outputId = "nsSPsBtnPop"),
          plotOutput(outputId = "nsSMPsBtnPop"),
          
          h3("Probabilities"),
          h4("Within surveys"),
          plotOutput(outputId = "probSMPsWtnPop"),
          # h4("Between surveys"),
          # plotOutput(outputId = "probSPsBtnPop"),

          h3("Numbers among sampled animals"),
          h4("Within surveys"),
          plotOutput(outputId = "nsPOPsWtn")
          # plotOutput(outputId = "nsSMPsWtn"),
          # plotOutput(outputId = "nsHSPsWtn"),
          # h4("Between surveys"),
          # plotOutput(outputId = "nsSPsBtn"),
          # plotOutput(outputId = "nsPOPsBtn")
        )
      )
    ), 
    # Model tab ----
    tabPanel(
      title = "Analyze model performance",
      value = "model.tab",
      h2("Analyze model performance"),
      h3("Parameter values"),
      tableOutput(outputId = "modParVals"),
      h3("Simulation options"),
      tableOutput(outputId = "modSimOpts"),
      h3("Model fitting success rates"),
      tableOutput(outputId = "modStats"),
      h3("Estimates"),
      plotOutput(outputId = "modComp"),
      h3("95% confidence interval coverage"),
      tableOutput(outputId = "CICov"),
      h3("Results for first study"),
      tableOutput(outputId = "firstResults"),
      h3("95% confidence intervals for lambda"),
      plotOutput(outputId = "CIPlot")
    ),
    # Save/load tab ----
    tabPanel(
      title = "Save/load",
      value = "save_load",
      h2("Save results"),
      downloadButton("downloadData", "Download"),
      h2("Load results"),
      fileInput("file", "Choose File", accept = ".Rdata"),
      textOutput("upMsg")
    )
    # ----
  )
)
