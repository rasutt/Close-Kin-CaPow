# Load shiny package
library(shiny)

# Function to make table outputs for unknown parents modules
unknPrntsTbl <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Unknown parents"),
    p("Average percentage of individuals with unknown parents, from start
    of simulation."),
    tableOutput(outputId = ns("unknPrntsWtn")),
    tableOutput(outputId = ns("unknPrntsBtn"))
  )
}

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
            value = "2000, 2020"
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
        selected = "POPs.tab",
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
          h2("Population sizes"),
          p("Numbers of individuals that are alive in the population."),
          plotOutput(outputId = "checkExpPop"),
          h3("Biases"),
          tableOutput(outputId = "biasNsPopWtn"),
          h3("Error distributions"),
          plotOutput(outputId = "nsWtnPop"),
        ),
        # All pairs ----
        tabPanel(
          title = "All pairs",
          value = "all.pairs",
          h2("All pairs"),
          p("Total numbers of pairs of individuals."),
          h3("Biases"),
          tableOutput(outputId = "biasAPsPopWtn"),
          tableOutput(outputId = "biasAPsPopBtn"),
          h3("Error distributions"),
          plotOutput(outputId = "nsAPsWtnPop"),
          plotOutput(outputId = "nsAPsBtnPop")
        ),
        # Self-pairs ----
        tabPanel(
          title = "Self-pairs",
          value = "SPs.tab",
          h2("Self-pairs"),
          p("Numbers of pairs of individuals that are the same individual in
          different survey-years."),
          h3("Biases"),
          tableOutput(outputId = "biasSPsPop"),
          h3("Error distributions"),
          plotOutput(outputId = "nsSPsPop"),
        ),
        # Parent-offspring pairs ----
        tabPanel(
          title = "Parent-offspring pairs",
          value = "POPs.tab",
          h2("Parent-offspring pairs"),
          p("Numbers of pairs of individuals that are parent and offspring."),
          unknPrntsTbl("POPs.tab"),
          h3("Biases"),
          tableOutput(outputId = "biasPOPsPopWtn"),
          tableOutput(outputId = "biasPOPsPopBtn"),
          h3("Error distributions"),
          plotOutput(outputId = "nsPOPsWtnPop"),
          plotOutput(outputId = "nsPOPsBtnPop")
        ),
        # Same-mother pairs ----
        tabPanel(
          title = "Same-mother pairs",
          value = "SMPs.tab",
          h2("Same-mother pairs"),
          p("Numbers of pairs of individuals with the same mothers."),
          unknPrntsTbl("SMPs.tab"),
          h3("Biases"),
          p("Numbers in the final year, with one born in
            the year indicated, and one born in the final year."),
          tableOutput(outputId = "biasSMPsAgeKnwn"),
          p("Numbers in and between survey-years, ages unknown."),
          tableOutput(outputId = "biasSMPsPopWtn"),
          tableOutput(outputId = "biasSMPsPopBtn"),
          h3("Error distributions"),
          plotOutput(outputId = "nsSMPsAgeKnwn"),
          plotOutput(outputId = "nsSMPsWtnPop"),
          plotOutput(outputId = "nsSMPsBtnPop")
        ),
        # Same-father pairs ----
        tabPanel(
          title = "Same-father pairs",
          value = "SFPs.tab",
          h2("Same-father pairs"),
          p("Numbers of pairs of individuals with the same fathers."),
          unknPrntsTbl("SFPs.tab"),
          h3("Biases"),
          p("Numbers in the final year, with one born in
            the year indicated, and one born in the final year."),
          tableOutput(outputId = "biasSFPsAgeKnwn"),
          p("Numbers in the final year, both born in the year indicated."),
          tableOutput(outputId = "biasSFPsSameAge"),
          p("Numbers in and between survey-years, ages unknown."),
          tableOutput(outputId = "biasSFPsPopWtn"),
          # tableOutput(outputId = "biasSFPsPopBtn"),
          h3("Error distributions"),
          plotOutput(outputId = "nsSFPsAgeKnwn"),
          plotOutput(outputId = "nsSFPsSameAge"),
          plotOutput(outputId = "nsSFPsWtnPop"),
          # plotOutput(outputId = "nsSFPsBtnPop")
        ),
        # Sibling-pairs ----
        tabPanel(
          title = "Sibling-pairs",
          value = "SiPs.tab",
          h2("Sibling-pairs"),
          p("Pairs of individuals that share one or two parents."),
          unknPrntsTbl("SibPs.tab"),
          h3("Biases"),
          h4("Full-sibling pairs"),
          p("Pairs of individuals that share two parents."),
          tableOutput(outputId = "biasFSPsPopWtn"),
          # tableOutput(outputId = "biasFSPsPopBtn"),
          h4("Half-sibling pairs"),
          p("Pairs of individuals that share one parent."),
          tableOutput(outputId = "biasHSPsPopWtn"),
          # tableOutput(outputId = "biasSFPsPopBtn"),
          h3("Error distributions"),
          plotOutput(outputId = "nsFSPsWtnPop"),
          plotOutput(outputId = "nsHSPsWtnPop")
          # plotOutput(outputId = "nsSFPsWtnPop"),
          # plotOutput(outputId = "nsSFPsBtnPop")
        ),
        # Biases ----
        tabPanel(
          title = "Overall biases",
          value = "bias.tab",
          h2("Overall biases"),
          p("Average estimator biases over all years."),
          h3("Unknown parents"),
          p("Average percentage of individuals with unknown parents, from start
            of simulation, in survey-years, and pairs of survey-years."),
          tableOutput(outputId = "percUnknPrnts"),
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
          
          # h3("Numbers among sampled animals"),
          # h4("Within surveys"),
          # tableOutput(outputId = "biasNsKPsCapWtn"),
          # h4("Between surveys")
          # tableOutput(outputId = "biasNsKPsCapBtn")
        ),
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
