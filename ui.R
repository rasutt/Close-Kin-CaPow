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
        selected = "demo.tab",
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
          tableOutput(outputId = "firstNsKPsWtnPop"),
          h4("Between surveys"),
          tableOutput(outputId = "firstNsKPsBtnPop"),
          
          h3("Estimated numbers of kin-pairs in whole population"),
          h4("Within surveys"),
          tableOutput(outputId = "firstEstNsKPsWtnPop"),
          h4("Between surveys"),
          tableOutput(outputId = "firstEstNsKPsBtnPop")
        ),
        # Populations ----
        tabPanel(
          title = "Population sizes",
          value = "populations",
          h2("Population sizes"),
          p("Numbers of individuals that are alive in the population."),
          plotOutput(outputId = "checkExpPop"),
          h3("Biases"),
          tableOutput(outputId = "bsNsWtnPop"),
          h3("Error distributions"),
          plotOutput(outputId = "errsNsWtnPop"),
        ),
        # Demographics ----
        tabPanel(
          title = "Demographics",
          value = "demo.tab",
          h2("Demographics"),
          p("Individual survival and population growth rates."),
          h3("Biases"),
          h4("Survival rate"),
          tableOutput(outputId = "bsPhi"),
          h3("Error distributions"),
          plotOutput(outputId = "errsPhi"),
        ),
        # All pairs ----
        tabPanel(
          title = "All pairs",
          value = "all.pairs",
          h2("All pairs"),
          p("Total numbers of pairs of individuals."),
          h3("Biases"),
          tableOutput(outputId = "bsAPsWtnPop"),
          tableOutput(outputId = "bsAPsBtnPop"),
          h3("Error distributions"),
          plotOutput(outputId = "errsAPsWtnPop"),
          plotOutput(outputId = "errsAPsBtnPop")
        ),
        # Self-pairs ----
        tabPanel(
          title = "Self-pairs",
          value = "SPs.tab",
          h2("Self-pairs"),
          p("Numbers of pairs of individuals that are the same individual in
          different survey-years."),
          h3("Biases"),
          tableOutput(outputId = "bsSPsPop"),
          h3("Error distributions"),
          plotOutput(outputId = "errsSPsPop")
        ),
        # Parent-offspring pairs ----
        tabPanel(
          title = "Parent-offspring pairs",
          value = "POPs.tab",
          h2("Parent-offspring pairs"),
          p("Numbers of pairs of individuals that are parent and offspring."),
          unknPrntsTbl("POPs.tab"),
          h3("Biases"),
          tableOutput(outputId = "bsPOPsWtnPop"),
          tableOutput(outputId = "bsPOPsBtnPop"),
          h3("Error distributions"),
          plotOutput(outputId = "errsPOPsWtnPop"),
          plotOutput(outputId = "errsPOPsBtnPop")
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
          tableOutput(outputId = "bsSMPsAgeKnwn"),
          p("Numbers in and between survey-years, ages unknown."),
          tableOutput(outputId = "bsSMPsWtnPop"),
          tableOutput(outputId = "bsSMPsBtnPop"),
          p("Numbers between survey-years, one born five years before first, one
            born in last."),
          tableOutput(outputId = "bsSMPsBtnAgeKnwnPop"),
          h3("Error distributions"),
          plotOutput(outputId = "errsSMPsAgeKnwn"),
          plotOutput(outputId = "errsSMPsWtnPop"),
          plotOutput(outputId = "errsSMPsBtnPop"),
          plotOutput(outputId = "errsSMPsBtnAgeKnwnPop")
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
          tableOutput(outputId = "bsSFPsAgeKnwn"),
          p("Numbers in the final year, both born in the year indicated."),
          tableOutput(outputId = "bsSFPsSameAge"),
          p("Numbers in survey-years, ages unknown."),
          tableOutput(outputId = "bsSFPsWtnPop"),
          # tableOutput(outputId = "bsSFPsBtnPop"),
          h3("Error distributions"),
          plotOutput(outputId = "errsSFPsAgeKnwn"),
          plotOutput(outputId = "errsSFPsSameAge"),
          plotOutput(outputId = "errsSFPsWtnPop"),
          # plotOutput(outputId = "errsSFPsBtnPop")
        ),
        # Sibling-pairs ----
        tabPanel(
          title = "Sibling-pairs",
          value = "SibPs.tab",
          h2("Sibling-pairs"),
          p("Pairs of individuals that share one or two parents."),
          unknPrntsTbl("SibPs.tab"),
          h3("Biases"),
          h4("Full-sibling pairs"),
          p("Pairs of individuals that share two parents."),
          tableOutput(outputId = "bsFSPsWtnPop"),
          # tableOutput(outputId = "bsFSPsBtnPop"),
          h4("Half-sibling pairs"),
          p("Pairs of individuals that share one parent."),
          tableOutput(outputId = "bsHSPsWtnPop"),
          # tableOutput(outputId = "bsSFPsBtnPop"),
          h3("Error distributions"),
          plotOutput(outputId = "errsFSPsWtnPop"),
          plotOutput(outputId = "errsHSPsWtnPop")
          # plotOutput(outputId = "errsSFPsWtnPop"),
          # plotOutput(outputId = "errsSFPsBtnPop")
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
          tableOutput(outputId = "bsNsKPsTemp"),
          h4("Within surveys"),
          tableOutput(outputId = "bsNsKPsWtnPop"),
          h4("Between surveys"),
          tableOutput(outputId = "bsNsKPsBtnPop"),
          
          # h3("Probabilities"),
          # h4("Within surveys"),
          # tableOutput(outputId = "bsProbsKPsWtn"),
          # h4("Between surveys"),
          # tableOutput(outputId = "bsProbsKPsBtn"),
          
          # h3("Numbers among sampled animals"),
          # h4("Within surveys"),
          # tableOutput(outputId = "bsNsKPsCapWtn"),
          # h4("Between surveys")
          # tableOutput(outputId = "bsNsKPsCapBtn")
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
