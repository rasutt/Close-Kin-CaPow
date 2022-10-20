# Load shiny package
library(shiny)

# Function to make outputs for values, predictions, and errors
VPE.ui <- function(id, title, desc, types) {
  ns <- NS(id)
  tagList(
    h2(title),
    p(desc),
    h3("Values"),
    lapply(paste0("vals", types), function(v.p) plotOutput(outputId = ns(v.p))),
    h3("Biases"),
    lapply(paste0("bss", types), function(b.t) tableOutput(outputId = ns(b.t))),
    h3("Errors"),
    lapply(paste0("errs", types), function(e.p) plotOutput(outputId = ns(e.p)))
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
        selected = "fst.std.tab",
        # First study ----
        tabPanel(
          title = "First study",
          value = "fst.std.tab",
          h2("First study simulated"),
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
          value = "N.tab",
          plotOutput(outputId = "checkExpPop"),
          VPE.ui(
            "N", "Population sizes", 
            "Numbers of individuals that are alive in the population.", "WtnPop"
          )
        ),
        # Demographics ----
        tabPanel(
          title = "Survival",
          value = "phi.tab",
          VPE.ui("phi", "Survival rates", "Individual survival rates.", "All")
        ),
        # Unknown parents ----
        tabPanel(
          title = "Unknown parents",
          value = "UPs.tab",
          h3("Unknown parents"),
          p("The first generation is simulated with unknown parents. 
          These individuals affect the observed numbers of kin-pairs, as parent
          -offspring pairs among them, and sibling-pairs including them, are 
          both unknown. This causes the appearance of prediction error, so it
          should be taken into account when evaluating predictor performance."),
          p("Below are the average percentages of individuals with unknown
          parents, alive within, and bewteen pairs of, survey-years.  The
          observed numbers of pairs including one or more of these individuals
          will be added soon."),
          tableOutput(outputId = "UPsWtn"),
          tableOutput(outputId = "UPsBtn"),
        ),
        # All pairs ----
        tabPanel(
          title = "All pairs",
          value = "APs.tab",
          VPE.ui(
            "APs", "All pairs", 
            "Total numbers of pairs of individuals.", c("WtnPop", "BtnPop")
          )
        ),
        # Self-pairs ----
        tabPanel(
          title = "Self-pairs",
          value = "SPs.tab",
          VPE.ui(
            "SPs", "Self-pairs", 
            "Numbers of pairs of individuals that are the same individual in
            different survey-years.", c("BtnPop", "PrntsKwn")
          )
        ),
        # Parent-offspring pairs ----
        tabPanel(
          title = "Parent-offspring pairs",
          value = "POPs.tab",
          VPE.ui(
            "POPs", "Parent-offspring pairs", 
            "Numbers of pairs of individuals that are parent and offspring.",
            c("WtnPop", "BtnPop")
          )
        ),
        # Same-mother pairs ----
        tabPanel(
          title = "Same-mother pairs",
          value = "SMPs.tab",
          p("Numbers in the final year, with one born in
            the year indicated, and one born in the final year."),
          p("Numbers in and between survey-years, ages unknown."),
          p("Numbers between survey-years, one born five years before first, one
            born in last."),
          VPE.ui(
            "SMPs", "Same-mother pairs", 
            "Numbers of pairs of individuals with the same mothers.",
            c("AgeKnwn", "WtnPop", "BtnPop", "BtnAgeKnwnPop")
          )
        ),
        # Same-father pairs ----
        tabPanel(
          title = "Same-father pairs",
          value = "SFPs.tab",
          p("Numbers in the final year, with one born in
            the year indicated, and one born in the final year."),
          p("Numbers in the final year, both born in the year indicated."),
          p("Numbers in and between survey-years, ages unknown."),
          VPE.ui(
            "SFPs", "Same-father pairs", 
            "Numbers of pairs of individuals with the same fathers.",
            c("AgeKnwn", "SameAge", "WtnPop", "BtnPop")
          )
        ),
        # Sibling-pairs ----
        tabPanel(
          title = "Sibling-pairs",
          value = "SibPs.tab",
          VPE.ui(
            "FSPs", "Full-sibling pairs",
            "Numbers of pairs of individuals that share two parents.",
            c("WtnPop", "BtnPop")
          ),
          VPE.ui(
            "HSPs", "Half-sibling pairs",
            "Numbers of pairs of individuals that share exactly one parent.",
            c("WtnPop", "BtnPop")
          )
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
