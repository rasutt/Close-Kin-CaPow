# Load shiny package
library(shiny)

# Function to make outputs for values, biases, and errors
VBE.ui <- function(id, types) {
  ns <- NS(id)
  tagList(
    h3("Values"),
    lapply(paste0("vals", types), function(v.p) plotOutput(outputId = ns(v.p))),
    h3("Biases"),
    lapply(paste0("bss", types), function(b.t) tableOutput(outputId = ns(b.t))),
    h3("Errors"),
    lapply(paste0("errs", types), function(e.p) plotOutput(outputId = ns(e.p)))
  )
}

# Function to make outputs for titles, descriptions, values, biases, and errors
TDVBE.ui <- function(id, title, desc, types) {
  tagList(
    h2(title),
    p(desc),
    VBE.ui(id, types)
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
        selected = "pars.opts.tab",
        # Parameters and options ----
        tabPanel(
          title = "Parameters and options",
          value = "pars.opts.tab",
          h2("Parameters and options"),
          p("Parameter values, simulation options, and biological scenario
            selected for simulation."),
          h3("Parameter values"),
          tableOutput(outputId = "lastParVals"),
          h3("Simulation options"),
          tableOutput(outputId = "lastSimOpts"),
          h3("Biological scenario"),
          tableOutput(outputId = "lastBioScen")
        ),
        # First study ----
        tabPanel(
          title = "First study",
          value = "fst.std.tab",
          h2("First study simulated"),
          p("Some data from the first population and sampling study
            simulated, as an example.  Note that individual simulations can be
            quite different from each other, but the average over many
            simulations should approach the expected values."),
          h3("First sample-histories"),
          p("The data for the first few animals sampled
            in the study, by earliest birth year."),
          tableOutput(outputId = "firstSampHists"),
          h3("Numbers of pairs within survey-years (whole population)"),
          p("Numbers of pairs of individuals with given relationships,
            where both are alive in the given survey-year.  Population
            sizes are included for reference."),
          h4("Predicted"),
          tableOutput(outputId = "firstEstNsKPsWtnPop"),
          h4("Simulated"),
          tableOutput(outputId = "firstNsKPsWtnPop"),
          h3("Numbers of pairs between survey-years (whole population)"),
          p("Numbers of pairs of individuals with given relationships,
            where one individual is alive in each of the given pair of 
            survey-years."),
          h4("Predicted"),
          tableOutput(outputId = "firstEstNsKPsBtnPop"),
          h4("Simulated"),
          tableOutput(outputId = "firstNsKPsBtnPop")
        ),
        # Populations ----
        tabPanel(
          title = "Population sizes",
          value = "N.tab",
          h2("Population sizes"),
          p("Numbers of individuals that are alive in the population in each
            year."),
          plotOutput(outputId = "checkExpPop"),
          VBE.ui("N", "WtnPop")
        ),
        # # Survival ----
        # tabPanel(
        #   title = "Survival",
        #   value = "phi.tab",
        #   TDVBE.ui(
        #     "phi", "Survival rates", "Individual annual survival rates.", "All"
        #   )
        # ),
        # Unknown parents ----
        tabPanel(
          title = "Unknown parents",
          value = "UPs.tab",
          h3("Unknown parents"),
          p("The first generation is simulated with unknown parents. 
          These individuals affect the observed numbers of kin-pairs, as 
          parent-offspring pairs among them, and sibling-pairs including them,
          are both unknown. This causes the appearance of prediction error, so
          it should be taken into account when evaluating predictor 
          performance."),
          p("Below are the average percentages of individuals with unknown
          parents, that are alive in survey-years, and pairs of survey-years. 
          The observed numbers of pairs of individuals, including one or more 
          with unknown parents, will be added soon."),
          tableOutput(outputId = "UPsWtn"),
          tableOutput(outputId = "UPsBtn"),
        ),
        # All pairs ----
        tabPanel(
          title = "All pairs",
          value = "APs.tab",
          TDVBE.ui(
            "APs", "All pairs", 
            "Total numbers of pairs of individuals.  The first numbers represent
            pairs in which both individuals are alive in the same survey-year,
            and the second represent pairs where one is alive in each of two
            different survey-years.", 
            c("WtnPop", "BtnPop")
          )
        ),
        # Self-pairs ----
        tabPanel(
          title = "Self-pairs",
          value = "SPs.tab",
          TDVBE.ui(
            "SPs", "Self-pairs", 
            "Numbers of pairs of individuals that are the same individual alive
            in two different survey-years.  The first numbers represent all such
            pairs, and the second represent only those with known parents, which
            are used in counting sibling-pairs (see unknown parents tab).", 
            c("BtnPop", "PrntsKwn")
          )
        ),
        # Parent-offspring pairs ----
        tabPanel(
          title = "Parent-offspring pairs",
          value = "POPs.tab",
          TDVBE.ui(
            "POPs", "Parent-offspring pairs", 
            "Numbers of pairs of individuals that are parent and offspring. 
            The first numbers represent pairs in which both individuals are
            alive in the same survey-year, and the second represent pairs where
            one is alive in each of two different survey-years.",
            c("WtnPop", "BtnPop")
          )
        ),
        # Same-mother pairs ----
        tabPanel(
          title = "Same-mother pairs",
          value = "SMPs.tab",
          # p("Numbers in the final year, with one born in
          #   the year indicated, and one born in the final year."),
          # p("Numbers between survey-years, one born five years before first, 
          # one born in last."),
          TDVBE.ui(
            "SMPs", "Same-mother pairs", 
            "Numbers of pairs of individuals with the same mother.  The first
            numbers represent pairs in which both individuals are
            alive in the same survey-year, and the second represent pairs where
            one is alive in each of two different survey-years.",
            # c("AgeKnwn", "WtnPop", "BtnPop", "BtnAgeKnwnPop")
            c("WtnPop", "BtnPop")
          )
        ),
        # Same-father pairs ----
        tabPanel(
          title = "Same-father pairs",
          value = "SFPs.tab",
          # p("Numbers in the final year, with one born in
          #   the year indicated, and one born in the final year."),
          # p("Numbers in the final year, both born in the year indicated."),
          TDVBE.ui(
            "SFPs", "Same-father pairs", 
            "Numbers of pairs of individuals with the same father.  The first
            numbers represent pairs in which both individuals are
            alive in the same survey-year, and the second represent pairs where
            one is alive in each of two different survey-years.",
            # c("AgeKnwn", "SameAge", "WtnPop", "BtnPop")
            c("WtnPop", "BtnPop")
          )
        ),
        # Full-sibling pairs ----
        tabPanel(
          title = "Full-sibling pairs",
          value = "FSPs.tab",
          TDVBE.ui(
            "FSPs", "Full-sibling pairs",
            "Numbers of pairs of individuals with the same parents.  The first
            numbers represent pairs in which both individuals are
            alive in the same survey-year, and the second represent pairs where
            one is alive in each of two different survey-years.",
            c("WtnPop", "BtnPop")
          )
        ),
        # Half-sibling pairs ----
        tabPanel(
          title = "Half-sibling pairs",
          value = "HSPs.tab",
          TDVBE.ui(
            "HSPs", "Half-sibling pairs",
            "Numbers of pairs of individuals that share exactly one parent.  The
            first numbers represent pairs in which both individuals are
            alive in the same survey-year, and the second represent pairs where
            one is alive in each of two different survey-years.",
            c("WtnPop", "BtnPop")
          )
        ),
        # Biases ----
        tabPanel(
          title = "Overall biases",
          value = "bias.tab",
          h2("Overall biases"),
          p("Average proportional differences between numbers simulated and 
            predicted, over survey-years and pairs of survey-years.  Some of
            these differences are affected by individuals with unknown parents
            (see unknown parents tab), so these are also reported below.  The 
            prediction for the total number of pairs of individuals (all pairs)
            seems to show a consistent difference.  The number of full-sibling
            pairs is much smaller than the others so the proportional difference
            from the prediction seems to be more variable."),
          h3("Unknown parents"),
          p("Average percentages of individuals with unknown parents, that are
            alive in survey-years, and pairs of survey-years."),
          tableOutput(outputId = "percUnknPrnts"),
          # h4("Temporal estimates"),
          # textOutput("tempEstBiasNote"),
          # tableOutput(outputId = "bsNsKPsTemp"),
          h3("Numbers of pairs within survey-years (whole population)"),
          p("Numbers of pairs of individuals with given relationships,
            where both are alive in the same survey-year.  Population
            sizes are included for reference."),
          tableOutput(outputId = "bsNsKPsWtnPop"),
          h3("Numbers of pairs between survey-years (whole population)"),
          p("Numbers of pairs of individuals with given relationships,
            where one individual is alive in each of a pair of survey-years."),
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
      h3("Biological scenario"),
      tableOutput(outputId = "modBioScen"),
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
