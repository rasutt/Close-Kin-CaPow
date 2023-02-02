# Load shiny package
library(shiny)

# Load high-precision numbers package
# library(Rmpfr)

# Function to make outputs for values, biases, and errors
VBE.ui <- function(id, types, descs) {
  ns <- NS(id)
  tagList(lapply(1:length(types), function(i) {
    type = types[i]
    list(
      h3(type),
      p(descs[i]),
      h4("Values"),
      plotOutput(ns(paste0("vals", type))),
      h4("Biases"),
      tableOutput(ns(paste0("bss", type))),
      h4("Errors"),
      plotOutput(ns(paste0("errs", type)))
    )
  }))
}

# Function to make outputs for titles, descriptions, values, biases, and errors
TDVBE.ui <- function(
  id, ttl_dsc, types = wtn_btn_headings, typ_dscs = wtn_btn_descs
) {
  tagList(h2(kp.nms[id]), p(ttl_dsc), VBE.ui(id, types, typ_dscs))
}

# Function to make kin-pair tab-panel
KP.tab.ui = function(id) {
  tabPanel(
    title = kp.nms[id],
    value = paste0(id, ".tab"),
    TDVBE.ui(id, rglr.kp.dscs[id])
  )
}

# Define UI for app
ui <- fluidPage(
  # Padding so that navigation bar doesn't cover title
  tags$style(type="text/css", "body {padding-top: 60px;}"),
  
  # Navigation bar with tabs
  navbarPage(
    title = "Close-kin CaPow!",
    id = "nav.tab",
    selected = "model.tab",
    position = "fixed-top",
    
    # Sim tab ----
    tabPanel(
      title = "Simulate studies",
      value = "sim.tab",
      
      # App title and description
      includeMarkdown("README.md"),
      
      # Side-bar layout separating inputs and outputs
      sidebarLayout(
        # Inputs ----
        sidebarPanel(
          sliderInput(
            inputId = "phi", label = HTML("Survival rate (&phi;):"),
            min_phi, max_phi, value = 0.95, step = step_phi
          ),
          sliderInput(
            inputId = "rho", label = HTML("Per capita birth rate (&rho;):"),
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
            value = "2010, 2015, 2020"
          ),
          sliderInput(
            inputId = "p", label = "Base level capture probability:",
            0, 0.2, value = 0.1, step = 0.01
          ),
          sliderInput(
            inputId = "L", label = "Number of SNP loci:",
            0, 1000, value = 10, step = 10
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
            1, 1000, value = 2
          ),
          actionButton(
            inputId = "simulate", label = "Simulate studies"
          )
        ),
        # Outputs ---- 
        mainPanel(
          h2("Next simulation features"),
          p("Expected population size and parameter values of next simulation
          implied by currently selected inputs."),
          h3("Expected population size"),
          plotOutput("nextExpPop"),
          h3("Implied parameters"),
          tableOutput("nextParsImpld")
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
        selected = "frst.lklhds.tb",
        # Simulation features ----
        tabPanel(
          title = "Simulation features",
          value = "sim.feat.tab",
          h2("Current simulation features"),
          p("Features of current simulation, selected and implied."),
          h3("Selected parameters"),
          tableOutput("currParsSltd"),
          h3("Implied parameters"),
          tableOutput("currParsImpld"),
          h3("Simulation options"),
          tableOutput("currSimOpts"),
          h3("Biological scenario"),
          tableOutput("currBioScen")
        ),
        # First study samples ----
        tabPanel(
          title = "First study samples",
          value = "frst.cps.tb",
          h2("First study samples"),
          p("Sample-histories from first population and study simulated."),
          
          h3("First sample-histories"),
          p("The data for the first few animals sampled
            in the study, by earliest birth year."),
          tableOutput("firstSampHists"),
          
          h3("Sufficient statistics"),
          p("Sufficient statistics for sample histories, used to fit popan 
            models."),
          tableOutput("firstSmpHstStats"),
          textOutput("firstNSmpHsts")
        ),
        # First study kin-pairs ----
        tabPanel(
          title = "First study kin-pairs",
          value = "frst.KPs.tb",
          h2("First study kin-pairs"),
          p("Numbers of kin-pairs predicted versus simulated in the first
            population and sampling study simulated, as an example.  Note that
            individual simulations can be quite different from each other, but 
            the average over many simulations should approach the expected 
            values."),
          
          h3("Kinpair probabilities predicted for first study from TMB"),
          p("Predictions given true parameter values"),
          h4("From TMB function"),
          tableOutput("firstKPPrbsTMB"),
          h4("From kin pair numbers checks"),
          tableOutput("firstKPPrbsR"),
          
          h3("Numbers of pairs within survey-years (whole population)"),
          p("Numbers of pairs of individuals with given relationships,
            where both are alive in the given survey-year.  Population
            sizes are included for reference."),
          h4("Predicted"),
          tableOutput("firstEstNsKPsWtnPop"),
          h4("Simulated"),
          tableOutput("firstNsKPsWtnPop"),
          
          h3("Numbers of pairs between survey-years (whole population)"),
          p("Numbers of pairs of individuals with given relationships,
            where one individual is alive in each of the given pair of 
            survey-years."),
          h4("Predicted"),
          tableOutput("firstEstNsKPsBtnPop"),
          h4("Simulated"),
          tableOutput("firstNsKPsBtnPop"),
        
          h3("Numbers of pairs within survey-years (sampled individuals)"),
          p("Numbers of pairs of individuals with given relationships,
            where both are sampled in the given survey-year.  Total numbers 
            sampled are included for reference."),
          h4("Predicted"),
          tableOutput("firstEstNsKPsWtnSmp"),
          h4("Simulated"),
          tableOutput("firstNsKPsWtnSmp"),
          
          h3("Numbers of pairs between survey-years (sampled individuals)"),
          p("Numbers of pairs of individuals with given relationships,
            where one individual is sampled in each of the given pair of 
            survey-years."),
          h4("Predicted"),
          tableOutput("firstEstNsKPsBtnSmp"),
          h4("Simulated"),
          tableOutput("firstNsKPsBtnSmp"),
        ),
        # First study genetics ----
        tabPanel(
          title = "First study genetics",
          value = "frst.gts.tb",
          h2("First study genetics"),
          p("Genetic analysis of samples from first population and study 
            simulated."),
          
          h3("First sample genotypes"),
          p("First few individuals sampled."),
          tableOutput("firstFSGs"),
          
          h3("Allele frequencies"),
          p("Relative frequencies, excluding samples from same animals in
             different surveys."),
          tableOutput("firstAFs"),
          
          h3("Genotype probabilities"),
          p("Probabilities of possible genotypes assuming that genes are 
            inherited with probabilities given by their frequencies above."),
          tableOutput("firstGPs"),
          
          h3("Genopair probabilities at first locus"),
          p("Probabilities of all possible genopairs given each kinship."),
          fluidRow(
            column(3, h4("Unrelated"), tableOutput("firstGPsUPL1")),
            column(3, h4("Half-siblings"), tableOutput("firstGPsHSPL1")),
            column(3, h4("Parent-offspring"), tableOutput("firstGPsPOPL1")),
            column(3, h4("Self-resample"), tableOutput("firstGPsSPL1"))
          ),
          
          h3("Genopair log-probabilities at first locus"),
          p("Log-probabilities of all possible genopairs given each kinship."),
          fluidRow(
            column(3, h4("Unrelated"), tableOutput("firstGLPsUPL1")),
            column(3, h4("Half-siblings"), tableOutput("firstGLPsHSPL1")),
            column(3, h4("Parent-offspring"), tableOutput("firstGLPsPOPL1")),
            column(3, h4("Self-resample"), tableOutput("firstGLPsSPL1"))
          ),
          
          h3("First sample genopair log-probabilities"),
          p("Observed values for first few sample-pairs, over all loci."),
          tableOutput("firstFSGLPs"),

          h3("All sample genopair log-probabilities"),
          p("Observed values for all sample-pairs, over all loci."),
          plotOutput("firstASGLPs"),
          
          h3("First sample genopair probabilities"),
          p("Possibly adjusted by a large positive factor to avoid underflow."),
          h4("All pairs"),
          p("Observed values for first few sample-pairs, over all loci."),
          tableOutput("firstFSGPsF"),

          p("Numbers of pairs with corresponding survey indices for first and
            second samples."),
          tableOutput("firstSYIPsF"),
          
          h4("Offset pairs"),
          p("Observed values for first few sample-pairs, over all loci (order is
          random to avoid bias due to age representation when individuals 
          repeated to include pairs between surveys with different numbers of 
            samples)."),
          tableOutput("firstFSGPsO"),
          p("Numbers of pairs with corresponding survey indices for first and
            second samples."),
          tableOutput("firstSYIPsO"),
          
          h3("Half-sibling vs unrelated pair PLODs"),
          p("Observed values for all sample-pairs, over all loci."),
          plotOutput("firstPLODs")
        ),
        # First study likelihoods ----
        tabPanel(
          title = "First study likelihoods",
          value = "frst.lklhds.tb",
          h2("First study likelihoods"),
          h3("Likelihood surface near true parameter values"),
          p("Negative log-likelihood over each parameter while others held at 
            true values."),
          plotOutput("firstNLLSurfs", height = 500)
        ),
        # First study estimates ----
        tabPanel(
          title = "First study estimates",
          value = "frst.ests.tb",
          h2("First study estimates"),
          p("Parameter estimates for first simulated study. Close-kin models 
             include self, parent-offspring, and half-sibling pairs."),
          h3("Results for first study"),
          tableOutput("firstResults"),

          h3("Kinpair probabilities estimated for first study"),
          h4("True kinship model"),
          tableOutput("firstKPPrbsTK")
        ),
        # Populations ----
        tabPanel(
          title = "Population sizes",
          value = "N.tab",
          h2("Population sizes"),
          p("Numbers of individuals that are alive in the population."),
          h3("Whole simulation"),
          plotOutput("checkExpPop"),
          VBE.ui("N", "In survey-years", "")
        ),
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
          tableOutput("UPsWtn"),
          tableOutput("UPsBtn"),
        ),
        # All pairs ----
        KP.tab.ui("APs"),
        # Self-pairs ----
        tabPanel(
          title = "Self-pairs",
          value = "SPs.tab",
          TDVBE.ui(
            "SPs", 
            "Numbers of pairs of individuals that are the same individual
              alive in two different survey-years.",
            c("All self-pairs", "Self-pairs with known parents"),
            c(
              "These include individuals with unknown parents (see unknown
              parents tab)",
              "These are excluded when counting sibling-pairs from different
              survey-years (see unknown parents tab)."
            )
          )
        ),
        # Parent-offspring pairs ----
        KP.tab.ui("POPs"),
        # Same-mother pairs ----
        KP.tab.ui("SMPs"),
        # Same-father pairs ----
        KP.tab.ui("SFPs"),
        # Full-sibling pairs ----
        KP.tab.ui("FSPs"),
        # Half-sibling pairs ----
        KP.tab.ui("HSPs"),
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
          tableOutput("percUnknPrnts"),
          # h4("Temporal estimates"),
          # textOutput("tempEstBiasNote"),
          # tableOutput("bsNsKPsTemp"),
          h3("Numbers of pairs within survey-years (whole population)"),
          p("Numbers of pairs of individuals with given relationships,
            where both are alive in the same survey-year.  Population
            sizes are included for reference."),
          tableOutput("bsNsKPsWtnPop"),
          h3("Numbers of pairs between survey-years (whole population)"),
          p("Numbers of pairs of individuals with given relationships,
            where one individual is alive in each of a pair of survey-years."),
          tableOutput("bsNsKPsBtnPop"),
          
          # h3("Probabilities"),
          # h4("Within surveys"),
          # tableOutput("bsProbsKPsWtn"),
          # h4("Between surveys"),
          # tableOutput("bsProbsKPsBtn"),
          
          # h3("Numbers among sampled animals"),
          # h4("Within surveys"),
          # tableOutput("bsNsKPsCapWtn"),
          # h4("Between surveys")
          # tableOutput("bsNsKPsCapBtn")
        ),
      )
    ), 
    # Model tab ----
    tabPanel(
      title = "Fit models",
      value = "model.tab",
      h2("Fit models"),
      sidebarLayout(
        sidebarPanel(
          checkboxGroupInput(
            inputId = "mdl.st", label = "Models to fit",
            choices = mdl.chcs, 
            # selected = mdl.chcs
            selected = c("Offset true kinship", "True kinship")
            # selected = c("True kinship", "Full genopair")
          ),
          checkboxGroupInput(
            inputId = "knshp.st", label = "Kinships to include",
            choices = knshp.chcs,
            selected = knshp.chcs
            # selected = c("Self", "Parent-offspring")
          ),
          actionButton(
            inputId = "fit", label = "Fit models"
          )
        ),
        mainPanel(
          h2("Compare model performance"),
          tableOutput("nDatasets"),
          tableOutput("knshpSt"),
          h3("Model fitting"),
          tableOutput("mdlFtRts"),
          h3("Estimator performance"),
          h4("Bias"),
          tableOutput("estBias"),
          h4("Coefficient of variation"),
          tableOutput("estCV"),
          h4("Estimates/errors"),
          plotOutput("mdlCmpnPlt"),
          h3("95% confidence interval coverage"),
          tableOutput("CICov"),
          h3("95% confidence intervals for lambda"),
          plotOutput("CIPlot")
        )
      )
    ),
    # Save/load tab ----
    tabPanel(
      title = "Save/Load",
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
