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
      plotOutput(outputId = ns(paste0("vals", type))),
      h4("Biases"),
      tableOutput(outputId = ns(paste0("bss", type))),
      h4("Errors"),
      plotOutput(outputId = ns(paste0("errs", type)))
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
            value = "2000, 2010, 2020"
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
          plotOutput(outputId = "nextExpPop"),
          h3("Implied parameters"),
          tableOutput(outputId = "nextParsImpld")
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
        selected = "frst.KPs.tb",
        # Simulation features ----
        tabPanel(
          title = "Simulation features",
          value = "sim.feat.tab",
          h2("Current simulation features"),
          p("Features of current simulation, selected and implied."),
          h3("Selected parameters"),
          tableOutput(outputId = "currParsSltd"),
          h3("Implied parameters"),
          tableOutput(outputId = "currParsImpld"),
          h3("Simulation options"),
          tableOutput(outputId = "currSimOpts"),
          h3("Biological scenario"),
          tableOutput(outputId = "currBioScen")
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
          tableOutput(outputId = "firstSampHists"),
          
          h3("Sufficient statistics"),
          p("Sufficient statistics for sample histories, used to fit popan 
            models."),
          tableOutput(outputId = "firstSmpHstStats"),
          textOutput(outputId = "firstNSmpHsts")
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
          tableOutput(outputId = "firstNsKPsBtnPop"),
        
          h3("Numbers of pairs within survey-years (sampled individuals)"),
          p("Numbers of pairs of individuals with given relationships,
            where both are sampled in the given survey-year.  Total numbers 
            sampled are included for reference."),
          h4("Predicted"),
          tableOutput(outputId = "firstEstNsKPsWtnSmp"),
          h4("Simulated"),
          tableOutput(outputId = "firstNsKPsWtnSmp"),
          
          h3("Numbers of pairs between survey-years (sampled individuals)"),
          p("Numbers of pairs of individuals with given relationships,
            where one individual is sampled in each of the given pair of 
            survey-years."),
          h4("Predicted"),
          tableOutput(outputId = "firstEstNsKPsBtnSmp"),
          h4("Simulated"),
          tableOutput(outputId = "firstNsKPsBtnSmp"),
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
          tableOutput(outputId = "firstGTs"),
          
          h3("Allele frequencies"),
          p("Relative frequencies, excluding samples from same animals in
             different surveys."),
          tableOutput(outputId = "firstAFs"),
          
          h3("Genopair probabilities given kinship"),
          p("Possible values at first locus."),
          fluidRow(
            column(3, h5("Unrelated"), tableOutput(outputId = "firstGPPsUP")),
            column(
              3, h5("Half-siblings"), tableOutput(outputId = "firstGPPsHSP")
            ),
            column(
              3, h5("Parent-offspring"), tableOutput(outputId = "firstGPPsPOP")
            ),
            column(
              3, h5("Self-resample"), tableOutput(outputId = "firstGPPsSP")
            )
          ),
          
          h3("Genopair log-probabilities given kinship"),
          p("Observed values for first few sample-pairs, over all loci."),
          tableOutput(outputId = "firstFewLGPPs"),
          p("Observed values for all sample-pairs, over all loci."),
          plotOutput(outputId = "firstLGPPs"),

          h3("Half-sibling vs unrelated pair PLODs"),
          p("Observed values for all sample-pairs, over all loci."),
          fluidRow(
            column(6, plotOutput(outputId = "firstPLODs")),
            column(6, plotOutput(outputId = "firstPLODsRare"))
          ),
          
          h3("Genopair probabilities given kinship"),
          h4("All pairs"),
          p("Observed values for first few sample-pairs, over all loci."),
          tableOutput(outputId = "firstFewGPPsFll"),
          p("Numbers of pairs with corresponding survey indices for first and
            second samples."),
          tableOutput(outputId = "frstSYIPCntsFll"),
          
          h4("Offset pairs"),
          p("Observed values for first few sample-pairs, over all loci (order is 
          random to avoid bias due to age representation when individuals 
          repeated to include pairs between surveys with different numbers of 
            samples)."),
          tableOutput(outputId = "firstFewGPPsOffst"),
          p("Numbers of pairs with corresponding survey indices for first and
            second samples."),
          tableOutput(outputId = "frstSYIPCntsOffst"),
          
          p("Observed values for all sample-pairs, over all loci."),
          plotOutput(outputId = "frstGpPs")
        ),
        # First study estimates ----
        tabPanel(
          title = "First study estimates",
          value = "frst.ests.tb",
          h2("First study estimates"),
          p("Parameter estimates for first simulated study."),

          h3("Likelihood surface near true parameter values"),
          p("Negative log-likelihood over each parameter while others held at 
            true values."),
          h4("Popan likelihood"),
          plotOutput(outputId = "firstPpnNLLSurfs"),
          h4("True kinship likelihood"),
          plotOutput(outputId = "firstCKNLLSurfs"),
          h4("Full genopair likelihood"),
          plotOutput(outputId = "firstFGPNLLSurfs"),
          h4("Offset genopair likelihood"),
          plotOutput(outputId = "firstOGPNLLSurfs"),
          
          h3("Results for first study"),
          tableOutput(outputId = "firstResults")
        ),
        # Populations ----
        tabPanel(
          title = "Population sizes",
          value = "N.tab",
          h2("Population sizes"),
          p("Numbers of individuals that are alive in the population."),
          h3("Whole simulation"),
          plotOutput(outputId = "checkExpPop"),
          VBE.ui("N", "In survey-years", "")
        ),
        # Genotypes ----
        tabPanel(
          title = "Genotypes",
          value = "gt.tab",
          h2("Genotypes"),
          p("Genotypes simulated.")
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
          tableOutput(outputId = "UPsWtn"),
          tableOutput(outputId = "UPsBtn"),
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
      title = "Fit models",
      value = "model.tab",
      h2("Fit models"),
      sidebarLayout(
        sidebarPanel(
          checkboxGroupInput(
            inputId = "mdl.st", label = "Models to fit",
            choices = mdl.chcs, 
            # selected = mdl.chcs
            selected = "Full genopair"
          ),
          checkboxGroupInput(
            inputId = "knshp.st", label = "Kinships to include",
            choices = knshp.chcs,
            selected = knshp.chcs
          ),
          actionButton(
            inputId = "fit", label = "Fit models"
          )
        ),
        mainPanel(
          h2("Compare model performance"),
          tableOutput(outputId = "nDatasets"),
          tableOutput(outputId = "knshpSt"),
          h3("Model fitting success rates"),
          tableOutput(outputId = "modStats"),
          h3("Estimates"),
          plotOutput(outputId = "modComp"),
          h3("95% confidence interval coverage"),
          tableOutput(outputId = "CICov"),
          h3("95% confidence intervals for lambda"),
          plotOutput(outputId = "CIPlot")
        )
      )
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
