# Load functions
funcs <- list.files("Functions")
for (i in 1:length(funcs)) source(paste0("Functions/", funcs[i]))

# Load TMB library and likelihood functions.  Have to restart R for compile and
# dyn.load to take effect sometimes!
library(TMB)
compile("TMB_objective_functions/POPANNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/POPANNLL"))
compile("TMB_objective_functions/CloseKinNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/CloseKinNLL"))

# Define server logic for app
server <- function(input, output) {
  # Reactive variables ----
  # Population growth rate
  lambda.rct <- reactive(input$rho + input$phi) 
  # Survey years
  srvy.yrs.rct = reactive({
    eval(parse(text = paste0("c(", sort(input$srvy.yrs), ")")))
  })
  # Number of surveys 
  k.rct = reactive(length(srvy.yrs.rct()))
  # Final survey year 
  f.year.rct = reactive(tail(srvy.yrs.rct(), 1))  
  # Survey gaps
  srvy.gaps.rct <- reactive(as.integer(diff(srvy.yrs.rct())))
  # Expected population size over time
  exp.N.t.rct = reactive({
    input$exp.N.base * lambda.rct()^(f.year.rct() - input$base.yr) / 
      lambda.rct()^((input$hist.len - 1):0)
  })
  # Expected super-population size
  exp.Ns.rct = reactive({
    exp.N.srvy.yrs = 
      exp.N.t.rct()[input$hist.len - f.year.rct() + srvy.yrs.rct()]
    sum(exp.N.srvy.yrs) - sum(exp.N.srvy.yrs[-length(srvy.yrs.rct())] * 
      input$phi^(srvy.gaps.rct()))
  })
  # Models to fit
  models = reactive(input$models) 
  # ----
  
  # Variables bound to simulate button ----
  # Individual survival rate
  phi <- bindEvent(reactive(input$phi), input$simulate, ignoreNULL = F)
  # Birthrate
  rho <- bindEvent(reactive(input$rho), input$simulate, ignoreNULL = F)
  # Population growth rate
  lambda <- bindEvent(lambda.rct, input$simulate, ignoreNULL = F)
  # Survey years
  srvy.yrs = bindEvent(srvy.yrs.rct, input$simulate, ignoreNULL = F)
  # Length of simulation
  hist.len = bindEvent(reactive(input$hist.len), input$simulate, ignoreNULL = F)
  # Number of simulations
  n_sims = bindEvent(reactive(input$n_sims), input$simulate, ignoreNULL = F)
  # Survey gaps
  srvy.gaps <- bindEvent(srvy.gaps.rct, input$simulate, ignoreNULL = F)
  # Number of surveys
  k <- bindEvent(k.rct, input$simulate, ignoreNULL = F)
  # Number of pairs of surveys
  n.srvy.prs <- reactive(choose(k(), 2))
  # Final survey year
  f.year <- bindEvent(f.year.rct, input$simulate, ignoreNULL = F) 
  # Expected population size over time
  exp.N.t = bindEvent(exp.N.t.rct, input$simulate, ignoreNULL = F) 
  # Expected super-population size
  exp.Ns = bindEvent(exp.Ns.rct, input$simulate, ignoreNULL = F) 
  # ----

  # Load functions and outputs for simulating studies, checking simulations, and
  # analyzing model performance
  source("Tabs/sim_tab.R", local = T)
  source("Tabs/check_tab.R", local = T)
  source("Tabs/model_tab.R", local = T)
}