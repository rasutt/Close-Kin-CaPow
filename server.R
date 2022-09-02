# Load functions
funcs <- list.files("Functions")
for (i in 1:length(funcs)) source(paste0("Functions/", funcs[i]))

# Load TMB library and likelihood functions.  Have to restart R for compile and
# dyn.load to take effect sometimes.  Also sometimes need to update R so that
# Matrix package matches, and set the path for Rtools again.
library(TMB)
compile("TMB_objective_functions/POPANNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/POPANNLL"))
compile("TMB_objective_functions/CloseKinNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/CloseKinNLL"))

# Define server logic for app
server <- function(input, output) {
  # Reactive variables (for next simulation) ----
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
  # Parameter names
  par.names.rct = reactive({
    c("lambda", "phi", "Exp_N_final", "Exp_Ns", paste0("p", srvy.yrs.rct()))
  })
  # Parameter values
  par.vals.rct = reactive({
    c(
      lambda.rct(), input$phi, exp.N.t.rct()[input$hist.len], exp.Ns.rct(), 
      rep(input$p, k.rct())
    )
  })
  # Simulation options 
  sim.opts.rct = reactive({
    df = data.frame(
      input$n_sims, input$hist.len, input$clvng.ints, input$clvng.p, 
      input$tmp.emgn, input$alpha
    )
    names(df) = c(
      "Number of studies", "Population history length", 
      "Female time-order breeding", "Additional calving-capture probability", 
      "Male absense probability", "Age of sexual maturity"
    )
    df
  })
  # Models to fit
  models = reactive(input$models) 
  # ----
  
  # Variables bound to simulate button (for last simulation) ----
  # Individual survival rate
  phi <- bindEvent(reactive(input$phi), input$simulate, ignoreNULL = F)
  # Birthrate
  rho <- bindEvent(reactive(input$rho), input$simulate, ignoreNULL = F)
  # Population growth rate
  lambda <- bindEvent(lambda.rct, input$simulate, ignoreNULL = F)
  # Base year for expected population size
  base.yr = bindEvent(reactive(input$base.yr), input$simulate, ignoreNULL = F)
  # Expected population size in base year
  exp.N.base = bindEvent({
    reactive(input$exp.N.base)
  }, input$simulate, ignoreNULL = F)
  # Survey years
  srvy.yrs = bindEvent(srvy.yrs.rct, input$simulate, ignoreNULL = F)
  # Capture probability
  p = bindEvent(reactive(input$p), input$simulate, ignoreNULL = F)
  # Calving intervals
  clvng.ints = bindEvent({
    reactive(input$clvng.ints)
  }, input$simulate, ignoreNULL = F)
  # Additional capture probability when calving
  clvng.p = bindEvent(reactive(input$clvng.p), input$simulate, ignoreNULL = F)
  # Probability of males being away from survey area
  tmp.emgn = bindEvent(reactive(input$tmp.emgn), input$simulate, ignoreNULL = F)
  # Age of sexual maturity
  alpha = bindEvent(reactive(input$alpha), input$simulate, ignoreNULL = F)
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
  # Indices of survey years within population histories
  s.yr.inds <- bindEvent({
    reactive(hist.len() + srvy.yrs() - f.year())
  }, input$simulate, ignoreNULL = F) 
  # Expected population size over time
  exp.N.t = bindEvent(exp.N.t.rct, input$simulate, ignoreNULL = F)
  # Expected final population size
  exp.N.fin = bindEvent({
    reactive(exp.N.t()[hist.len()])
  }, input$simulate, ignoreNULL = F)
  # Expected super-population size
  exp.Ns = bindEvent(exp.Ns.rct, input$simulate, ignoreNULL = F) 
  # Simulation parameter names
  par.names = reactiveVal()
  par.names = bindEvent(par.names.rct, input$simulate, ignoreNULL = F) 
  # Parameter values
  par.vals = bindEvent(par.vals.rct, input$simulate, ignoreNULL = F) 
  # Estimate parameter names
  est.par.names = bindEvent({
    reactive(c("lambda", "phi", "N_final", "Ns", paste0("p", srvy.yrs())))
  }, input$simulate, ignoreNULL = F) 
  # Simulation options 
  sim.opts = bindEvent(sim.opts.rct, input$simulate, ignoreNULL = F) 
  # ----

  # Function to make data frame of parameter values for display
  par_vals_df = function(par_vals, par_names) {
    par_mat = matrix(par_vals, nrow = 1)
    colnames(par_mat) = par_names
    par_df = data.frame(par_mat)
    par_df[, 3] = as.integer(par_df[, 3])
    par_df[, 4] = as.integer(par_df[, 4])
    par_df
  }
  
  # Function to prepare proportion to print as percentage
  perc = function(prpn) paste0(round(prpn * 100, 1), "%")

  # Load data, functions, and outputs for simulating studies, checking
  # simulations, and analyzing model performance
  load("ckc_saved_objs.Rdata")
  par.names = reactiveVal(saved_objs$par.names)
  par.vals = reactiveVal(saved_objs$par.vals)
  sim.opts = reactiveVal(saved_objs$sim.opts)
  sim.lst = reactiveVal(saved_objs$sim.lst)
  checks.lst = reactiveVal(saved_objs$checks.lst)
  
  source("Tabs/sim_tab.R", local = T)
  source("Tabs/check_tab.R", local = T)
  source("Tabs/Check_sub_tabs/first_study.R", local = T)
  source("Tabs/Check_sub_tabs/kin_pairs.R", local = T)
  source("Tabs/Check_sub_tabs/temp_ests.R", local = T)
  source("Tabs/Check_sub_tabs/err_dists.R", local = T)
  source("Tabs/model_tab.R", local = T)
  source("Tabs/save_load_tab.R", local = T)
}