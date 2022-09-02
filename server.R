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
  
  # Load saved objects ----
  load("ckc_saved_objs.Rdata")
  phi = reactiveVal(saved_objs$phi)
  rho = reactiveVal(saved_objs$rho)
  lambda = reactiveVal(saved_objs$lambda)
  base.yr = reactiveVal(saved_objs$base.yr)
  exp.N.base = reactiveVal(saved_objs$exp.N.base)
  srvy.yrs = reactiveVal(saved_objs$srvy.yrs)
  p = reactiveVal(saved_objs$p)
  clvng.ints = reactiveVal(saved_objs$clvng.ints)
  clvng.p = reactiveVal(saved_objs$clvng.p)
  tmp.emgn = reactiveVal(saved_objs$tmp.emgn)
  alpha = reactiveVal(saved_objs$alpha)
  hist.len = reactiveVal(saved_objs$hist.len)
  n_sims = reactiveVal(saved_objs$n_sims)
  srvy.gaps = reactiveVal(saved_objs$srvy.gaps)
  k = reactiveVal(saved_objs$k)
  n.srvy.prs = reactiveVal(saved_objs$n.srvy.prs)
  f.year = reactiveVal(saved_objs$f.year)
  s.yr.inds = reactiveVal(saved_objs$s.yr.inds)
  exp.N.t = reactiveVal(saved_objs$exp.N.t)
  exp.N.fin = reactiveVal(saved_objs$exp.N.fin)
  exp.Ns = reactiveVal(saved_objs$exp.Ns)
  par.names = reactiveVal(saved_objs$par.names)
  par.vals = reactiveVal(saved_objs$par.vals)
  est.par.names = reactiveVal(saved_objs$est.par.names)
  sim.opts = reactiveVal(saved_objs$sim.opts)
  sim.lst = reactiveVal(saved_objs$sim.lst)
  checks.lst = reactiveVal(saved_objs$checks.lst)
  # ----

  # Variables bound to simulate button (for last simulation) ----
  observeEvent(input$simulate, {
    # Individual survival rate
    phi(input$phi)
    # Birthrate
    rho(input$rho)
    # Population growth rate
    lambda(lambda.rct())
    # Base year for expected population size
    base.yr(input$base.yr)
    # Expected population size in base year
    exp.N.base(input$exp.N.base)
    # Survey years
    srvy.yrs(srvy.yrs.rct())
    # Capture probability
    p(input$p)
    # Calving intervals
    clvng.ints(input$clvng.ints)
    # Additional capture probability when calving
    clvng.p(input$clvng.p)
    # Probability of males being away from survey area
    tmp.emgn(input$tmp.emgn)
    # Age of sexual maturity
    alpha(input$alpha)
    # Length of simulation
    hist.len(input$hist.len)
    # Number of simulations
    n_sims(input$n_sims)
    # Survey gaps
    srvy.gaps(srvy.gaps.rct())
    # Number of surveys
    k(k.rct())
    # Number of pairs of surveys
    n.srvy.prs(choose(k(), 2))
    # Final survey year
    f.year(f.year.rct())
    # Indices of survey years within population histories
    s.yr.inds(hist.len() + srvy.yrs() - f.year())
    # Expected population size over time
    exp.N.t(exp.N.t.rct())
    # Expected final population size
    exp.N.fin(exp.N.t()[hist.len()])
    # Expected super-population size
    exp.Ns(exp.Ns.rct())
    # Parameter names
    par.names(par.names.rct())
    # Parameter values
    par.vals(par.vals.rct())
    # Estimate parameter names
    est.par.names(c("lambda", "phi", "N_final", "Ns", paste0("p", srvy.yrs())))
    # Simulation options 
    sim.opts(sim.opts.rct())
  })
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
  
  # Load functions, and outputs for simulating studies, checking
  # simulations, and analyzing model performance
  source("Tabs/sim_tab.R", local = T)
  source("Tabs/check_tab.R", local = T)
  source("Tabs/Check_sub_tabs/first_study.R", local = T)
  source("Tabs/Check_sub_tabs/kin_pairs.R", local = T)
  source("Tabs/Check_sub_tabs/temp_ests.R", local = T)
  source("Tabs/Check_sub_tabs/err_dists.R", local = T)
  source("Tabs/model_tab.R", local = T)
  source("Tabs/save_load_tab.R", local = T)
}