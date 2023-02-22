# Load functions
funcs <- list.files("Functions")
for (i in 1:length(funcs)) source(paste0("Functions/", funcs[i]))

# Load TMB library and likelihood functions.  Have to restart R for compile and
# dyn.load to take effect sometimes.  Also sometimes need to update R so that
# Matrix package matches, or actually delete and recompile files after
# reinstalling TMB
library(TMB)
compile("TMB_objective_functions/UnifiedNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/UnifiedNLL"))

# Load parallel library for multicore processing
library(parallel)

# Define server logic for app
server <- function(input, output) {
  # Reactive variables (for next simulation).  Can ignore warnings for invalid
  # inputs while typing survey years ---- 
  # Population growth rate
  lambda.rct <- reactive(input$rho + input$phi) 
  # Implied birth rate for mature females surviving to birth year
  beta.rct = reactive({
    2 * (1 - input$phi / lambda.rct()) * (lambda.rct() / input$phi)^input$alpha
  })
  # Survey years
  srvy.yrs.rct = reactive({
    sort(eval(parse(text = paste0("c(", input$srvy.yrs, ")"))))
  })
  # Number of surveys 
  k.rct = reactive(length(srvy.yrs.rct()))
  # Survey pairs
  srvy.prs.rct = reactive({
    apply(combn(srvy.yrs.rct(), 2), 2, paste, collapse = "-")
  })
  # Number of pairs of surveys
  n.srvy.prs.rct = reactive(length(srvy.prs.rct()))
  # Final survey/simulation year 
  fnl.year.rct = reactive(tail(srvy.yrs.rct(), 1))  
  # First simulation year
  fst.year.rct = reactive(fnl.year.rct() - input$hist.len + 1)
  # Simulation years
  sim.yrs.rct = reactive(fst.year.rct():fnl.year.rct())
  # Indices of survey years within population histories
  s.yr.inds.rct = reactive(input$hist.len + srvy.yrs.rct() - fnl.year.rct())
  # Survey gaps
  srvy.gaps.rct <- reactive(as.integer(diff(srvy.yrs.rct())))
  # Expected population size over time
  exp.N.t.rct = reactive({
    input$exp.N.base * lambda.rct()^(fnl.year.rct() - input$base.yr) / 
      lambda.rct()^((input$hist.len - 1):0)
  })
  # Initial population size
  N.init.rct = reactive(round(exp.N.t.rct()[1]))
  # Expected super-population size
  exp.Ns.rct = reactive({
    exp.N.srvy.yrs = 
      exp.N.t.rct()[input$hist.len - fnl.year.rct() + srvy.yrs.rct()]
    sum(exp.N.srvy.yrs) - sum(exp.N.srvy.yrs[-length(srvy.yrs.rct())] * 
                                input$phi^(srvy.gaps.rct()))
  })
  # Parameter names
  par.names.rct = reactive({
    c("lambda", "phi", "Exp N final", "Exp Ns", paste0("p", srvy.yrs.rct()))
  })
  # Parameter values
  par.vals.rct = reactive({
    c(
      lambda.rct(), input$phi, exp.N.t.rct()[input$hist.len], exp.Ns.rct(), 
      rep(input$p, k.rct())
    )
  })

  # Load saved objects ----
  load("Datasets/ckc_saved_objs.Rdata")
  phi = reactiveVal(saved.objs$phi)
  rho = reactiveVal(saved.objs$rho)
  lambda = reactiveVal(saved.objs$lambda)
  beta = reactiveVal(saved.objs$beta)
  exp.N.base = reactiveVal(saved.objs$exp.N.base)
  base.yr = reactiveVal(saved.objs$base.yr)
  srvy.yrs = reactiveVal(saved.objs$srvy.yrs)
  p = reactiveVal(saved.objs$p)
  L = reactiveVal(saved.objs$L)
  clvng.ints = reactiveVal(saved.objs$clvng.ints)
  clvng.p = reactiveVal(saved.objs$clvng.p)
  tmp.emgn = reactiveVal(saved.objs$tmp.emgn)
  alpha = reactiveVal(saved.objs$alpha)
  hist.len = reactiveVal(saved.objs$hist.len)
  n.sims.rqd = reactiveVal(saved.objs$n.sims.rqd)
  n.sims = reactiveVal(saved.objs$n.sims)
  srvy.gaps = reactiveVal(saved.objs$srvy.gaps)
  k = reactiveVal(saved.objs$k)
  n.srvy.prs = reactiveVal(saved.objs$n.srvy.prs)
  srvy.prs = reactiveVal(saved.objs$srvy.prs)
  fnl.year = reactiveVal(saved.objs$fnl.year)
  fst.year = reactiveVal(saved.objs$fst.year)
  sim.yrs = reactiveVal(saved.objs$sim.yrs)
  n.yrs.chk.t = reactiveVal(saved.objs$n.yrs.chk.t)
  yrs.chk.t = reactiveVal(saved.objs$yrs.chk.t)
  s.yr.inds = reactiveVal(saved.objs$s.yr.inds)
  exp.N.t = reactiveVal(saved.objs$exp.N.t)
  N.init = reactiveVal(saved.objs$N.init)
  exp.N.fin = reactiveVal(saved.objs$exp.N.fin)
  exp.Ns = reactiveVal(saved.objs$exp.Ns)
  par.names = reactiveVal(saved.objs$par.names)
  par.vals = reactiveVal(saved.objs$par.vals)
  est.par.names = reactiveVal(saved.objs$est.par.names)
  sim.opts.bio.scen = reactiveVal(saved.objs$sim.opts.bio.scen)
  sim.lst = reactiveVal(saved.objs$sim.lst)
  checks.lst = reactiveVal(saved.objs$checks.lst)
  fit.lst = reactiveVal(saved.objs$fit.lst)
  N.t.mat = reactiveVal(saved.objs$N.t.mat)
  avg.phi.obs = reactiveVal(saved.objs$avg.phi.obs)
  ns.SPs = reactiveVal(saved.objs$ns.SPs)
  ns.POPs = reactiveVal(saved.objs$ns.POPs)
  ns.SMPs = reactiveVal(saved.objs$ns.SMPs)
  ns.SFPs = reactiveVal(saved.objs$ns.SFPs)
  ns.SMPs.t = reactiveVal(saved.objs$ns.SMPs.t)
  ns.SFPs.t = reactiveVal(saved.objs$ns.SFPs.t)
  ns.SibPs = reactiveVal(saved.objs$ns.SibPs)
  pns.UPs = reactiveVal(saved.objs$pns.UPs)
  frst.fglps = reactiveVal(saved.objs$frst.fglps)
  mdl.st = reactiveVal(saved.objs$mdl.st)
  knshp.st = reactiveVal(saved.objs$knshp.st)
  osisyips = reactiveVal(saved.objs$osisyips)
  
  # Create empty multicore cluster object, as needs to be updated when new
  # datasets loaded/simulated, and new models fit
  cl = reactiveVal(NULL)
  
  # Variables bound to simulate button (for last simulation) ----
  observeEvent(input$simulate, {
    # Individual survival rate
    phi(input$phi)
    # Birthrate
    rho(input$rho)
    # Population growth rate
    lambda(lambda.rct())
    # Implied birth rate for mature females surviving to birth year
    beta(beta.rct())
    # Expected population size in base year
    exp.N.base(input$exp.N.base)
    # Base year for expected population size
    base.yr(input$base.yr)
    # Survey years
    srvy.yrs(srvy.yrs.rct())
    # Capture probability
    p(input$p)
    # Number of SNP loci
    L(input$L)
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
    # Survey gaps
    srvy.gaps(srvy.gaps.rct())
    # Number of surveys
    k(k.rct())
    # Survey pairs
    srvy.prs(srvy.prs.rct())
    # Number of pairs of surveys
    n.srvy.prs(n.srvy.prs.rct())
    # Final survey/simulation year
    fnl.year(fnl.year.rct())
    # First simulation year
    fst.year(fst.year.rct())
    # Simulation years
    sim.yrs(sim.yrs.rct())
    # Number of years to check temporal estimates
    n.yrs.chk.t(min(hist.len() - 1, n.yrs.try.chk.t))
    # Years to check temporal estimates
    yrs.chk.t(sim.yrs()[(hist.len() - n.yrs.chk.t()):(hist.len() - 1)])
    # Indices of survey years within population histories
    s.yr.inds(s.yr.inds.rct())
    # Expected population size over time
    exp.N.t(exp.N.t.rct())
    # Initial population size
    N.init(N.init.rct())
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
    # Simulation options and biological scenario
    sim.opts.bio.scen({
      df = data.frame(
        input$n.sims.rqd, input$hist.len, input$clvng.ints, input$clvng.p, 
        input$tmp.emgn, input$alpha
      )
      names(df) = c(
        "Number of studies requested", "Population history length", 
        "Female time-order breeding", "Additional calving-capture probability", 
        "Male absense probability", "Age of sexual maturity"
      )
      df
    })
  })
  
  # Load functions, and outputs for simulating studies, checking
  # simulations, and analyzing model performance ----
  source("Tabs/1_sim_tab.R", local = T)
  source("Tabs/2_1_check_tab.R", local = T)
  source("Tabs/2_2_preds_and_errs.R", local = T)
  source("Tabs/Check_sub_tabs/0_sim_feats.R", local = T)
  source("Tabs/Check_sub_tabs/1_0_first_caps.R", local = T)
  source("Tabs/Check_sub_tabs/1_1_first_KPs.R", local = T)
  source("Tabs/Check_sub_tabs/1_2_first_gts.R", local = T)
  source("Tabs/Check_sub_tabs/1_3_first_lklhds.R", local = T)
  source("Tabs/Check_sub_tabs/1_4_first_ests.R", local = T)
  source("Tabs/Check_sub_tabs/2_pops_gts_and_UPs.R", local = T)
  source("Tabs/Check_sub_tabs/3_kin_pair_tabs.R", local = T)
  source("Tabs/Check_sub_tabs/4_bias.R", local = T)
  source("Tabs/3_1_fglps.R", local = T)
  source("Tabs/3_2_model_fits.R", local = T)
  source("Tabs/3_3_model_outputs.R", local = T)
  source("Tabs/4_save_load_tab.R", local = T)
}
