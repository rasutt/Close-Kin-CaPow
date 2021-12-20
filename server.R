# Load functions
funcs <- list.files("Functions")
for (i in 1:length(funcs)) source(paste0("Functions/", funcs[i]))

# Define server logic for app
server <- function(input, output) {
  # Reactive variables ----
  # Survey years
  srvy.yrs.rct = reactive({
    eval(parse(text = paste0("c(", sort(input$srvy.yrs), ")")))
  })
  # Number of surveys 
  k.rct = reactive(length(srvy.yrs.rct()))
  # Final survey year 
  f.year.rct = reactive(tail(srvy.yrs.rct(), 1))
  # Models to fit
  models = reactive(input$models) 
  # ----
  
  # Variables bound to simulate button ----
  # Birthrate
  rho <- bindEvent(reactive(input$rho), input$simulate, ignoreNULL = F)
  # Individual survival rate
  phi <- bindEvent(reactive(input$phi), input$simulate, ignoreNULL = F)
  # Population growth rate
  lambda <- reactive(rho() + phi()) 
  # Survey years
  srvy.yrs = bindEvent(srvy.yrs.rct, input$simulate, ignoreNULL = F)
  # Length of simulation
  hist.len = bindEvent(reactive(input$hist.len), input$simulate, ignoreNULL = F)
  # Number of simulations
  n_sims = bindEvent(reactive(input$n_sims), input$simulate, ignoreNULL = F)
  # Survey gaps
  srvy.gaps <- reactive(as.integer(diff(srvy.yrs())))
  # Number of surveys
  k <- bindEvent(k.rct, input$simulate, ignoreNULL = F)
  # Number of pairs of surveys
  n.srvy.prs <- reactive(choose(k(), 2))
  # Final survey year
  f.year <- bindEvent(f.year.rct, input$simulate, ignoreNULL = F) 
  # Expected population size over time
  exp.N.t = reactive({
    # Expected final population size
    phi.gaps <- phi()^srvy.gaps()
    pents = pent_func(lambda()^srvy.gaps(), phi.gaps, k())
    exp.N.fin <- sum(pents * exp.Ns * prod(phi.gaps) / cumprod(c(1, phi.gaps)))
    
    exp.N.fin / lambda()^((hist.len() - 1):0)
  })
  # ----

  # Load functions and outputs for simulating studies, checking simulations, and
  # analyzing model performance
  source("Tabs/sim_tab.R", local = T)
  source("Tabs/check_tab.R", local = T)
  source("Tabs/model_tab.R", local = T)
}