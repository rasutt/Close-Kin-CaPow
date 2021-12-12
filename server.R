# Define server logic for app
server <- function(input, output) {
  # Simulate population and capture histories
  hists.lst = reactive({
    # Population growth rate
    lambda <- input$lambda 
    
    # "Population birthrate"
    rho <- lambda - phi 
    
    # Load simulation and POPAN functions
    source("Functions/SimPopStud.R", local = T)
    source("Functions/PopanFuncs.R", local = T)
    
    # Expected final population size.  gaps variables only used here
    lambda.gaps <- lambda^srvy.gaps # Population growth rate between surveys
    phi.gaps <- phi^srvy.gaps # Individual survival rate between surveys
    exp.N.fin <- sum(pent_func(lambda.gaps, phi.gaps) * exp.Ns * 
                       prod(phi.gaps) / cumprod(c(1, phi.gaps)))
    
    # Expected population size over time
    exp.N.t <- exp.N.fin / lambda^((hist.len - 1):0)
    N.init <- round(exp.N.t[1]) # Initial population size
    
    # Create list for population and capture histories
    hists.lst <- vector("list", n.stds.sim)
    
    # Display progress
    cat("Simulating study: ")
    
    # Loop over histories
    for (hist.ind in 1:n.stds.sim) {
      # Display progress
      if (hist.ind %% 100 == 1) cat(hist.ind, "")
      
      # Simulate family and capture histories of population of animals over time
      hists.lst[[hist.ind]] <- SimPopStud()
    }
    
    hists.lst
  })
  
  # Print head of first study
  output$dataHead <- renderTable(head(data.frame(hists.lst()[[1]])))
  
  # Plot first population size over time
  output$popPlot <- renderPlot({
    # Population growth rate
    lambda <- input$lambda 
    
    # "Population birthrate"
    rho <- lambda - phi 
    
    # Load POPAN functions
    source("Functions/PopanFuncs.R", local = T)
    
    # Expected final population size.  gaps variables only used here
    lambda.gaps <- lambda^srvy.gaps # Population growth rate between surveys
    phi.gaps <- phi^srvy.gaps # Individual survival rate between surveys
    exp.N.fin <- sum(pent_func(lambda.gaps, phi.gaps) * exp.Ns * 
                       prod(phi.gaps) / cumprod(c(1, phi.gaps)))
    
    # Expected population size over time
    exp.N.t <- exp.N.fin / lambda^((hist.len - 1):0)

    # Plot first population size over time
    plot((f.year - 79):f.year, attributes(hists.lst()[[1]])$N.t.vec, t = "l", 
         main = "First population size over time",
         ylab = "E(Nt)", xlab = "Year")
    lines((f.year - 79):f.year, exp.N.t, col = 2)
  })
  
  # Plot negative log-likelihood surface for first study
  output$NLLPlot <- renderPlot({
    # Source NLL and kin pairs functions
    source("Functions/CloseKinNLL.R", local = T)
    source("Functions/FindNsKinPairs.R", local = T)
    
    # Get simulated family and capture histories of population of animals over
    # time
    pop.cap.hist <- hists.lst()[[1]]
    
    # Get numbers of animals captured in each survey
    ns.caps <- attributes(pop.cap.hist)$ns.caps
    
    # Find numbers of kin pairs
    ns.kps.lst <- FindNsKinPairs()
    
    # Create grids of parameter and NLL values
    nll_grid = par_grid = seq(min_lambda, max_lambda, step_lambda)
    
    # MLEs
    params = mod.ests.lst()[[1]][1, 1:3]
    
    # Find NLL over grid of parameter values
    for (i in seq_along(nll_grid)) {
      # Rho = lambda - phi
      params[1] = par_grid[i] - params[2]
      nll_grid[i] = CloseKinNLL(params)
    }
    
    # Plot NLL
    plot(
      par_grid, nll_grid,
      main = "Negative log-likelihood for first study",
      xlab = "lambda", ylab = "NLL", type = 'l'
    )
    abline(v = input$lambda, col = 2)
    abline(v = mod.ests.lst()[[1]][1, 1], col = 4)
    # abline(v = cis()[1, 1], col = 4, lty = 2)
    # abline(v = cis()[2, 1], col = 4, lty = 2)
  })
  
  # Find parameter estimates
  mod.ests.lst = reactive({
    # Population growth rate
    lambda <- input$lambda 
    
    # "Population birthrate"
    rho <- lambda - phi 
    
    # Load functions
    funcs <- list.files("Functions")
    for (i in 1:length(funcs)) source(paste0("Functions/",funcs[i]), local = T)
    
    # Create general optimizer starting-values and bounds, NAs filled in below
    ck.start <- c(rho, phi, NA)
    ck.lwr <- c(0, 0.75, NA)
    ck.upr <- c(0.35, 1, Inf)
    
    # Create vectors for superpopulation and final population sizes
    Ns.vec <- N.fin.vec <- numeric(n.stds.fit)
    
    # Create matrices for estimates
    ck.tmb.ests <- matrix(nrow = n.stds.fit, ncol = 5 + k, dimnames = list(
      NULL, c("lambda", "phi", "N_final", "Ns", paste0("p", 1:k), "cnvg")))
    
    # Loop over histories
    for (hist.ind in 1:n.stds.fit) {
      # Display progress
      cat("History:", hist.ind, "\n")
      
      # Get simulated family and capture histories of population of animals over
      # time
      pop.cap.hist <- hists.lst()[[hist.ind]]
      
      # Store superpopulation and final population size
      Ns.vec[hist.ind] <- attributes(pop.cap.hist)$Ns
      N.fin.vec[hist.ind] <- attributes(pop.cap.hist)$N.t.vec[hist.len]
      
      # Get numbers of animals captured in study and each survey
      n.cap.hists <- nrow(pop.cap.hist)
      ns.caps <- attributes(pop.cap.hist)$ns.caps
      
      # Find numbers of kin pairs
      ns.kps.lst <- FindNsKinPairs()
      
      # Update optimiser starting-values and bounds
      ck.start[3] <- N.fin.vec[hist.ind]
      ck.lwr[3] <- ns.caps[k]
      
      # Try to fit models
      ck.tmb.ests[hist.ind, -(5:(4 + k))] <- TryCloseKinTMB()
    }
    
    # Combine model estimates as list
    mod.ests.lst <- list(ck_tmb = ck.tmb.ests)  
    mod.ests.lst
  })
  
  # Print first few estimates
  output$firstEsts <- renderTable(
    lapply(mod.ests.lst(), function(mod.ests) head(round(mod.ests, 3)))[[1]]
  )
  
  # Plot estimates using model comparison plot function (though only one model
  # atm)
  output$modComp <- renderPlot({
    # Find optimization attempts that converged
    cvgd.ests.lst <- lapply(mod.ests.lst(), function(ests.mat)
      ests.mat[!ests.mat[, "cnvg"], ])
    
    # Plot estimates from all models side-by-side
    
    # Set four plots per page
    par(mfrow = c(2, 2))
    
    # Population growth rate
    lambda <- input$lambda 
    
    # Load POPAN functions
    source("Functions/PopanFuncs.R", local = T)
    
    # Expected final population size.  gaps variables only used here
    lambda.gaps <- lambda^srvy.gaps # Population growth rate between surveys
    phi.gaps <- phi^srvy.gaps # Individual survival rate between surveys
    exp.N.fin <- sum(pent_func(lambda.gaps, phi.gaps) * exp.Ns * 
                       prod(phi.gaps) / cumprod(c(1, phi.gaps)))
    
    # Expected population size over time
    exp.N.t <- exp.N.fin / lambda^((hist.len - 1):0)

    # Load comparison plot function
    source("Functions/ComparisonPlot.R", local = T)
    
    # Plot estimates for lambda
    ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 1]),
                   "Lambda", lambda)
    
    # Plot estimates for Phi
    ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 2]),
                   "Phi", phi)
    
    # Plot estimates of N_2020
    ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 3]),
                   "N_2020", exp.N.t[hist.len])
    
    # Plot estimates of superpopulation size
    ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 4]),
                   "Ns", exp.Ns)
  })
  
  # # Plot first sample
  # output$scatterPlot <- renderPlot({
  #   switch (
  #     input$rv,
  #     "Bernoulli" = barplot(table(y()[, 1])),
  #     "Poisson" = barplot(table(y()[, 1])),
  #     "Normal" = hist(y()[, 1], main = "", xlab = "", ylab = "")
  #   )
  #   title(main = "Value counts for first sample", ylab = "Count", 
  #         xlab = "Value")
  # })
  # 
  # # Set parameter name for plot
  # par_name = reactive(switch(
  #   input$rv,
  #   "Bernoulli" = "Probability",
  #   "Poisson" = "Rate",
  #   "Normal" = "Mean"
  # ))
  # 
  # # Find maximum likelihood estimates and confidence intervals
  # mles = reactive(colMeans(y()))
  # cis = reactive({
  #   var_est = switch (
  #     input$rv,
  #     "Bernoulli" = mles() * (1 - mles()),
  #     "Poisson" = mles(),
  #     "Normal" = apply(y(), 2, var)
  #   )
  #   rbind(mles(), mles()) + c(-1, 1) * 1.96 * 
  #     sqrt(matrix(var_est, 2, n_samps, T) / n)
  # })
  # 
  # # Plot MLEs
  # output$MLEplot = renderPlot({
  #   boxplot(mles(), main = "Maximum likelihood estimates for all samples", 
  #           ylab = par_name())
  #   abline(h = true_val(), col = 2)
  # })
  # 
  # # Find parameter bounds
  # lb = reactive({
  #   switch(
  #     input$rv,
  #     "Bernoulli" = 0,
  #     "Poisson" = 0,
  #     "Normal" = -Inf,
  #   )
  # })
  # ub = reactive({
  #   switch(
  #     input$rv,
  #     "Bernoulli" = 1,
  #     "Poisson" = Inf,
  #     "Normal" = Inf,
  #   )
  # })
  # 
  # 
  # # Check CI doesn't cross parameter bounds
  # ci_ok = reactive(cis()[1, ] > lb() & cis()[2, ] < ub())
  # 
  # # Check CI coverage
  # ci_cov = reactive({
  #   ci_ok() & true_val() > cis()[1, ] & true_val() < cis()[2, ]
  # })
  # 
  # # Plot confidence intervals
  # output$CIPlot = renderPlot({
  #   plot(rep(1:n, 2), cis(), main = "Confidence intervals for all samples",
  #        ylab = par_name(), xlab = "Sample", type = 'n')
  #   arrows(1:n, cis()[1, ], 1:n, cis()[2, ], code = 3, length = 0.02, 
  #          angle = 90, lwd = 1 + !ci_cov())
  #   abline(h = true_val(), col = 2)
  #   abline(h = c(lb(), ub()))
  # })
  # 
  # # Print CI coverage
  # output$ciCov = renderText(
  #   paste("Confidence interval coverage over all samples:", mean(ci_cov()))
  # )
}