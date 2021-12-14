# Define server logic for app
server <- function(input, output) {
  # Population growth rate
  lambda <- reactive(input$lambda)
  
  # "Population birthrate"
  rho <- reactive(lambda() - phi)
  
  # Expected population size over time
  exp.N.t = reactive({
    # Load POPAN functions
    source("Functions/PopanFuncs.R", local = T)
    
    # Expected final population size.  gaps variables only used here
    lambda.gaps <- lambda()^srvy.gaps # Population growth rate between surveys
    phi.gaps <- phi^srvy.gaps # Individual survival rate between surveys
    exp.N.fin <- sum(pent_func(lambda.gaps, phi.gaps) * exp.Ns * 
                       prod(phi.gaps) / cumprod(c(1, phi.gaps)))
    
    # Expected population size over time
    exp.N.t <- exp.N.fin / lambda()^((hist.len - 1):0)
    
    exp.N.t
  })
  
  # Simulate population and capture histories
  hists.lst = reactive({
    # Population growth rate
    lambda = lambda()
    
    # Set number of studies to simulate
    n.stds.sim <- input$n_sims
    
    # Initial population size
    N.init <- round(exp.N.t()[1]) 
    
    # Load simulation and POPAN functions
    source("Functions/SimPopStud.R", local = T)
    source("Functions/PopanFuncs.R", local = T)
    
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
  
  # Check simulated studies
  checks.lst = reactive({
    lambda = lambda()
    rho = rho()
    
    # Set number of studies to check
    n.stds.chk <- input$n_sims
    
    # Expected final population size
    exp.N.fin <- exp.N.t()[hist.len]
    
    # Load kin pair functions
    source("Functions/FindNsKinPairs.R", local = T)
    source("Functions/FindExpNsKPs.R", local = T)
    
    # Create matrices for population trajectories (as columns), numbers
    # captured, and expected and observed number of kin pairs
    N.t.mat <- matrix(nrow = n.stds.chk, ncol = hist.len)
    ns.caps.mat <- ns.clvng.caps.mat <- ns.clvng.mat <- ns.POPs.wtn.mat <- 
      ns.HSPs.wtn.mat <- exp.ns.HSPs.wtn.mat <- 
      exp.ns.POPs.wtn.mat <- matrix(nrow = n.stds.chk, ncol = k)
    ns.POPs.btn.mat <- ns.SPs.btn.mat <- exp.ns.POPs.btn.mat <- 
      exp.ns.SPs.btn.mat <- matrix(nrow = n.stds.chk, ncol = n.srvy.prs)
    
    # Create vectors for proportions with unknown parents, and superpopulation
    # sizes
    prpn.prnts.unkn.vec <- Ns.vec <- numeric(n.stds.chk)
    
    # Display progress
    cat("Checking study: ")
    
    # Loop over histories
    for (hist.ind in 1:n.stds.chk) {
      # Display progress
      if (hist.ind %% 100 == 1) cat(hist.ind, "")
      
      # Get simulated family and capture histories of population of animals over
      # time
      pop.cap.hist <- hists.lst()[[hist.ind]]
      
      # Record population curve
      N.t.mat[hist.ind, ] <- attributes(pop.cap.hist)$N.t.vec
      
      # Record superpopulation size
      Ns.vec[hist.ind] <- attributes(pop.cap.hist)$Ns
      
      # Get numbers captured and calving in each survey
      ns.caps <- attributes(pop.cap.hist)$ns.caps
      ns.caps.mat[hist.ind, ] <- ns.caps
      ns.clvng.mat[hist.ind, ] <- attributes(pop.cap.hist)$ns.clvng
      ns.clvng.caps.mat[hist.ind, ] <- 
        colSums(pop.cap.hist[, 4:(3 + k)] * pop.cap.hist[, (4 + k):(3 + 2 * k)])
      
      # Find proportion captured with unknown parents
      prpn.prnts.unkn.vec[hist.ind] <- mean(is.na(pop.cap.hist$mum))
      
      # Find numbers of known kin pairs
      ns.kps.lst <- FindNsKinPairs()
      
      # Record in matrices
      ns.POPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.wtn
      ns.HSPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.HSPs.wtn
      ns.POPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.btn
      ns.SPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.SPs.btn
      
      # Find expected numbers of kin pairs
      exp.ns.kps.lst <- FindExpNsKPs()
      
      # Record in matrices
      exp.ns.POPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.wtn
      exp.ns.HSPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.HSPs.wtn
      exp.ns.POPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.btn
      exp.ns.SPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.SPs.btn
    }
    
    list(
      N.t.mat = N.t.mat, 
      ns.POPs.wtn.mat = ns.POPs.wtn.mat,
      exp.ns.POPs.wtn.mat = exp.ns.POPs.wtn.mat,
      prpn.prnts.unkn.vec = prpn.prnts.unkn.vec,
      ns.caps.mat = ns.caps.mat
    )
  })
  
  # Plot population sizes over time
  output$popPlot <- renderPlot({
    # Plot population trajectories and expected value
    matplot(
      (f.year - 79):f.year, t(checks.lst()$N.t.mat), type = 'l', 
      col = rgb(0, 0, 0, alpha = 0.1), lty = 1, 
      xlab = 'Year', ylab = 'Nt', main = "Population sizes over time"
    )
    lines((f.year - 79):f.year, exp.N.t(), col = 'red', lwd = 2)
    
    # Surveys
    abline(v = srvy.yrs, lty = 2)
    
    # Add legend
    legend(
      "topleft", 
      legend = c("Population sizes", "Expected population size", 
                 "Survey years"),
      col = c(rgb(0, 0, 0, alpha = 0.1), 2, 1),
      lwd = c(1, 2, 1),
      lty = c(1, 1, 2)
    )
  })
  
  # Print head of first study
  output$dataHead <- renderTable(head(data.frame(hists.lst()[[1]])))
  
  # Plot numbers of parent-offspring pairs within samples
  output$nPOPsPlot <- renderPlot({
    # srvy.inds = hist.len + srvy.yrs - f.year
    # prpns = checks.lst()$ns.POPs.wtn.mat / choose(checks.lst()$ns.caps.mat, 2)
    # exp.prpns = 
    #   checks.lst()$exp.ns.POPs.wtn.mat / choose(exp.N.t()[srvy.inds] * p, 2)
    diffs = checks.lst()$ns.POPs.wtn.mat - checks.lst()$exp.ns.POPs.wtn.mat
    boxplot(
      diffs,
      main = "Expected vs observed numbers of parent-offspring pairs within 
      samples", 
      # sub = paste(
      #   "Average difference over all surveys:",
      #   signif(mean(diffs), 3)
      # ),
      xlab = "Survey", 
      ylab = "Observed - expected numbers"
    )
    abline(h = 0, col = 'red')
    abline(h = mean(diffs), col = 'blue')
    legend(
      "topleft", 
      legend = c(
        "Expected difference over all surveys (zero)", 
        "Average difference over all surveys"
      ),
      col = c(2, 4),
      lty = 1
    )
  })
  
  # # Display proportion of captures for which the parents are unknown
  # output$nUnknPrnts <- renderText({
  #   nPOPs = checks.lst()$ns.POPs.wtn.mat
  # 
  #   # Pairs are lost quadratically with animals
  #   prpns.pairs.lost = 1 - (1 - checks.lst()$prpn.prnts.unkn.vec)^2
  # 
  #   paste(
  #     "Expected number of parent-offspring pairs within samples 
  #     lost due to unknown parents:",
  #     signif(mean(prpns.pairs.lost * checks.lst()$ns.POPs.wtn.mat), 3),
  #     "\n"
  #   )
  # })
  
  # Find parameter estimates
  mod.ests.lst = reactive({
    # Set number of studies to fit models to (reduce for faster testing)
    n.stds.fit <- input$n_sims
    
    # Load functions
    funcs <- list.files("Functions")
    for (i in 1:length(funcs)) source(paste0("Functions/", funcs[i]), local = T)
    
    # Create general optimizer starting-values and bounds, NAs filled in below
    ck.start <- c(rho(), phi, NA)
    ck.lwr <- c(0, 0.75, NA)
    ck.upr <- c(0.35, 1, Inf)
    
    # Create vectors for superpopulation and final population sizes
    Ns.vec <- N.fin.vec <- numeric(n.stds.fit)
    
    # Create matrices for estimates
    ck.ests <- ck.tmb.ests <- matrix(
      nrow = n.stds.fit, ncol = 5 + k, dimnames = list(
        NULL, c("lambda", "phi", "N_final", "Ns", paste0("p", 1:k), "cnvg")
      )
    )
    
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
      print(n.cap.hists)
      
      # Find numbers of kin pairs
      ns.kps.lst <- FindNsKinPairs()
      
      # Update optimiser starting-values and bounds
      ck.start[3] <- N.fin.vec[hist.ind]
      ck.lwr[3] <- ns.caps[k]

      # Try to fit models
      # ck.ests[hist.ind, -(4:(4 + k))] <- TryCloseKin()
      ck.tmb.ests[hist.ind, -(5:(4 + k))] <- TryCloseKinTMB()
    }
    
    # # Find close kin estimates of Ns without TMB
    # ck.ests[, 4] <- Ns_vec_func(ck.ests[, 3], ck.ests[, 1], ck.ests[, 2])
    
    # Combine model estimates as list
    mod.ests.lst <- list(
      # ck = ck.ests
      ck.tmb = ck.tmb.ests
    )  
    mod.ests.lst
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
      main = "Negative log-likelihood at MLEs for first study",
      xlab = "lambda", ylab = "NLL", type = 'l'
    )
    abline(v = lambda(), col = 2)
    abline(v = mod.ests.lst()[[1]][1, 1], col = 4)
    legend(
      "topleft", 
      legend = c("Negative log likelihood", "True lambda", "MLE of lambda"),
      col = c(1, 2, 4),
      lty = 1
    )
    # abline(v = cis()[1, 1], col = 4, lty = 2)
    # abline(v = cis()[2, 1], col = 4, lty = 2)
  })
  
  # Print first few estimates
  output$firstEsts <- renderTable(
    lapply(mod.ests.lst(), function(mod.ests) head(round(mod.ests, 3)))[[1]]
  )
  
  # Plot estimates using model comparison plot function (though only one model
  # atm)
  output$modComp <- renderPlot({
    # Expected population size over time
    exp.N.t <- exp.N.t()
    
    # Find optimization attempts that converged
    cvgd.ests.lst <- lapply(mod.ests.lst(), function(ests.mat)
      ests.mat[!ests.mat[, "cnvg"], ])
    
    # Plot estimates from all models side-by-side
    
    # Set four plots per page
    par(mfrow = c(2, 2))

    # Load comparison plot function
    source("Functions/ComparisonPlot.R", local = T)
    
    # Plot estimates for lambda
    ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 1]),
                   "Lambda", lambda())
    
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