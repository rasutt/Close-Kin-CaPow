# True parameter values
true_vals = reactive({
  c(lambda(), phi(), exp.N.t()[hist.len()], exp.Ns(), rep(p, k()))
})
# Parameter names
par_names = reactive(c("lambda", "phi", "N_final", "Ns", paste0("p", 1:k())))
# Number of parameters
n_pars = reactive(length(par_names()))
# Names of models requested
mod_names = reactive(mod_choices[c(input$popan, input$close_kin)])
# Number of models requested
n_mods = reactive(input$popan + input$close_kin)
# Function to prepare proportion to print as percentage
perc = function(prpn) paste0(round(prpn * 100, 1), "%")

# Fit close-kin model
fit.ck = reactive(if (input$close_kin) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), NA)
  ck.lwr <- c(0, 0.75, NA)
  ck.upr <- c(0.35, 1, Inf)
  
  # Create vector for model convergences
  ck.tmb.cnvg <- numeric(n_sims())
  
  # Create matrices for estimates and standard errors
  ck.tmb.ests <- ck.tmb.ses <- 
    matrix(nrow = n_sims(), ncol = 4 + k(), dimnames = list(NULL, par_names()))
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n_sims()) {
      # Display progress
      cat("History:", hist.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hist.ind]]
      
      # Get numbers of animals captured in each survey
      ns.caps <- attributes(pop.cap.hist)$ns.caps
      
      # Find numbers of kin pairs
      ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
      
      # Update optimiser starting-values and bounds
      ck.start[3] <- attributes(pop.cap.hist)$N.t.vec[hist.len()]
      ck.lwr[3] <- ns.caps[k()]
      
      # Try to fit model
      ck.tmb.res <- TryCloseKinTMB(
        k(), srvy.gaps(), f.year(), srvy.yrs(), ns.caps, ns.kps.lst, 
        ck.start, ck.lwr, ck.upr
      )
      ck.tmb.ests[hist.ind, -(5:(4 + k()))] <- ck.tmb.res$est.se.df[, 1]
      ck.tmb.ses[hist.ind, -(5:(4 + k()))] <- ck.tmb.res$est.se.df[, 2]
      ck.tmb.cnvg[hist.ind] = ck.tmb.res$cnvg
      
      incProgress(1/n_sims())
    }
  }, value = 0, message = "Fitting models")
  
  # Combine model estimates, standard errors, and convergences, as lists and
  # return those requested
  list(ests = ck.tmb.ests, ses = ck.tmb.ses, cnvgs = !ck.tmb.cnvg)
})

# Fit popan model
fit.ppn = reactive(if (input$popan) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), NA)
  ck.lwr <- c(0, 0.75, NA)
  ck.upr <- c(0.35, 1, Inf)
  ppn.start <- cbd.start <- c(ck.start, rep(p, k()))
  ppn.lwr <- cbd.lwr <- c(ck.lwr, rep(0, k()))
  ppn.upr <- cbd.upr <- c(ck.upr, rep(1, k()))
  
  # Create vector model convergences
  ppn.tmb.cnvg <- numeric(n_sims())
  
  # Create matrices for estimates and standard errors
  ppn.tmb.ests <- ppn.tmb.ses <- 
    matrix(nrow = n_sims(), ncol = 4 + k(), dimnames = list(NULL, par_names()))
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n_sims()) {
      # Display progress
      cat("History:", hist.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hist.ind]]
      
      # Get numbers of animals captured in study
      n.cap.hists <- nrow(pop.cap.hist)
      
      # Summarise data for POPAN model
      pop.sum <- FindPopSum(k(), pop.cap.hist, n.cap.hists)
      
      # Update optimiser starting-values and bounds
      ppn.start[3] <- attributes(pop.cap.hist)$Ns
      ppn.lwr[3] <- n.cap.hists
      
      # Try to fit models
      ppn.tmb.res <- TryPOPANTMB(
        k(), srvy.gaps(), n.cap.hists, pop.sum, ppn.start, ppn.lwr, ppn.upr
      )
      ppn.tmb.ests[hist.ind, ] <- ppn.tmb.res$est.se.df[, 1]
      ppn.tmb.ses[hist.ind, ] <- ppn.tmb.res$est.se.df[, 2]
      ppn.tmb.cnvg[hist.ind] = ppn.tmb.res$cnvg
      
      incProgress(1/n_sims())
    }
  }, value = 0, message = "Fitting models")
  
  # Combine model estimates, standard errors, and convergences, as lists and
  # return those requested
  list(ests = ppn.tmb.ests, ses = ppn.tmb.ses, cnvgs = !ppn.tmb.cnvg)
})

# Combine model estimates
fit.lst = reactive({
  # Boolean for models requested 
  mod.bool = c(input$popan, input$close_kin)
  
  list(
    ests = list(popan = fit.ppn()$ests, close_kin = fit.ck()$ests)[mod.bool],
    ses = list(popan = fit.ppn()$ses, close_kin = fit.ck()$ses)[mod.bool],
    cnvgs = list(popan = fit.ppn()$cnvgs, close_kin = fit.ck()$cnvgs)[mod.bool]
  )
})

# Check when optimizer converged and standard errors calculable
check.ests = reactive({
  # Lists for retained model estimate and standard error matrices
  ests = ses = ses_ok = cis_ok = lcbs = ucbs = ci_cov = N.fin.errs = Ns.errs =
    list(n_mods())
  # Vectors for model stats
  prpn_cnvgd = prpn_ses_ok = prpn_cis_ok = numeric(n_mods())
  # Matrix for confidence interval coverage
  prpn_ci_cov = matrix(NA, n_mods(), n_pars())
  
  # Loop over models requested
  for (i in 1:n_mods()) {
    # Find where model fit successfully
    ses_ok[[i]] = rowSums(is.na(fit.lst()$ses[[i]][, 1:4])) == 0
    cis_ok[[i]] = fit.lst()$cnvgs[[i]] & ses_ok[[i]]
    ests[[i]] = fit.lst()$ests[[i]][cis_ok[[i]], ]
    ses[[i]] = fit.lst()$ses[[i]][cis_ok[[i]], ]
    
    # Find differences between population parameter estimates and true values
    N.fin.errs[[i]] = ests[[i]][, 3] / sim.lst()$N.fin.vec[cis_ok[[i]]] - 1
    Ns.errs[[i]] = ests[[i]][, 4] / sim.lst()$Ns.vec[cis_ok[[i]]] - 1

    # Confidence intervals (creates matrices of correct size)
    radius = 1.96 * fit.lst()$ses[[i]]
    lcbs[[i]] = fit.lst()$ests[[i]] - radius
    ucbs[[i]] = fit.lst()$ests[[i]] + radius
    
    # Overwrite with log-normal CI's for population parameters
    l_vars = log(1 + (fit.lst()$ses[[i]][, 3:4] / fit.lst()$ests[[i]][, 3:4])^2)
    fctr <- exp(1.959964 * sqrt(l_vars))
    lcbs[[i]][, 3:4] <- fit.lst()$ests[[i]][, 3:4] / fctr
    ucbs[[i]][, 3:4] <- fit.lst()$ests[[i]][, 3:4] * fctr
    
    # Bounds are matrices for studies x parameters, true_vals is a vector for
    # parameters, and cis_ok is a vector for studies
    ci_cov[[i]] = 
      t(true_vals() > t(lcbs[[i]]) & true_vals() < t(ucbs[[i]])) & cis_ok[[i]]
    
    # Overwrite for population parameters
    true_pops = cbind(sim.lst()$N.fin.vec, sim.lst()$Ns.vec)
    ci_cov[[i]][, 3:4] = 
      true_pops > lcbs[[i]][, 3:4] & true_pops < ucbs[[i]][, 3:4] & cis_ok[[i]]
    
    # Find proportions
    prpn_ci_cov[i, ] = colMeans(ci_cov[[i]])
    prpn_cnvgd[i] = mean(fit.lst()$cnvgs[[i]])
    prpn_ses_ok[i] = mean(ses_ok[[i]])
    prpn_cis_ok[i] = mean(cis_ok[[i]])
  }
  
  names(ests) = names(ses) = names(cis_ok) = names(lcbs) = names(ucbs) = 
    names(ci_cov) = names(N.fin.errs) = names(Ns.errs) = names(fit.lst()$ests)
  colnames(prpn_ci_cov) = par_names()
  
  list(
    ests = ests, ses = ses, lcbs = lcbs, ucbs = ucbs, ci_cov = ci_cov,
    ses_ok = ses_ok, cis_ok = cis_ok, prpn_ci_cov = prpn_ci_cov,
    prpn_cnvgd = prpn_cnvgd, prpn_ses_ok = prpn_ses_ok, 
    prpn_cis_ok = prpn_cis_ok, N.fin.errs = N.fin.errs, Ns.errs = Ns.errs
  )
})

# Show convergence, standard error acceptance, and confidence interval
# acceptance rates for all models requested
output$modStats = renderTable({
  perc = function(stat) paste0(round(stat * 100, 1), "%")
  data.frame(
    model = mod_names(), 
    optimizer_converged = perc(check.ests()$prpn_cnvgd), 
    standard_errors_found = perc(check.ests()$prpn_ses_ok), 
    fit_successful = perc(check.ests()$prpn_cis_ok)
  )
})

# Print first few estimates for each model
output$firstEsts <- renderTable({
  # Setup matrix and true parameter values
  rows = matrix(NA, n_mods() + 1, n_pars())
  rows[1, ] = true_vals()
  rows[1, 3:4] = c(sim.lst()$N.fin.vec[1], sim.lst()$Ns.vec[1])
  colnames(rows) = par_names()
  
  # Add model estimates
  for (i in 1:n_mods()) {
    rows[i + 1, ] = fit.lst()$ests[[i]][1, ]
  }
  
  # Format and output table
  ests.cnvg = data.frame(rows, row.names = c("True values", mod_names()))
  ests.cnvg[, 4] = as.integer(ests.cnvg[, 4])
  ests.cnvg[, 5] = as.integer(ests.cnvg[, 5])
  ests.cnvg
}, digits = 3, rownames = T)

# Plot estimates using model comparison plot function
output$modComp <- renderPlot({
  # Set four plots per figure
  par(mfrow = c(2, 2))
  
  # Plot estimates from all models side-by-side
  ComparisonPlot(
    lapply(check.ests()$ests, function(ests.mat) ests.mat[, 1]), 
    "Population growth rate", lambda()
  )
  ComparisonPlot(
    lapply(check.ests()$ests, function(ests.mat) ests.mat[, 2]),
    "Survival rate", phi()
  )
  ComparisonPlot(
    check.ests()$N.fin.errs, "Final population size proportional errors", 0
  )
  ComparisonPlot(
    check.ests()$Ns.errs, "Super-population size proportional errors", 0
  )
})

# Print CI coverage
output$CICov = renderTable({
  data.frame(
    model = mod_names(), 
    matrix(
      perc(check.ests()$prpn_ci_cov), n_mods(), n_pars(), 
      dimnames = list(NULL, par_names())
    )
  )
})

# Plot confidence intervals for lambda
output$CIPlot = renderPlot({
  par(mfrow = c(n_mods(), 1), mar = c(3.1, 4.1, 2.1, 2.1))
  # Loop over models requested
  for (m in 1:n_mods()) {
    # Loop over parameters, just lambda for now
    for (p in 1) {
      ord = order(fit.lst()$ests[[m]][, p])
      # Setup plot
      plot(
        rep(1:n_sims(), 2), 
        c(check.ests()$lcbs[[m]][, p], check.ests()$ucbs[[m]][, p]), 
        main = mod_names()[m], ylab = par_names()[p], xlab = "", type = 'n'
      )
      # Plot estimates
      points(
        1:n_sims(), fit.lst()$ests[[m]][ord, p], pch = "-", 
        col = 1 + !check.ests()$ci_cov[[m]][ord, p]
      )
      # Plot intervals
      arrows(
        1:n_sims(), check.ests()$lcbs[[m]][ord, p], 
        1:n_sims(), check.ests()$ucbs[[m]][ord, p], 
        code = 3, length = 0.02, angle = 90, 
        col = 1 + !check.ests()$ci_cov[[m]][ord, p]
      )
      # True parameter value
      abline(h = true_vals()[p], col = 2)
      # abline(h = c(lb(), ub()))
      # legend(
      #   "topleft", col = 1:2, lwd = c(1, 1), 
      #   legend = c(
      #     "Optimizer converged and CI covers true value", 
      #     "Did not converge and/or does not cover"
      #   )
      # )
    }
  }
})

# Print first few estimates and convergences for first model
output$firstCIs <- renderTable({
  ests.cnvg = data.frame(
    rbind(
      fit.lst()$ests[[1]][1, ], fit.lst()$ses[[1]][1, ],
      check.ests()$lcbs[[1]][1, ], check.ests()$ucbs[[1]][1, ]
    ), 
    row.names = c(
      "estimate", "standard_error", "lower_confidence_bound",
      "upper_confidence_bound"
    )
  )
}, digits = 3, rownames = T)

# Plot negative log-likelihood surface for first study
output$NLLPlot <- renderPlot({
  # Get simulated family and capture histories of population of animals over
  # time
  pop.cap.hist <- sim.lst()$hists.lst[[1]]
  # Get numbers of animals captured in each survey
  ns.caps <- attributes(pop.cap.hist)$ns.caps
  # Find numbers of kin pairs
  ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
  # MLEs for rho, phi, and Ns
  params = fit.lst()$ests[["close_kin"]][1, 1:3]
  # Create grids of rho and NLL values
  nll_grid = rho_grid = seq(min_rho, max_rho, step_rho) + phi() - params[2]
  
  # Find NLL over grid of rho values
  for (i in seq_along(nll_grid)) {
    params[1] = rho_grid[i]
    nll_grid[i] = CloseKinNLL(
      params, k(), f.year(), srvy.yrs(), ns.kps.lst, ns.caps
    )
  }
  
  # Plot NLL over grid of lambda values
  plot(
    rho_grid + params[2], nll_grid,
    main = "Negative log-likelihood at MLEs from close kin model
      for first study",
    xlab = "Lambda", ylab = "NLL", type = 'l'
  )
  abline(v = lambda(), col = 2)
  abline(v = fit.lst()$ests[["close_kin"]][1, 1], col = 4)
  legend(
    "topleft", 
    legend = c("Negative log likelihood", "True lambda", "MLE of lambda"),
    col = c(1, 2, 4),
    lty = 1
  )
  # abline(v = cis()[1, 1], col = 4, lty = 2)
  # abline(v = cis()[2, 1], col = 4, lty = 2)
})