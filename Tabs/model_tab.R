# Find parameter estimates, standard errors, and model convergences
fit.lst = reactive({
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), NA)
  ck.lwr <- c(0, 0.75, NA)
  ck.upr <- c(0.35, 1, Inf)
  ppn.start <- cbd.start <- c(ck.start, rep(p, k()))
  ppn.lwr <- cbd.lwr <- c(ck.lwr, rep(0, k()))
  ppn.upr <- cbd.upr <- c(ck.upr, rep(1, k()))
  
  # Create vectors for superpopulation and final population sizes, and model
  # convergences
  Ns.vec <- N.fin.vec <- ppn.tmb.cnvg <- ck.tmb.cnvg <- numeric(n_sims())
  
  # Create matrices for estimates and standard errors
  ppn.tmb.ests <- ck.tmb.ests <- ppn.tmb.ses <- ck.tmb.ses <- matrix(
    nrow = n_sims(), ncol = 4 + k(), 
    dimnames = list(
      NULL, c("lambda", "phi", "N_final", "Ns", paste0("p", 1:k()))
    )
  )
  
  # Boolean for models requested 
  mod.bool = c("POPAN", "Close kin") %in% models()
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n_sims()) {
      # Display progress
      cat("History:", hist.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- hists.lst()[[hist.ind]]
      
      # Store superpopulation and final population size
      Ns.vec[hist.ind] <- attributes(pop.cap.hist)$Ns
      N.fin.vec[hist.ind] <- attributes(pop.cap.hist)$N.t.vec[hist.len()]
      
      # Get numbers of animals captured in study and each survey
      n.cap.hists <- nrow(pop.cap.hist)
      ns.caps <- attributes(pop.cap.hist)$ns.caps
      
      # Summarise data for POPAN model
      pop.sum <- FindPopSum(k(), pop.cap.hist, n.cap.hists)
      
      # Find numbers of kin pairs
      ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
      
      # Update optimiser starting-values and bounds
      ppn.start[3] <- attributes(pop.cap.hist)$Ns
      ppn.lwr[3] <- n.cap.hists
      ck.start[3] <- N.fin.vec[hist.ind]
      ck.lwr[3] <- ns.caps[k()]
      
      # Try to fit models
      if (mod.bool[1]) {
        ppn.tmb.res <- TryPOPANTMB(
          k(), srvy.gaps(), n.cap.hists, pop.sum, ppn.start, ppn.lwr, ppn.upr
        )
        ppn.tmb.ests[hist.ind, ] <- ppn.tmb.res$est.se.df[, 1]
        ppn.tmb.ses[hist.ind, ] <- ppn.tmb.res$est.se.df[, 2]
        ppn.tmb.cnvg[hist.ind] = ppn.tmb.res$cnvg
      }
      if (mod.bool[2]) {
        ck.tmb.res <- TryCloseKinTMB(
          k(), srvy.gaps(), f.year(), srvy.yrs(), ns.caps, ns.kps.lst, 
          ck.start, ck.lwr, ck.upr
        )
        ck.tmb.ests[hist.ind, -(5:(4 + k()))] <- ck.tmb.res$est.se.df[, 1]
        ck.tmb.ses[hist.ind, -(5:(4 + k()))] <- ck.tmb.res$est.se.df[, 2]
        ck.tmb.cnvg[hist.ind] = ck.tmb.res$cnvg
      }
      
      incProgress(1/n_sims())
    }
  }, value = 0, message = "Fitting models")
  
  # Combine model estimates, standard errors, and convergences, as lists and
  # return those requested
  list(
    ests = list(popan = ppn.tmb.ests, close_kin = ck.tmb.ests)[mod.bool],
    ses = list(popan = ppn.tmb.ses, close_kin = ck.tmb.ses)[mod.bool],
    cnvgs = list(popan = !ppn.tmb.cnvg, close_kin = !ck.tmb.cnvg)[mod.bool]
  )
})

# Check when optimizer converged and standard errors calculable
check.ests = reactive({
  # Number of models requested
  n_mods = length(models())
  
  # Lists for retained model estimate and standard error matrices
  ests = ses = ses_ok = cis_ok = list(n_mods)
  
  # Vectors for model stats
  prpn_cnvgd = prpn_ses_ok = prpn_cis_ok = numeric(n_mods)
  
  # Loop over models requested
  for (i in 1:n_mods) {
    ses_ok[[i]] = rowSums(is.na(fit.lst()$ses[[i]][, 1:4])) == 0
    cis_ok[[i]] = fit.lst()$cnvgs[[i]] & ses_ok[[i]]
    ests[[i]] = fit.lst()$ests[[i]][cis_ok[[i]], ]
    ses[[i]] = fit.lst()$ses[[i]][cis_ok[[i]], ]
    prpn_cnvgd[i] = mean(fit.lst()$cnvgs[[i]])
    prpn_ses_ok[i] = mean(ses_ok[[i]])
    prpn_cis_ok[i] = mean(cis_ok[[i]])
  }
  names(ests) = names(ses) = names(fit.lst()$ests)
  list(
    ests = ests, ses = ses, ses_ok = ses_ok, cis_ok = cis_ok, 
    prpn_cnvgd = prpn_cnvgd, prpn_ses_ok = prpn_ses_ok, 
    prpn_cis_ok = prpn_cis_ok
  )
})
ests = reactive(check.ests()$ests)

# Show convergence, standard error acceptance, and confidence interval
# acceptance rates for all models requested
output$modStats = renderTable({
  perc = function(stat) paste0(round(stat * 100, 1), "%")
  data.frame(
    model = models(), 
    percent_converged = perc(check.ests()$prpn_cnvgd), 
    percent_SEs_OK = perc(check.ests()$prpn_ses_ok), 
    percent_CIs_OK = perc(check.ests()$prpn_cis_ok)
  )
})

cnvgd = reactive(fit.lst()$cnvgs$close_kin)
ses_ok = reactive(rowSums(is.na(fit.lst()$ses$close_kin[, 1:4])) == 0)
cis_ok = reactive(cnvgd() & ses_ok())

# Plot negative log-likelihood surface for first study
output$NLLPlot <- renderPlot({
  # Get simulated family and capture histories of population of animals over
  # time
  pop.cap.hist <- hists.lst()[[1]]
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

# Print first few estimates and convergences for first model
output$firstEsts <- renderTable({
  ests.cnvg = data.frame(
    round(fit.lst()$ests[[1]], 3), cnvgd = fit.lst()$cnvgs[[1]]
  )
  head(ests.cnvg)
})

# Plot estimates using model comparison plot function
output$modComp <- renderPlot({
  # Plot estimates from all models side-by-side
  
  # Set four plots per page
  par(mfrow = c(2, 2))
  
  # Plot estimates for lambda
  ComparisonPlot(
    lapply(ests(), function(ests.mat) ests.mat[, 1]), 
    "Population growth rate", lambda()
  )
  
  # Plot estimates for Phi
  ComparisonPlot(
    lapply(ests(), function(ests.mat) ests.mat[, 2]),
    "Survival rate", phi()
  )
  
  # Plot estimates of final population size
  ComparisonPlot(
    lapply(ests(), function(ests.mat) ests.mat[, 3]),
    "Final population size", exp.N.t()[hist.len()]
  )
  
  # Plot estimates of superpopulation size
  ComparisonPlot(
    lapply(ests(), function(ests.mat) ests.mat[, 4]),
    "Super-population size", exp.Ns
  )
})

# Find confidence intervals for close kin model
lcbs = reactive({
  fit.lst()$ests[["close_kin"]] - 1.96 * fit.lst()$ses[["close_kin"]]
})
ucbs = reactive({
  fit.lst()$ests[["close_kin"]] + 1.96 * fit.lst()$ses[["close_kin"]]
})

# # Check CI doesn't cross parameter bounds
# ci_ok = reactive(cis()[1, ] > lb() & cis()[2, ] < ub())

# Check CI coverage
ci_cov = reactive({
  # ci_ok() & true_val() > cis()[1, ] & true_val() < cis()[2, ]
  lambda() > lcbs()[, 1] & lambda() < ucbs()[, 1]
})

# Plot confidence intervals for close kin model for lambda
output$CIPlot = renderPlot({
  ord = order(fit.lst()$ests[["close_kin"]][, 1])
  plot(
    rep(1:n_sims(), 2), c(lcbs()[, 1], ucbs()[, 1]), 
    main = "Confidence intervals for close kin model for lambda",
    ylab = "Lambda", xlab = "", type = 'n'
  )
  arrows(
    1:n_sims(), lcbs()[, 1][ord], 
    1:n_sims(), ucbs()[, 1][ord], 
    code = 3, length = 0.02, angle = 90,
    lwd = 1 + !ci_cov()[ord]
  )
  abline(h = lambda(), col = 2)
  # abline(h = c(lb(), ub()))
})

# Print CI coverage
output$CICov = renderText({
  paste0(
    "95% confidence interval coverage for close kin model for lambda: ", 
    round(mean(ci_cov()) * 100, 1), "%"
  )
})