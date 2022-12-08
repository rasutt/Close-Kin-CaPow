# Number of parameters
n.pars = reactive(length(est.par.names()))
# Names of models requested
mod.names = reactive({
  mod.choices[c(input$popan, input$close.kin, input$genopair)]
})
# Number of models requested
n.mods = reactive(input$popan + input$close.kin + input$genopair)

# Fit close-kin model
fit.gp = reactive(if (input$genopair) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), NA)
  ck.lwr <- c(0, 0.75, NA)
  ck.upr <- c(0.35, 1, Inf)
  
  # Create vector for model convergences
  gp.tmb.cnvg <- numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  gp.tmb.ests <- gp.tmb.ses <- matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hist.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hist.ind]]
      
      # Get numbers of animals captured in each survey
      ns.caps <- attributes(pop.cap.hist)$ns.caps
      
      # Update optimiser starting-values and bounds
      ck.start[3] <- attributes(pop.cap.hist)$N.t.vec[hist.len()]
      ck.lwr[3] <- ns.caps[k()]
      
      smp.hsts = as.matrix(pop.cap.hist[, 4:(3 + k())])
      smp.inds = row(smp.hsts)[as.logical(smp.hsts)]
      smp.gts = attributes(pop.cap.hist)$unq.smp.gts[, , smp.inds]
      n.smps = length(smp.inds)
      
      # 2 x n_pairs matrix of indices of samples in each pair to include in
      # likelihood, possibly all pairs or just consecutive pairs
      # smp.pr.inds = combn(n.smps, 2)
      # smp.pr.inds = rbind(1:n.smps, c(2:n.smps, 1))
      smp.pr.inds = rbind(1:(n.smps - 1), 2:n.smps)
      
      # Frequencies of 1-coded SNP alleles are means over both alleles at each
      # locus for each sample
      ale.frqs.1 = apply(attributes(pop.cap.hist)$unq.smp.gts, 2, mean)
      # Combine with frequencies for 0-coded alleles
      ale.frqs = rbind(1 - ale.frqs.1, ale.frqs.1)
      # Possible genopair probabilities at each locus given each kinship
      pss.gp.prbs.KPs = find.pss.gp.prbs.KPs(pss.gts, n.pss.gts, ale.frqs, L())
      # Genopair log-probabilities over all loci given each kinship, for each
      # pair to include in likelihood
      lg.gp.prbs.KPs = 
        FindLogGPProbsKP(smp.gts, smp.pr.inds, L(), pss.gp.prbs.KPs)
      
      # Get genopair probabilities (by excluding probabilities given half-sibs
      # for now) and check for pairs where all probabilities underflow to zero
      lg.gpp.slct = lg.gp.prbs.KPs[, -2]
      gpp.slct = exp(lg.gpp.slct)
      all_undrflw = rowSums(gpp.slct) == 0
      
      # If there is underflow adjust log-probabilities by factor giving equal
      # weight to smallest and largest values to avoid both under and overflow
      if (any(all_undrflw)) {
        cat("Proportion of pairs for which probabilities given all kinships 
        underflow to zero:", mean(all_undrflw), "\n")
        
        # Want smallest maximum kinship probability and largest probability to
        # be equally far from one
        adj = mean(c(min(apply(lg.gpp.slct, 1, max)), max(lg.gpp.slct)))
        lg.gpp.adj = lg.gpp.slct - adj
        gpp.adj = exp(lg.gpp.adj)
        
        # Show adjustment and results
        cat("Probabilities adjusted by factor of exp(", adj, ")\n", sep = "")
        print("Adjusted log-probabilities and probabilities of genopairs:")
        print(summary(lg.gpp.adj))
        print(summary(gpp.adj))
        cat("Proportion of pairs for which adjusted probabilities given all 
        kinships underflow to zero:", mean(rowSums(gpp.adj) == 0), "\n")
      }
      
      # Indices of survey-years for each sample
      smp.yr.inds = col(smp.hsts)[as.logical(smp.hsts)] - 1
      # n_pairs x 2 matrix of survey-years for each sample in each pair
      smp.yr.ind.prs = matrix(smp.yr.inds[as.vector(t(smp.pr.inds))], ncol = 2)
      
      print(table(smp.yr.ind.prs[, 1], smp.yr.ind.prs[, 2]))
      print(str(gpp.slct))
      
      # Try to fit genopair likelihood model
      gp.tmb.res = TryGenopairTMB(
        if (any(all_undrflw)) gpp.adj else gpp.slct, smp.yr.ind.prs, 
        k(), srvy.gaps(), fnl.year(), srvy.yrs(), ck.start, ck.lwr, ck.upr,
        alpha()
      )
      
      # If optimiser did not give error
      if(!all(is.na(gp.tmb.res))) {
        gp.tmb.ests[hist.ind, -(5:(4 + k()))] <- gp.tmb.res$est.se.df[, 1]
        gp.tmb.ses[hist.ind, -(5:(4 + k()))] <- gp.tmb.res$est.se.df[, 2]
        gp.tmb.cnvg[hist.ind] = gp.tmb.res$cnvg
      }
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting genopair model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = gp.tmb.ests, ses = gp.tmb.ses, cnvgs = !gp.tmb.cnvg)
})

# Fit close-kin model
fit.ck = reactive(if (input$close.kin) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), NA)
  ck.lwr <- c(0, 0.75, NA)
  ck.upr <- c(0.35, 1, Inf)
  
  # Create vector for model convergences
  ck.tmb.cnvg <- numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  ck.tmb.ests <- ck.tmb.ses <- matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n.sims()) {
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
        ns.kps.lst, k(), srvy.gaps(), fnl.year(), srvy.yrs(), ns.caps, 
        ck.start, ck.lwr, ck.upr, alpha()
      )
      ck.tmb.ests[hist.ind, -(5:(4 + k()))] <- ck.tmb.res$est.se.df[, 1]
      ck.tmb.ses[hist.ind, -(5:(4 + k()))] <- ck.tmb.res$est.se.df[, 2]
      ck.tmb.cnvg[hist.ind] = ck.tmb.res$cnvg
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting close-kin model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = ck.tmb.ests, ses = ck.tmb.ses, cnvgs = !ck.tmb.cnvg)
})

# Fit popan model
fit.ppn = reactive(if (input$popan) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), NA)
  ck.lwr <- c(0, 0.75, NA)
  ck.upr <- c(0.35, 1, Inf)
  ppn.start <- cbd.start <- c(ck.start, rep(p(), k()))
  ppn.lwr <- cbd.lwr <- c(ck.lwr, rep(0, k()))
  ppn.upr <- cbd.upr <- c(ck.upr, rep(1, k()))
  
  # Create vector model convergences
  ppn.tmb.cnvg <- numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  ppn.tmb.ests <- ppn.tmb.ses <- matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n.sims()) {
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
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting Popan model")
  
  # Combine model estimates, standard errors, and convergences, as lists and
  # return those requested
  list(ests = ppn.tmb.ests, ses = ppn.tmb.ses, cnvgs = !ppn.tmb.cnvg)
})

# Nullify model-fits when new datasets simulated
observeEvent(input$simulate, {
  fit.lst(NULL)
})

# Combine model estimates and update in fit.lst reactive value
observeEvent(input$nav.tab, {
  if (input$nav.tab == "model.tab" && is.null(fit.lst())) {
    # Boolean for models requested 
    mod.bool = c(input$popan, input$close.kin, input$genopair)
    
    fit.lst(list(
      ests = list(
        popan = fit.ppn()$ests, close.kin = fit.ck()$ests,
        genopair = fit.gp()$ests
      )[mod.bool],
      ses = list(
        popan = fit.ppn()$ses, close.kin = fit.ck()$ses,
        genopair = fit.gp()$ses
      )[mod.bool],
      cnvgs = list(
        popan = fit.ppn()$cnvgs, close.kin = fit.ck()$cnvgs,
        genopair = fit.gp()$cnvgs
      )[mod.bool]
    ))
  }
})

# Check when optimizer converged and standard errors calculable
check.ests = reactive({
  # Lists for retained model estimate and standard error matrices
  ests = ses = ses.ok = cis.ok = lcbs = ucbs = ci.cov = N.fin.errs = Ns.errs =
    list(n.mods())
  # Vectors for model stats
  prpn.cnvgd = prpn.ses.ok = prpn.cis.ok = numeric(n.mods())
  # Matrix for confidence interval coverage
  prpn.ci.cov = matrix(NA, n.mods(), n.pars())
  
  # Loop over models requested
  for (i in 1:n.mods()) {
    # Find where model fit successfully
    ses.ok[[i]] = rowSums(is.na(fit.lst()$ses[[i]][, 1:4])) == 0
    cis.ok[[i]] = fit.lst()$cnvgs[[i]] & ses.ok[[i]]
    ests[[i]] = fit.lst()$ests[[i]][cis.ok[[i]], ]
    ses[[i]] = fit.lst()$ses[[i]][cis.ok[[i]], ]
    
    # Find differences between population parameter estimates and true values
    N.fin.errs[[i]] = ests[[i]][, 3] / sim.lst()$N.fin.vec[cis.ok[[i]]] - 1
    Ns.errs[[i]] = ests[[i]][, 4] / sim.lst()$Ns.vec[cis.ok[[i]]] - 1
    
    # Confidence intervals (creates matrices of correct size)
    radius = 1.96 * fit.lst()$ses[[i]]
    lcbs[[i]] = fit.lst()$ests[[i]] - radius
    ucbs[[i]] = fit.lst()$ests[[i]] + radius
    
    # Overwrite with log-normal CI's for population parameters
    l.vars = log(1 + (fit.lst()$ses[[i]][, 3:4] / fit.lst()$ests[[i]][, 3:4])^2)
    fctr <- exp(1.959964 * sqrt(l.vars))
    lcbs[[i]][, 3:4] <- fit.lst()$ests[[i]][, 3:4] / fctr
    ucbs[[i]][, 3:4] <- fit.lst()$ests[[i]][, 3:4] * fctr
    
    # Bounds are matrices for studies x parameters, par.vals is a vector for
    # parameters, and cis.ok is a vector for studies
    ci.cov[[i]] = 
      t(par.vals() > t(lcbs[[i]]) & par.vals() < t(ucbs[[i]])) & cis.ok[[i]]
    
    # Overwrite for population parameters
    true.pops = cbind(sim.lst()$N.fin.vec, sim.lst()$Ns.vec)
    ci.cov[[i]][, 3:4] = 
      true.pops > lcbs[[i]][, 3:4] & true.pops < ucbs[[i]][, 3:4] & cis.ok[[i]]
    
    # Find proportions
    prpn.ci.cov[i, ] = colMeans(ci.cov[[i]])
    prpn.cnvgd[i] = mean(fit.lst()$cnvgs[[i]])
    prpn.ses.ok[i] = mean(ses.ok[[i]])
    prpn.cis.ok[i] = mean(cis.ok[[i]])
  }
  
  names(ests) = names(ses) = names(cis.ok) = names(lcbs) = names(ucbs) = 
    names(ci.cov) = names(N.fin.errs) = names(Ns.errs) = names(fit.lst()$ests)
  colnames(prpn.ci.cov) = est.par.names()
  
  list(
    ests = ests, ses = ses, lcbs = lcbs, ucbs = ucbs, ci.cov = ci.cov,
    ses.ok = ses.ok, cis.ok = cis.ok, prpn.ci.cov = prpn.ci.cov,
    prpn.cnvgd = prpn.cnvgd, prpn.ses.ok = prpn.ses.ok, 
    prpn.cis.ok = prpn.cis.ok, N.fin.errs = N.fin.errs, Ns.errs = Ns.errs
  )
})

# Show number of datasets models fit to
output$nDatasets = renderTable({
  df = data.frame(n.sims())
  names(df) = "Number of datasets"
  df
})

# Show convergence, standard error acceptance, and confidence interval
# acceptance rates for all models requested
output$modStats = renderTable({
  perc = function(stat) paste0(round(stat * 100, 1), "%")
  df = data.frame(
    mod.names(), 
    perc(check.ests()$prpn.cnvgd), 
    perc(check.ests()$prpn.ses.ok), 
    perc(check.ests()$prpn.cis.ok)
  )
  names(df) = c(
    "Model", 
    "Optimizer converged", 
    "Standard errors found", 
    "Fit successful"
  )
  df
})

# Plot estimates using model comparison plot function
output$modComp <- renderPlot({
  # Set four plots per figure
  par(mfrow = c(2, 2))
  
  # If any estimates OK
  if(check.ests()$prpn.cis.ok > 0) {
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
  }
})

# Print CI coverage
output$CICov = renderTable({
  df = data.frame(
    mod.names(), matrix(perc(check.ests()$prpn.ci.cov), n.mods(), n.pars())
  )
  names(df) = c("Model", est.par.names())
  df
})

# Print results for first study
output$firstResults <- renderTable({
  # True parameter values
  res.mat = matrix(par.vals(), nrow = 1)
  res.mat[1, 3:4] = c(sim.lst()$N.fin.vec[1], sim.lst()$Ns.vec[1])
  colnames(res.mat) = est.par.names()
  
  # Model estimates
  for (i in 1:n.mods()) {
    res.mat = rbind(
      res.mat, 
      fit.lst()$ests[[i]][1, ], fit.lst()$ses[[i]][1, ],
      check.ests()$lcbs[[i]][1, ], check.ests()$ucbs[[i]][1, ]
    )
  }
  
  # Result and model labels
  res.vals = c(
    "estimate", "standard_error", "lower_confidence_bound",
    "upper_confidence_bound"
  )
  res.df = data.frame(
    model = c("true_values", rep(mod.names(), each = 4)), 
    value = c("true_values", rep(res.vals, n.mods())),
    res.mat  
  )
  
  # Formatting
  res.df[, 5] = as.integer(res.df[, 5])
  res.df[, 6] = as.integer(res.df[, 6])
  res.df
}, digits = 3)

# Plot confidence intervals for lambda
output$CIPlot = renderPlot({
  par(mfrow = c(n.mods(), 1), mar = c(3.1, 4.1, 2.1, 2.1))
  # Loop over models requested
  for (m in 1:n.mods()) {
    # Loop over parameters, just lambda for now
    for (p in 1) {
      ord = order(fit.lst()$ests[[m]][, p])
      
      # If no valid estimates create empty plot
      if (all(is.na(fit.lst()$ses[[m]][, p]))) {
        plot(
          1:n.sims(), 
          rep(par.vals()[p], n.sims()), 
          main = mod.names()[m], ylab = est.par.names()[p], xlab = "", 
          type = 'n'
        )
      } 
      # Otherwise plot CIs
      else {
        # Setup plot
        plot(
          rep(1:n.sims(), 2), 
          c(check.ests()$lcbs[[m]][, p], check.ests()$ucbs[[m]][, p]), 
          main = mod.names()[m], ylab = est.par.names()[p], xlab = "", 
          type = 'n'
        )
      }
      # Plot estimates
      points(
        1:n.sims(), fit.lst()$ests[[m]][ord, p], pch = "-", 
        col = 1 + !check.ests()$ci.cov[[m]][ord, p]
      )
      # Plot intervals
      arrows(
        1:n.sims(), check.ests()$lcbs[[m]][ord, p], 
        1:n.sims(), check.ests()$ucbs[[m]][ord, p], 
        code = 3, length = 0.02, angle = 90, 
        col = 1 + !check.ests()$ci.cov[[m]][ord, p]
      )
      # True parameter value
      abline(h = par.vals()[p], col = 2)
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