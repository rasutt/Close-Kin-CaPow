# Number of parameters
n.pars = reactive(length(est.par.names()))
# Number of models requested
n.mods = reactive(length(mdl.st()))


# When new datasets simulated
observeEvent(input$simulate, {
  # Nullify objects
  knshp.st(NULL)
  fit.lst(NULL)
  fll.SI.SY.IPs.lst(NULL)
  offst.SYIPs.lst(NULL)
  
  # If started multi-core cluster
  if (!is.null(cl())) {
    # Stop R sessions on other nodes
    stopCluster(cl())
  }
})

# When fit models button clicked
observeEvent(input$fit, {
  # Update reactive values for model and kinship sets
  mdl.st(input$mdl.st)
  knshp.st(input$knshp.st)
  
  # If fitting full genopair model for first time since simulation/loading
  if ("Full genopair" %in% mdl.st() && is.null(fll.SI.SY.IPs.lst())) {
    # Full set of sample-individual and sample-year index pairs for genopair
    # models, 2 x n_pairs x 2(??), representing individual and year that each
    # sample came from. Sample-years start from zero for TMB C++ objective
    # function
    fll.SI.SY.IPs.lst({
      # Make lists
      SIIPs.lst = SYIPs.lst = vector("list", n.sims())
      
      # Show progress-bar
      withProgress(
        # Loop over histories
        for (hst.ind in 1:n.sims()) {
          # Get simulated family and capture histories of population of
          # animals over time
          pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
          
          # Sample history matrix, n_individuals x n_surveys, rows ordered by
          # individual ID
          smp.hsts.bln = pop.cap.hist[, 4:(3 + k())] == 1
          
          # Sample-individual and sample-year index pairs
          SIIPs.lst[[hst.ind]] = combn(row(smp.hsts.bln)[smp.hsts.bln], 2)
          SYIPs.lst[[hst.ind]] = 
            t(combn(col(smp.hsts.bln)[smp.hsts.bln] - 1, 2))
          
          # Update progress-bar
          incProgress(1/n.sims())
        }, 
        value = 0, 
        message = "Finding all sample-individual and sample-year index pairs"
      )
      
      # Return lists
      list(SIIPs.lst = SIIPs.lst, SYIPs.lst = SYIPs.lst)
    })
    
    # Start new R sessions on separate "logical" CPU cores (nodes) for finding
    # genopair log-probabilities. Using 6 of my 12 cores seems to be optimal
    cl(makeCluster(6))
    
    # Evaluate reactive objects and pass values to new R sessions. Passing large
    # objects to parLapply does not improve performance (to my surprise). Can't
    # index in parLapply for some reason
    hst.lst.prll = sim.lst()$hists.lst
    SIIPs.lst.prll = fll.SI.SY.IPs.lst()$SIIPs.lst
    L.prll = L()
    k.prll = k()
    clusterExport(
      cl(),
      list(
        "hst.lst.prll", "SIIPs.lst.prll", "L.prll", "k.prll", "FindLogGPProbsKP"
      ), 
      environment()
    )
  }
  
  # If fitting offset genopair model for first time find offset survey-year
  # index pairs
  if (
    F
    # is.null(offst.SYIPs.lst()) &&
    # "Offset genopair" %in% mdl.st()
  ) {
    offst.SYIPs.lst.tmp = vector("list", n.sims())
    
    # Loop over histories
    withProgress({
      for (hst.ind in 1:n.sims()) {
        # Sample-year indices, columns in sample history matrix, representing
        # the survey that each sample came from, ordered by survey-year then
        # individual ID.  Counting from zero as will be passed to TMB objective
        # function
        smp.yr.inds = SYIs.lst()[[hst.ind]]
        
        # Sample index pairs for just consecutive pairs
        smp.ind.prs.offst = FindSIPsOffset(k(), smp.yr.inds)
        
        # Sample-year index pairs, n_pairs x 2, representing survey-year of each
        # sample in each offset pair, counting from zero as passing into TMB C++
        # objective function
        offst.SYIPs.lst.tmp[[hst.ind]] = 
          matrix(smp.yr.inds[as.vector(t(smp.ind.prs.offst))], ncol = 2)
        
        incProgress(1/n.sims())
      }
    }, value = 0, message = "Finding offset sample-year index pairs")
    
    offst.SYIPs.lst(offst.SYIPs.lst.tmp)
  }
  
  # Boolean for models requested 
  mod.bool = mdl.chcs %in% mdl.st()
  
  # Combine and keep only for selected models
  ests = list(fit.ppn()$ests, fit.ck()$ests, fit.gp()$ests, fit.os()$ests)[
    mod.bool
  ]
  ses = list(fit.ppn()$ses, fit.ck()$ses, fit.gp()$ses, fit.os()$ses)[mod.bool]
  cnvgs = list(fit.ppn()$cnvgs, fit.ck()$cnvgs, fit.gp()$cnvgs, fit.os()$cnvgs)[
    mod.bool
  ]
  names(ests) = names(ses) = names(cnvgs) = mdl.st()
  
  # Update list of model fits
  fit.lst(list(ests = ests, ses = ses, cnvgs = cnvgs))
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
  
  if(!is.null(fit.lst())) {
    # Loop over models requested
    for (i in 1:n.mods()) {
      # Find where model fit successfully
      ses.ok[[i]] = rowSums(is.na(fit.lst()$ses[[i]][, 1:4])) == 0
      cis.ok[[i]] = fit.lst()$cnvgs[[i]] & ses.ok[[i]]
      ests[[i]] = fit.lst()$ests[[i]][cis.ok[[i]], , drop = F]
      ses[[i]] = fit.lst()$ses[[i]][cis.ok[[i]], , drop = F]
      
      # Find differences between population parameter estimates and true values
      N.fin.errs[[i]] = ests[[i]][, 3] / sim.lst()$N.fin.vec[cis.ok[[i]]] - 1
      Ns.errs[[i]] = ests[[i]][, 4] / sim.lst()$Ns.vec[cis.ok[[i]]] - 1
      
      # Confidence intervals (creates matrices of correct size)
      radius = 1.96 * fit.lst()$ses[[i]]
      lcbs[[i]] = fit.lst()$ests[[i]] - radius
      ucbs[[i]] = fit.lst()$ests[[i]] + radius
      
      # Overwrite with log-normal CI's for population parameters
      l.vars = 
        log(1 + (fit.lst()$ses[[i]][, 3:4] / fit.lst()$ests[[i]][, 3:4])^2)
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
        true.pops > lcbs[[i]][, 3:4] & true.pops < ucbs[[i]][, 3:4] & 
        cis.ok[[i]]
      
      # Find proportions
      prpn.ci.cov[i, ] = colMeans(ci.cov[[i]])
      prpn.cnvgd[i] = mean(fit.lst()$cnvgs[[i]])
      prpn.ses.ok[i] = mean(ses.ok[[i]])
      prpn.cis.ok[i] = mean(cis.ok[[i]])
    }
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

# Show number of datasets models fit to
output$knshpSt = renderTable({
  if(!is.null(knshp.st())) {
    df = data.frame(knshp.st())
    names(df) = "Kinships included"
    df
  }
})

# Show convergence, standard error acceptance, and confidence interval
# acceptance rates for all models requested
output$mdlFtRts = renderTable({
  perc = function(stat) paste0(round(stat * 100, 1), "%")
  df = data.frame(
    mdl.st(), 
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

# Find CV and bias for estimates
est.bias.cv = reactive({
  # Make matrices
  ests.bias = ests.cv = matrix(NA, n.mods(), n.pars())
  
  # Find empirical CVs and biases
  ests.lst = check.ests()$ests
  ests.means = sapply(ests.lst, colMeans)
  ests.sds = sapply(ests.lst, apply, 2, sd)
  
  # Loop over parameters
  for (i in 1:n.pars()) {
    true.val = par.vals()[i]
    ests.mean = ests.means[i, ]
    ests.sd = ests.sds[i, ]
    
    # If the true value for this parameter is zero
    if (true.val == 0) {
      # Find proportional differences
      ests.bias[, i] <- ests.mean * 100
      ests.cv[, i] <- ests.sd * 100
    } else {
      # Find standard estimates
      ests.bias[, i] <- (ests.mean - true.val) / true.val * 100
      ests.cv[, i] <- ests.sd / ests.mean * 100
    }
  }
  
  # Return as list
  list(bias = ests.bias, cv = ests.cv)
})

# Function to make estimate performance table for output
makeEstPrfTbl = function(vals) {
  if (length(check.ests()$ests) > 0 && !is.null(nrow(check.ests()$ests[[1]]))) {
    df = data.frame(cbind(
      mdl.st(), matrix(paste0(round(vals, 1), "%"), n.mods())
    ))
    names(df) = c("Model", par.names())
    df
  }
}

# Show estimator performance in terms of bias and coefficient of variation
output$estBias = renderTable(makeEstPrfTbl(est.bias.cv()$bias))

# Show estimator performance in terms of bias and coefficient of variation
output$estCV = renderTable(makeEstPrfTbl(est.bias.cv()$cv))

# Function to compare estimates from POPAN and close kin models
ComparisonPlot <- function(ests.lst, par, true.val) {
  # Plot estimates
  boxplot(ests.lst, main = par, show.names = T)
  abline(h = true.val, col = 'red')
}

# Plot estimates using model comparison plot function
output$mdlCmpnPlt <- renderPlot({
  # If any estimates OK
  if(any(check.ests()$prpn.cis.ok > 0)) {
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
  }
})

# Print CI coverage
output$CICov = renderTable({
  df = data.frame(
    mdl.st(), matrix(perc(check.ests()$prpn.ci.cov), n.mods(), n.pars())
  )
  names(df) = c("Model", est.par.names())
  df
})

# Plot confidence intervals for lambda
output$CIPlot = renderPlot({
  if(!is.null(fit.lst())) {
    # Set space for multiple plots and set their margins
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
            main = mdl.st()[m], ylab = est.par.names()[p], xlab = "", 
            type = 'n'
          )
        } 
        # Otherwise plot CIs
        else {
          # Setup plot
          plot(
            rep(1:n.sims(), 2), 
            c(check.ests()$lcbs[[m]][, p], check.ests()$ucbs[[m]][, p]), 
            main = mdl.st()[m], ylab = est.par.names()[p], xlab = "", 
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
  }
})