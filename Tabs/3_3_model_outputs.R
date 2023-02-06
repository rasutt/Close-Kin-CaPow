# Number of parameters
n.pars = reactive(length(est.par.names()))
# Number of models requested
n.mods = reactive(length(mdl.st()))

# When new datasets simulated
observeEvent(input$simulate, {
  # Nullify objects
  knshp.st(NULL)
  fit.lst(NULL)
  osisyips.lst(NULL)
  
  # If started multi-core cluster
  if (!is.null(cl())) {
    # Stop R sessions on other nodes
    stopCluster(cl())
    
    # Nullify cluster reactive value
    cl(NULL)
  }
})

# When fit models button clicked
observeEvent(input$fit, {
  # Update reactive values for model and kinship sets
  mdl.st(input$mdl.st)
  knshp.st(input$knshp.st)
  
  # If the numbers of studies and/or loci are large, and multicore cluster not
  # setup for current dataset, and full genopair model being fit
  if (n.sims() * L() > 1e3 && 
      is.null(cl()) && 
      "Full genopair" %in% mdl.st()
  ) {
    # Start new R sessions on separate "logical" CPU cores (nodes) for finding
    # genopair log-probabilities. Using 6 of my 12 cores seems to be optimal
    cl.tmp = makeCluster(detectCores() %/% 2)
    
    # Evaluate reactive objects and pass values to new R sessions. Passing
    # large objects to parLapply does not improve performance (to my
    # surprise). Can't index in parLapply for some reason
    hst.lst.prll = sim.lst()$hists.lst
    siips.lst.prll = fsisyips.lst()$siips.lst
    L.prll = L()
    k.prll = k()
    clusterExport(
      cl.tmp,
      list(
        "hst.lst.prll", "siips.lst.prll", "L.prll", "k.prll", 
        "FindGLPs"
      ), 
      environment()
    )
    
    # Update cluster reactive value
    cl(cl.tmp)
  } 
  
  # If fitting offset model for first time since simulation/loading
  if (
    is.null(osisyips.lst()) &&
    any(c("Offset true kinship", "Offset genopair") %in% mdl.st())
  ) {
    osyips = osiips = vector("list", n.sims())
    
    # Loop over histories
    withProgress(
      {
        for (hst.ind in 1:n.sims()) {
          # Sample-year indices
          syis = sisyis.lst()$syis[[hst.ind]]
          
          # Sample index pairs for just consecutive pairs
          osips = FindSIPsOffset(k(), syis)
          
          # Sample-individual and sample-year index pairs, n_pairs x 2,
          # representing individual and survey-year of each sample in each
          # offset pair, counting from zero for survey-years as passing into TMB
          # C++ objective function
          osyips[[hst.ind]] = FindSISYIPs(syis, osips)
          osiips[[hst.ind]] = FindSISYIPs(sisyis.lst()$siis[[hst.ind]], osips)
          
          incProgress(1/n.sims())
        }
      }, value = 0, 
      message = "Finding offset sample-individual and sample-year index pairs"
    )
    
    osisyips.lst(list(osiips = osiips, osyips = osyips))
  }
  
  # Boolean for models requested 
  mod.bool = mdl.chcs %in% mdl.st()

  # Combine and keep only for selected models
  fit.lst.tmp = list(fit.ppn(), fit.ftk(), fit.otk(), fit.fg(), fit.og())
  ests = lapply(fit.lst.tmp, function(fit) fit$ests)[mod.bool]
  ses = lapply(fit.lst.tmp, function(fit) fit$ses)[mod.bool]
  cnvgs = lapply(fit.lst.tmp, function(fit) fit$cnvgs)[mod.bool]
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