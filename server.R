# Load functions
funcs <- list.files("Functions")
for (i in 1:length(funcs)) 
  source(paste0("Functions/", funcs[i]))

# Define server logic for app
server <- function(input, output) {
  # Reactive global variables ----
  
  # Survey years
  srvy.yrs.rct = reactive({
    eval(parse(text = paste0("c(", input$srvy.yrs, ")")))
  })
  # Number of surveys 
  k.rct = reactive(length(srvy.yrs.rct()))
  # Final survey year 
  f.year.rct = reactive(tail(srvy.yrs.rct(), 1))
  # Models to fit
  models = reactive(input$models) 
  
  # Reactive implied parameter outputs ----
  
  # Population growth rate
  output$lambda = renderText({
    paste("Population growth rate (lambda):", input$rho + input$phi)
  })
  # Number of surveys
  output$k = renderText(paste("Number of surveys (k):", k.rct()))
  # Final survey year
  output$f.year = renderText(paste("Final survey year:", f.year.rct()))
  
  # Global variables bound to simulate button ----
  
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
  
  # Functions and outputs bound to simulate button ----
  
  # Simulate population and capture histories
  hists.lst = reactive({
    # Initial population size
    N.init = round(exp.N.t()[1])
    
    # Create list for population and capture histories
    hists.lst <- vector("list", n_sims())

    # Display progress
    cat("Simulating study: ")
    
    # Loop over histories
    withProgress({
      for (hist.ind in 1:n_sims()) {
        # Display progress
        if (hist.ind %% 100 == 1) cat(hist.ind, "")
        
        # Simulate family and capture histories of population of animals over
        # time
        hists.lst[[hist.ind]] <- SimPopStud(
          phi(), lambda(), N.init, hist.len(), srvy.yrs(), k(), f.year()
        )
        
        # Update progress. Unexplained "Error in as.vector: object 'x' not
        # found" seen 19/12/2021 coming from incProgress...
        incProgress(1/n_sims())
      }
    }, value = 0, message = "Simulating populations")
    
    hists.lst
  })
  
  # Check simulated studies
  checks.lst = reactive({
    # Expected final population size
    exp.N.fin <- exp.N.t()[hist.len()]
    
    # Create matrices for population trajectories (as columns), numbers
    # captured, and expected and observed numbers of kin pairs
    N.t.mat <- matrix(nrow = n_sims(), ncol = hist.len())
    ns.caps.mat <- ns.clvng.caps.mat <- ns.clvng.mat <- ns.POPs.wtn.mat <- 
      ns.HSPs.wtn.mat <- exp.ns.HSPs.wtn.mat <- exp.ns.POPs.wtn.mat <- 
      matrix(nrow = n_sims(), ncol = k())
    ns.POPs.btn.mat <- ns.SPs.btn.mat <- exp.ns.POPs.btn.mat <- 
      exp.ns.SPs.btn.mat <- matrix(nrow = n_sims(), ncol = n.srvy.prs())
    
    # Create vectors for proportions with unknown parents, and
    # super-population sizes
    prpn.prnts.unkn.vec <- Ns.vec <- numeric(n_sims())
    
    # Display progress
    cat("Checking study: ")
    
    # Loop over histories
    withProgress({
      for (hist.ind in 1:n_sims()) {
        # Display progress
        if (hist.ind %% 100 == 1) cat(hist.ind, "")
        
        # Get simulated family and capture histories of population of animals
        # over time
        pop.cap.hist <- hists.lst()[[hist.ind]]
        
        # Record population curve
        N.t.mat[hist.ind, ] <- attributes(pop.cap.hist)$N.t.vec
        
        # Record superpopulation size
        Ns.vec[hist.ind] <- attributes(pop.cap.hist)$Ns
        
        # Get numbers captured and calving in each survey
        ns.caps <- attributes(pop.cap.hist)$ns.caps
        ns.caps.mat[hist.ind, ] <- ns.caps
        ns.clvng.mat[hist.ind, ] <- attributes(pop.cap.hist)$ns.clvng
        ns.clvng.caps.mat[hist.ind, ] <- colSums(
          pop.cap.hist[, 4:(3 + k())] * pop.cap.hist[, (4 + k()):(3 + 2 * k())]
        )
        
        # Find proportion captured with unknown parents
        prpn.prnts.unkn.vec[hist.ind] <- mean(is.na(pop.cap.hist$mum))
        
        # Find numbers of known kin pairs
        ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
        
        # Record in matrices
        ns.POPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.wtn
        ns.HSPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.HSPs.wtn
        ns.POPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.btn
        ns.SPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.SPs.btn
        
        # Find expected numbers of kin pairs
        exp.ns.kps.lst <- FindExpNsKPs(
          k(), n.srvy.prs(), exp.N.fin, lambda(), f.year(), srvy.yrs(), phi(), 
          rho(), ns.caps
        )
        
        # Record in matrices
        exp.ns.POPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.wtn
        exp.ns.HSPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.HSPs.wtn
        exp.ns.POPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.btn
        exp.ns.SPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.SPs.btn
        
        incProgress(1/n_sims())
      }
    }, value = 0, message = "Checking simulations")
    
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
      (f.year() - hist.len() + 1):f.year(), t(checks.lst()$N.t.mat), type = 'l',
      col = rgb(0, 0, 0, alpha = 0.1), lty = 1, 
      xlab = 'Year', ylab = 'Nt', main = "Population sizes over time"
    )
    lines((f.year() - hist.len() + 1):f.year(), exp.N.t(), col = 'red', lwd = 2)
    
    # Surveys
    abline(v = srvy.yrs(), lty = 2)
    
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
    # prpns = checks.lst()$ns.POPs.wtn.mat / 
    #   choose(checks.lst()$ns.caps.mat, 2)
    # exp.prpns = 
    #   checks.lst()$exp.ns.POPs.wtn.mat / choose(exp.N.t()[srvy.inds] * p, 2)
    diffs = checks.lst()$ns.POPs.wtn.mat - checks.lst()$exp.ns.POPs.wtn.mat
    boxplot(
      diffs,
      main = 
        "Expected vs observed numbers of parent-offspring pairs within samples",
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
  
  # Display percentage of animals captured for which the parents are unknown
  output$percUnknPrnts <- renderText({
    # nPOPs = checks.lst()$ns.POPs.wtn.mat
    # 
    # # Pairs are lost quadratically with animals
    # prpns.pairs.lost = 1 - (1 - checks.lst()$prpn.prnts.unkn.vec)^2
    # 
    # paste(
    #   "Expected number of parent-offspring pairs within samples
    #   lost due to unknown parents:",
    #   signif(mean(prpns.pairs.lost * checks.lst()$ns.POPs.wtn.mat), 3),
    #   "\n"
    # )
    paste0(
      "Percentage of captured animals with unknown parents (1DP): ", 
      round(mean(checks.lst()$prpn.prnts.unkn.vec) * 100, 1), "%"
    )
  })
  
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
          ppn.tmb.ests[hist.ind, ] <- ppn.tmb.res[["est.se.df"]][, 1]
          ppn.tmb.ses[hist.ind, ] <- ppn.tmb.res[["est.se.df"]][, 2]
          ppn.tmb.cnvg[hist.ind] = ppn.tmb.res[["cnvg"]]
        }
        if (mod.bool[2]) {
          ck.tmb.res <- TryCloseKinTMB(
            k(), srvy.gaps(), f.year(), srvy.yrs(), ns.caps, ns.kps.lst, 
            ck.start, ck.lwr, ck.upr
          )
          ck.tmb.ests[hist.ind, -(5:(4 + k()))] <- 
            ck.tmb.res[["est.se.df"]][, 1]
          ck.tmb.ses[hist.ind, -(5:(4 + k()))] <- ck.tmb.res[["est.se.df"]][, 2]
          ck.tmb.cnvg[hist.ind] = ck.tmb.res[["cnvg"]]
        }
        
        incProgress(1/n_sims())
      }
    }, value = 0, message = "Fitting models")

    # Combine model estimates, standard errors, and convergences, as lists and
    # return those requested
    list(
      ests = list(popan = ppn.tmb.ests, close_kin = ck.tmb.ests)[mod.bool],
      ses = list(popan = ppn.tmb.ses, close_kin = ck.tmb.ses)[mod.bool],
      cnvgs = list(popan = ppn.tmb.cnvg, close_kin = ck.tmb.cnvg)[mod.bool]
    )
  })
  ests.lst = reactive(fit.lst()[["ests"]])
  ses.lst = reactive(fit.lst()[["ses"]])
  cnvgs.lst = reactive(fit.lst()[["cnvgs"]])
  
  # Report number of datasets for which the optimizer reported convergence
  output$cnvgRate = renderText({
    cnvgd = cnvgs.lst()[["close_kin"]] == 0
    paste0(
      "Optimizer reported convergence for ", round(mean(cnvgd) * 100, 1), 
      "% of datasets (", sum(cnvgd), "/", length(cnvgd), 
      ") for close kin model."
    )
  })
  
  # Report number of datasets for which standard errors were calculable
  output$sesRate = renderText({
    ses.clbl = rowSums(is.na(fit.lst()$ses$close_kin[, 1:4])) == 0
    paste0(
      "Standard errors calculable for ", round(mean(ses.clbl) * 100, 1), 
      "% of datasets (", sum(ses.clbl), "/", length(ses.clbl), 
      ") for close kin model."
    )
  })
  
  # Find estimates when optimizer converged.  A list of lists of matrices with
  # one for each model requested
  cvgd.fit.lst = reactive({
    cvgd.ests.lst = cvgd.ses.lst = list()
    # Loop over models requested
    for (i in 1:length(ests.lst())) {
      cvgd.ests.lst[i] = list(ests.lst()[[i]][!cnvgs.lst()[[i]], ])
      cvgd.ses.lst[i] = list(ses.lst()[[i]][!cnvgs.lst()[[i]], ])
    }
    names(cvgd.ests.lst) = names(cvgd.ses.lst) = names(ests.lst())
    list(cvgd.ests.lst = cvgd.ests.lst, cvgd.ses.lst = cvgd.ses.lst)
  })
  cvgd.ests.lst = reactive(cvgd.fit.lst()$cvgd.ests.lst)
  
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
    params = ests.lst()[["close_kin"]][1, 1:3]
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
    abline(v = ests.lst()[["close_kin"]][1, 1], col = 4)
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
    ests.cnvg = cbind(round(ests.lst()[[1]], 3), cnvg = cnvgs.lst()[[1]])
    head(ests.cnvg)
  })
  
  # Plot estimates using model comparison plot function
  output$modComp <- renderPlot({
    # Plot estimates from all models side-by-side
    
    # Set four plots per page
    par(mfrow = c(2, 2))
    
    # Plot estimates for lambda
    ComparisonPlot(
      lapply(cvgd.ests.lst(), function(ests.mat) ests.mat[, 1]), 
      "Population growth rate", lambda()
    )
    
    # Plot estimates for Phi
    ComparisonPlot(
      lapply(cvgd.ests.lst(), function(ests.mat) ests.mat[, 2]),
      "Survival rate", phi()
    )
    
    # Plot estimates of final population size
    ComparisonPlot(
      lapply(cvgd.ests.lst(), function(ests.mat) ests.mat[, 3]),
      "Final population size", exp.N.t()[hist.len()]
    )
    
    # Plot estimates of superpopulation size
    ComparisonPlot(
      lapply(cvgd.ests.lst(), function(ests.mat) ests.mat[, 4]),
      "Super-population size", exp.Ns
    )
  })
  
  # Find confidence intervals for close kin model
  lcbs = reactive(ests.lst()[["close_kin"]] - 1.96 * ses.lst()[["close_kin"]])
  ucbs = reactive(ests.lst()[["close_kin"]] + 1.96 * ses.lst()[["close_kin"]])
  
  # # Check CI doesn't cross parameter bounds
  # ci_ok = reactive(cis()[1, ] > lb() & cis()[2, ] < ub())
  
  # Check CI coverage
  ci_cov = reactive({
    # ci_ok() & true_val() > cis()[1, ] & true_val() < cis()[2, ]
    lambda() > lcbs()[, 1] & lambda() < ucbs()[, 1]
  })
  
  # Plot confidence intervals for close kin model for lambda
  output$CIPlot = renderPlot({
    ord = order(ests.lst()[["close_kin"]][, 1])
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
}