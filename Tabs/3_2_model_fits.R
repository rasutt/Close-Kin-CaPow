# Reactives to fit models when requested

# Kinship indicator variables, binary variable for each kinship offered
kivs = reactive(as.integer(knshp.chcs %in% knshp.st()))

# Fit popan model
fit.ppn = reactive(if ("Popan" %in% mdl.st()) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start = c(rho(), phi(), NA)
  ck.lwr = c(0, 0.75, NA)
  ck.upr = c(0.35, 1, Inf)
  ppn.start = cbd.start = c(ck.start, rep(p(), k()))
  ppn.lwr = cbd.lwr = c(ck.lwr, rep(0, k()))
  ppn.upr = cbd.upr = c(ck.upr, rep(1, k()))
  
  # Create vector model convergences
  ppn.cnvg = numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  ppn.ests = ppn.ses = matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      stdy = sim.lst()$hists.lst[[hst.ind]]
      
      # Get numbers of animals captured in study
      n.cap.hists = nrow(stdy)
      
      # Update optimizer starting-values and bounds
      ppn.start[3] = attributes(stdy)$Ns
      ppn.lwr[3] = n.cap.hists
      
      # Summarise data for POPAN model
      pop.sum = FindPopSum(k(), stdy, n.cap.hists)
      
      # Create TMB function
      ppn.obj = MakeTMBObj(
        ppn.start, "popan",
        k(), srvy.gaps(), 
        n_cap_hists = n.cap.hists, first_tab = pop.sum$first.tab, 
        last_tab = pop.sum$last.tab, caps = pop.sum$caps, 
        non_caps = pop.sum$non.caps, survives = pop.sum$survives[-k()]
      )
      
      # Try to fit models
      ppn.res = TryModelTMB(ppn.obj, ppn.lwr, ppn.upr, "popan")
      ppn.ests[hst.ind, ] = ppn.res$est.se.df[, 1]
      ppn.ses[hst.ind, ] = ppn.res$est.se.df[, 2]
      ppn.cnvg[hst.ind] = ppn.res$cnvg
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting Popan model")
  
  # Combine model estimates, standard errors, and convergences, as lists and
  # return those requested
  list(ests = ppn.ests, ses = ppn.ses, cnvgs = !ppn.cnvg)
})

# Fit full true kinship model
fit.ftk = reactive(if ("Full true kinship" %in% mdl.st()) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start = c(rho(), phi(), NA)
  ck.lwr = c(0, 0.75, NA)
  ck.upr = c(0.35, 1, Inf)
  
  # Create vector for model convergences
  ck.cnvg = numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  ck.ests = ck.ses = matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      stdy = sim.lst()$hists.lst[[hst.ind]]
      
      # Get numbers of animals captured in each survey
      ns.caps = attributes(stdy)$ns.caps
      
      # Find numbers of kin pairs
      ns.kps.lst = FindNsKinPairs(k(), n.srvy.prs(), stdy)
      
      # Update optimiser starting-values and bounds
      ck.start[3] = attributes(stdy)$N.t.vec[hist.len()]
      ck.lwr[3] = ns.caps[k()]
      
      # Create TMB function
      obj = MakeTMBObj(
        ck.start, "full true kinship",
        k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
        alpha = alpha(), knshp_st_bool = kivs(),
        ns_SPs_btn = ns.kps.lst$btn[1, ], ns_POPs_wtn = ns.kps.lst$wtn[1, ], 
        ns_POPs_btn = ns.kps.lst$btn[2, ], ns_HSPs_wtn = ns.kps.lst$wtn[2, ],
        ns_HSPs_btn = ns.kps.lst$btn[3, ], ns_caps = ns.caps
      )
      
      # Try to fit close-kin likelihood model
      ck.res = TryModelTMB(obj, ck.lwr, ck.upr, "full true kinship")
      
      # If optimiser did not give error
      if(!all(is.na(ck.res))) {
        # Store results separately
        ck.ests[hst.ind, -(5:(4 + k()))] = ck.res$est.se.df[, 1]
        ck.ses[hst.ind, -(5:(4 + k()))] = ck.res$est.se.df[, 2]
        ck.cnvg[hst.ind] = ck.res$cnvg
      }
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting full true kinship model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = ck.ests, ses = ck.ses, cnvgs = !ck.cnvg)
})

# Fit offset true kinships model
fit.otk = reactive(if ("Offset true kinship" %in% mdl.st()) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start = c(rho(), phi(), NA)
  ck.lwr = c(0, 0.75, NA)
  ck.upr = c(0.35, 1, Inf)
  
  # Create vector for model convergences
  ck.cnvg = numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  ck.ests = ck.ses = matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      stdy = sim.lst()$hists.lst[[hst.ind]]
      
      # Offset sample-year index pairs
      osyips = osisyips()$osyips[[hst.ind]]
      
      # Find numbers of kin pairs
      ns.kps.lst = FindNsOKPs(
        k(), n.srvy.prs(), stdy, osisyips()$osiips[[hst.ind]], osyips
      )
      
      # Get numbers of pairs for each pair of surveys
      ns.pairs = table(data.frame(osyips))
      
      # Update optimiser starting-values and bounds
      ck.start[3] = attributes(stdy)$N.t.vec[hist.len()]
      ck.lwr[3] = attributes(stdy)$ns.caps[k()]
      
      # Create TMB function
      obj = MakeTMBObj(
        ck.start(), "offset true kinship",
        k(), srvy.gaps(), fnl.year(), srvy.yrs(),
        alpha = alpha(), knshp_st_bool = kivs(),
        ns_SPs_btn = ns.kps.lst$btn[1, ], ns_POPs_wtn = ns.kps.lst$wtn[1, ],
        ns_POPs_btn = ns.kps.lst$btn[2, ], ns_HSPs_wtn = ns.kps.lst$wtn[2, ],
        ns_HSPs_btn = ns.kps.lst$btn[3, ], 
        ns_pairs_wtn = diag(ns.pairs), 
        ns_pairs_btn = t(ns.pairs)[lower.tri(ns.pairs)]
      )

      # Try to fit close-kin likelihood model
      ck.res = TryModelTMB(obj, ck.lwr, ck.upr, "offset true kinship")
      
      # If optimizer did not give error
      if(!all(is.na(ck.res))) {
        # Store results separately
        ck.ests[hst.ind, -(5:(4 + k()))] = ck.res$est.se.df[, 1]
        ck.ses[hst.ind, -(5:(4 + k()))] = ck.res$est.se.df[, 2]
        ck.cnvg[hst.ind] = ck.res$cnvg
      }
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting offset true kinship model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = ck.ests, ses = ck.ses, cnvgs = !ck.cnvg)
})

# Fit genopair model
fit.fg = reactive(if ("Full genopair" %in% mdl.st()) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start = c(rho(), phi(), NA)
  ck.lwr = c(0, 0.75, NA)
  ck.upr = c(0.35, 1, Inf)
  
  # Create vector for model convergences
  gp.cnvg = numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  gp.ests = gp.ses = matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      stdy = sim.lst()$hists.lst[[hst.ind]]
      
      # Update optimizer starting-values and bounds
      ck.start[3] = attributes(stdy)$N.t.vec[hist.len()]
      ck.lwr[3] = attributes(stdy)$ns.caps[k()]
      
      # Full genopair log-probabilities given selected kinships, n_pairs x
      # n_kinships, combining those for unrelated pairs and selected kinships
      fglps = cbind(
        fglps.UPs()[[hst.ind]],
        if ("Half-sibling" %in% knshp.st()) fglps.HSPs()[[hst.ind]],
        if ("Parent-offspring" %in% knshp.st()) fglps.POPs()[[hst.ind]],
        if ("Self" %in% knshp.st()) fglps.SPs()[[hst.ind]]
      )
      colnames(fglps) = c("UP", "HSP", "POP", "SP")[c(T, rev(kivs()) == 1)]
      
      # Make objective function
      obj = MakeTMBObj(
        ck.start, "genopair",
        k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
        alpha = alpha(), knshp_st_bool = kivs(),
        gp_probs = FindGPs(fglps), 
        smp_yr_ind_prs = fsisyips()$syips[[hst.ind]]
      )
      
      # Try to fit genopair likelihood model
      gp.res = TryModelTMB(obj, ck.lwr, ck.upr, "full genopair")
      
      # If optimizer did not give error
      if(!all(is.na(gp.res))) {
        gp.ests[hst.ind, -(5:(4 + k()))] = gp.res$est.se.df[, 1]
        gp.ses[hst.ind, -(5:(4 + k()))] = gp.res$est.se.df[, 2]
        gp.cnvg[hst.ind] = gp.res$cnvg
      }
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting full genopair model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = gp.ests, ses = gp.ses, cnvgs = !gp.cnvg)
})

# Fit offset model
fit.og = reactive(if ("Offset genopair" %in% mdl.st()) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start = c(rho(), phi(), NA)
  ck.lwr = c(0, 0.75, NA)
  ck.upr = c(0.35, 1, Inf)
  
  # Create vector for model convergences
  og.cnvg = numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  og.ests = og.ses = matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      stdy = sim.lst()$hists.lst[[hst.ind]]
      
      # Update optimizer starting-values and bounds
      ck.start[3] = attributes(stdy)$N.t.vec[hist.len()]
      ck.lwr[3] = attributes(stdy)$ns.caps[k()]
      
      # Genopair probabilities given each kinship selected, n_pairs x
      # n_kinships
      gps = FindGPsMdl(
        stdy, L(), knshp.st(), osisyips()$osiips[[hst.ind]]
      )
      
      # Make objective function
      obj = MakeTMBObj(
        ck.start, "genopair",
        k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
        alpha = alpha(), knshp_st_bool = kivs(),
        gp_probs = gps, smp_yr_ind_prs = osisyips()$osyips[[hst.ind]]
      )
      
      # Try to fit genopair likelihood model
      og.res = TryModelTMB(obj, ck.lwr, ck.upr, "offset genopair")
      
      # If optimizer did not give error
      if(!all(is.na(og.res))) {
        og.ests[hst.ind, -(5:(4 + k()))] = og.res$est.se.df[, 1]
        og.ses[hst.ind, -(5:(4 + k()))] = og.res$est.se.df[, 2]
        og.cnvg[hst.ind] = og.res$cnvg
      }
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting offset genopair model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = og.ests, ses = og.ses, cnvgs = !og.cnvg)
})

