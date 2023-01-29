# Kinship indicator variables, binary variable for each kinship offered
kivs = reactive(as.integer(knshp.chcs %in% knshp.st()))

# Fit popan model
fit.ppn = reactive(if ("Popan" %in% mdl.st()) {
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
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
      
      # Get numbers of animals captured in study
      n.cap.hists <- nrow(pop.cap.hist)
      
      # Update optimiser starting-values and bounds
      ppn.start[3] <- attributes(pop.cap.hist)$Ns
      ppn.lwr[3] <- n.cap.hists
      
      # Summarise data for POPAN model
      pop.sum <- FindPopSum(k(), pop.cap.hist, n.cap.hists)
      
      # Create TMB function
      ppn.obj = MakeTMBObj(
        ppn.start, "popan",
        k(), srvy.gaps(), 
        n_cap_hists = n.cap.hists, first_tab = pop.sum$first.tab, 
        last_tab = pop.sum$last.tab, caps = pop.sum$caps, 
        non_caps = pop.sum$non.caps, survives = pop.sum$survives[-k()]
      )
      
      # Try to fit models
      ppn.tmb.res = TryModelTMB(ppn.obj, ppn.lwr, ppn.upr, "popan")
      ppn.tmb.ests[hst.ind, ] <- ppn.tmb.res$est.se.df[, 1]
      ppn.tmb.ses[hst.ind, ] <- ppn.tmb.res$est.se.df[, 2]
      ppn.tmb.cnvg[hst.ind] = ppn.tmb.res$cnvg
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting Popan model")
  
  # Combine model estimates, standard errors, and convergences, as lists and
  # return those requested
  list(ests = ppn.tmb.ests, ses = ppn.tmb.ses, cnvgs = !ppn.tmb.cnvg)
})

# Fit close-kin model
fit.ck = reactive(if ("True kinship" %in% mdl.st()) {
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
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
      
      # Get numbers of animals captured in each survey
      ns.caps <- attributes(pop.cap.hist)$ns.caps
      
      # Find numbers of kin pairs
      ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
      
      # Update optimiser starting-values and bounds
      ck.start[3] <- attributes(pop.cap.hist)$N.t.vec[hist.len()]
      ck.lwr[3] <- ns.caps[k()]
      
      # Create TMB function
      obj = MakeTMBObj(
        ck.start, "true kinship",
        k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
        alpha = alpha(), knshp_st_bool = kivs(),
        ns_SPs_btn = ns.kps.lst$btn[1, ], ns_POPs_wtn = ns.kps.lst$wtn[1, ], 
        ns_POPs_btn = ns.kps.lst$btn[2, ], ns_HSPs_wtn = ns.kps.lst$wtn[2, ],
        ns_HSPs_btn = ns.kps.lst$btn[3, ], ns_caps = ns.caps
      )
      
      # Try to fit close-kin likelihood model
      ck.tmb.res = TryModelTMB(obj, ck.lwr, ck.upr, "true kinship")
      
      # If optimiser did not give error
      if(!all(is.na(ck.tmb.res))) {
        # Store results separately
        ck.tmb.ests[hst.ind, -(5:(4 + k()))] <- ck.tmb.res$est.se.df[, 1]
        ck.tmb.ses[hst.ind, -(5:(4 + k()))] <- ck.tmb.res$est.se.df[, 2]
        ck.tmb.cnvg[hst.ind] = ck.tmb.res$cnvg
      }
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting close-kin model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = ck.tmb.ests, ses = ck.tmb.ses, cnvgs = !ck.tmb.cnvg)
})

# Fit genopair model
fit.gp = reactive(if ("Full genopair" %in% mdl.st()) {
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
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
      
      # Update optimizer starting-values and bounds
      ck.start[3] <- attributes(pop.cap.hist)$N.t.vec[hist.len()]
      ck.lwr[3] <- attributes(pop.cap.hist)$ns.caps[k()]
      
      # Full genopair log-probabilities given selected kinships, n_pairs x
      # n_kinships, combining those for unrelated pairs and selected kinships
      fglps.gvn.ks = cbind(
        fglps.UPs.lst()[[hst.ind]],
        if ("Half-sibling" %in% knshp.st()) fglps.HSPs.lst()[[hst.ind]],
        if ("Parent-offspring" %in% knshp.st()) fglps.POPs.lst()[[hst.ind]],
        if ("Self" %in% knshp.st()) fglps.SPs.lst()[[hst.ind]]
      )
      colnames(fglps.gvn.ks) = 
        c("UP", "HSP", "POP", "SP")[c(T, rev(kivs()) == 1)]
      
      # Full genopair probabilities given kinships, n_pairs x n_kinships, 
      # exponentiating log-probabilities, checking for underflow and trying to
      # adjust if necessary.  Would be good to raise a proper error if
      # adjustment impossible
      fgps.gvn.ks = FindGPsGvnKs(fglps.gvn.ks)
      cat("Full genopair probabilities given kinships \n")
      print(head(fgps.gvn.ks))
      cat("\n")

      # Survey-year index pairs, n_pairs x 2, starting from zero for C++
      # objective function
      cat("Survey-year index pairs \n")
      syips = data.frame(fsisyips.lst()$syips.lst[[hst.ind]])
      names(syips) = paste("Sample", 1:2)
      print(table((syips)))
      cat("\n")
      
      # print("fsisyips.lst()$siips.lst[[hst.ind]]")
      # print(str(fsisyips.lst()$siips.lst[[hst.ind]]))
      # 
      # print("pss.gt.prbs.ary()[, 1:2, hst.ind]")
      # print(pss.gt.prbs.ary()[, 1:2, hst.ind])
      # 
      # print("pss.gt.1.prbs.ary()[, , 1:2, hst.ind]")
      # print(pss.gt.1.prbs.ary()[, , 1:2, hst.ind])
      # 
      # print("pss.gp.prbs.UPs.ary()[, , 1:2, hst.ind]")
      # print(pss.gp.prbs.UPs.ary()[, , 1:2, hst.ind])
      # 
      # print("pss.gp.prbs.HSPs.ary()[, , 1:2, hst.ind]")
      # print(pss.gp.prbs.HSPs.ary()[, , 1:2, hst.ind])
      # 
      # print("pss.gp.prbs.POPs.ary()[, , 1:2, hst.ind]")
      # print(pss.gp.prbs.POPs.ary()[, , 1:2, hst.ind])
      # 
      # print("kivs")
      # print(kivs())
      
      # Make objective function
      obj = MakeTMBObj(
        ck.start, "genopair",
        k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
        alpha = alpha(), knshp_st_bool = kivs(),
        gp_probs = fgps.gvn.ks, 
        smp_yr_ind_prs = fsisyips.lst()$syips.lst[[hst.ind]]
      )
      
      # Try to fit genopair likelihood model
      gp.tmb.res = TryModelTMB(obj, ck.lwr, ck.upr, "genopair")
      
      # If optimiser did not give error
      if(!all(is.na(gp.tmb.res))) {
        gp.tmb.ests[hst.ind, -(5:(4 + k()))] <- gp.tmb.res$est.se.df[, 1]
        gp.tmb.ses[hst.ind, -(5:(4 + k()))] <- gp.tmb.res$est.se.df[, 2]
        gp.tmb.cnvg[hst.ind] = gp.tmb.res$cnvg
      }
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting genopair model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = gp.tmb.ests, ses = gp.tmb.ses, cnvgs = !gp.tmb.cnvg)
})

# Fit offset model
fit.os = reactive(if ("Offset genopair" %in% mdl.st()) {
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), NA)
  ck.lwr <- c(0, 0.75, NA)
  ck.upr <- c(0.35, 1, Inf)
  
  # Create vector for model convergences
  os.tmb.cnvg <- numeric(n.sims())
  
  # Create matrices for estimates and standard errors
  os.tmb.ests <- os.tmb.ses <- matrix(
    nrow = n.sims(), ncol = 4 + k(), dimnames = list(NULL, est.par.names())
  )
  
  # Loop over histories
  withProgress({
    for (hst.ind in 1:n.sims()) {
      # Display progress
      cat("History:", hst.ind, "\n")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
      
      # Update optimiser starting-values and bounds
      ck.start[3] <- attributes(pop.cap.hist)$N.t.vec[hist.len()]
      ck.lwr[3] <- attributes(pop.cap.hist)$ns.caps[k()]
      
      # Genopair model inputs, list of 2, genopair probabilities given kinship
      # set, n_pairs x n_kinships, and sample-year index pairs, n_pairs x 2
      Gp.Mdl.Inpts = 
        FindGpMdlInpts(pop.cap.hist, L(), k(), os.mdl = T, knshp.st())
      
      # Make objective function
      obj = MakeTMBObj(
        ck.start, "genopair",
        k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
        alpha = alpha(), knshp_st_bool = kivs(),
        gp_probs = Gp.Mdl.Inpts$GPPs, 
        smp_yr_ind_prs = Gp.Mdl.Inpts$smp.yr.ind.prs
      )
      
      # Try to fit genopair likelihood model
      os.tmb.res = TryModelTMB(obj, ck.lwr, ck.upr, "genopair")
      
      # If optimiser did not give error
      if(!all(is.na(os.tmb.res))) {
        os.tmb.ests[hst.ind, -(5:(4 + k()))] <- os.tmb.res$est.se.df[, 1]
        os.tmb.ses[hst.ind, -(5:(4 + k()))] <- os.tmb.res$est.se.df[, 2]
        os.tmb.cnvg[hst.ind] = os.tmb.res$cnvg
      }
      
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Fitting offset model")
  
  # Combine model estimates, standard errors, and convergences, and return as
  # lists
  list(ests = os.tmb.ests, ses = os.tmb.ses, cnvgs = !os.tmb.cnvg)
})

