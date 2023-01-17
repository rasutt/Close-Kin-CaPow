# Number of parameters
n.pars = reactive(length(est.par.names()))
# Number of models requested
n.mods = reactive(length(mdl.st()))
# Kinship set Boolean vector
knshp.st.bln = reactive(as.integer(knshp.chcs %in% knshp.st()))

# Allele frequencies for each study, 2 x n_loci X n_sims
ale.frqs.ary = reactive({
  # Make empty list
  ale.frqs.ary.tmp = array(dim = c(2, L(), n.sims()))
  
  # Loop over histories
  withProgress({
    for (hst.ind in 1:n.sims()) {
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
      
      # Allele frequencies, 2 x n_loci matrices, representing relative
      # frequencies of 0 and 1-coded SNP alleles at each locus
      ale.frqs.ary.tmp[, , hst.ind] = 
        FindAleFrqs(attributes(pop.cap.hist)$unq.smp.gts)
      
    }
  }, value = 0, message = "Finding allele frequencies")
  
  ale.frqs.ary.tmp
})

# Possible genotype probabilities for each study, 3 x n_loci x n_sims
pss.gt.prbs.ary = reactive({
  FindPssGtPrbsAry(ale.frqs.ary())
})

# Possible first genotype probabilities for sample pairs for each study, 3 x 3 x
# n_loci x n_sims, first genotypes are rows
pss.gt.1.prbs.ary = reactive({
  FindPssFrstGtPrbsAry(pss.gt.prbs.ary(), L(), n.sims())
})

# Possible genopair probabilities for each kinship for each study, 3 x 3 x
# n_loci x n_sims
pss.gp.prbs.UPs.ary = reactive({
  FindPssGpPsUPsAry(pss.gt.1.prbs.ary())
})
pss.gp.prbs.SPs.ary = reactive({
  FindPssGpPsSPsAry(pss.gt.prbs.ary(), L(), n.sims())
})
pss.gp.prbs.POPs.ary = reactive({
  FindPssGpPsPOPsAry(pss.gt.1.prbs.ary(), ale.frqs.ary(), L(), n.sims())
})
pss.gp.prbs.HSPs.ary = reactive({
  FindPssGpPsHSPsAry(pss.gp.prbs.UPs.ary(), pss.gp.prbs.POPs.ary())
})

# Find full set of genopair log-probabilities for each kinship

# Normal function, taking an array of possible genopair probabilities given a
# certain kinship, over all studies, as an input, and using reactive objects
FindFllLgGpPrbsKP = function(hst.ind, pss.gp.prbs.KP.ary, kp.tp) {
  # Get simulated family and capture histories of population of animals
  # over time
  pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
  
  # Sample history matrix, n_individuals x n_surveys, rows ordered by
  # individual ID
  smp.hsts = as.matrix(pop.cap.hist[, 4:(3 + k())])
  
  # Sample index pairs, 2 x n_pairs, representing indices of samples in each
  # pair
  smp.ind.prs = combn(sum(smp.hsts), 2)
  
  # Sample genotypes, extracted from matrix of individual genotypes,
  # n_samples x n_loci, rows ordered by survey-year then individual ID
  smp.gts = FindSmpGts(smp.hsts, attributes(pop.cap.hist)$unq.smp.gts)
  
  # Possible genopair probabilities over genopairs and loci for this study
  pss.gp.prbs.KP = pss.gp.prbs.KP.ary[, , , hst.ind]
  
  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  FindLogGPProbsKP(
    pss.gp.prbs.KP, smp.gts, smp.ind.prs, L(), sngl.knshp = T, kp.tp
  )
}

# Function using "parallel" library to use multiple CPU cores, giving > 2x
# speedup for datasets >= 100 loci and 100 studies. Sample histories, genotypes,
# and possible genopair probabilities are exported to the new R sessions in
# advance
FindFllLgGpPrbsKPPrll = function(hst.ind, kp.tp) {
  # Get simulated family and capture histories of population of animals
  # over time
  pop.cap.hist <- hst.lst.prll[[hst.ind]]
  
  # Sample history matrix, n_individuals x n_surveys, rows ordered by
  # individual ID
  smp.hsts = as.matrix(pop.cap.hist[, 4:(3 + k.prll)])
  
  # Sample index pairs, 2 x n_pairs, representing indices of samples in each
  # pair
  smp.ind.prs = combn(sum(attributes(pop.cap.hist)$ns.caps), 2)
  
  # Sample genotypes, extracted from matrix of individual genotypes,
  # n_samples x n_loci, rows ordered by survey-year then individual ID
  smp.gts = FindSmpGts(smp.hsts, attributes(pop.cap.hist)$unq.smp.gts)
  
  # Possible genopair probabilities over genopairs and loci for this study
  pss.gp.prbs.KP = pss.gp.prbs.KP.ary[, , , hst.ind]

  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  FindLogGPProbsKP(
    pss.gp.prbs.KP, smp.gts, smp.ind.prs, L.prll, sngl.knshp = T, kp.tp
  )
}

# Set up multicore processing and find log-probabilities given unrelated pairs 
MakeFLGPsKPRctv = function(kp.tp) {
  reactive({
    # Record time
    s = Sys.time()
    
    pss.gp.prbs.KP.ary = switch(
      kp.tp,
      "unrelated" = pss.gp.prbs.UPs.ary(),
      "self" = pss.gp.prbs.SPs.ary(),
      "parent-offspring" = pss.gp.prbs.POPs.ary(),
      "half-sibling" = pss.gp.prbs.HSPs.ary()
    )
    clusterExport(cl(), list("pss.gp.prbs.KP.ary"), environment())
      
    # Find log-probabilities for different chunks of studies on different nodes
    # and combine results when done
    withProgress(
      {
        lst = parLapply(cl(), 1:n.sims(), function(hst.ind) {
          FindFllLgGpPrbsKPPrll(hst.ind, kp.tp)
        })
      },
      value = 0,
      message = paste(
        "Finding genopair log-probabilities given", kp.tp, "pairs"
      ),
      detail = "Using multiple cores so no progress updates available"
    )
    
    # Show time taken
    cat(
      "Found log probabilities, given", kp.tp, "pairs, over",
      L(), "loci and",
      n.sims(), "studies in",
      Sys.time() - s, "\n\n"
    )
    
    # Return results
    lst
  })
}
fll.gp.lg.prbs.UP.lst = MakeFLGPsKPRctv("unrelated")
fll.gp.lg.prbs.SP.lst = MakeFLGPsKPRctv("self")
fll.gp.lg.prbs.POP.lst = MakeFLGPsKPRctv("parent-offspring")
fll.gp.lg.prbs.HSP.lst = MakeFLGPsKPRctv("half-sibling")

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
        alpha = alpha(), knshp_st_bool = knshp.st.bln(),
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
  
  # Start new R sessions on separate "logical" CPU cores (nodes) for finding
  # genopair log-probabilities. Using 6 of my 12 cores seems to be optimal
  cl(makeCluster(6))
  
  # Evaluate reactive objects and pass values to new R sessions. Passing large
  # objects to parLapply does not improve performance (to my surprise). Can't
  # index in parLapply for some reason
  hst.lst.prll = sim.lst()$hists.lst
  L.prll = L()
  k.prll = k()
  clusterExport(
    cl(),
    list("hst.lst.prll", "L.prll", "k.prll", "FindSmpGts", "FindLogGPProbsKP"), 
    environment()
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
      
      # Find genopair log-probabilities given selected kinships
      fll.gp.lg.prbs.KP = cbind(
        fll.gp.lg.prbs.UP.lst()[[hst.ind]],
        if ("Half-sibling" %in% knshp.st()) 
          fll.gp.lg.prbs.HSP.lst()[[hst.ind]],
        if ("Parent-offspring" %in% knshp.st()) 
          fll.gp.lg.prbs.POP.lst()[[hst.ind]],
        if ("Self" %in% knshp.st()) 
          fll.gp.lg.prbs.SP.lst()[[hst.ind]]
      )
      
      # Exponentiate genopair log-probabilities given kinship set, checking for
      # underflow and trying to adjust if necessary.  Would be good to raise a
      # proper error if adjustment impossible
      fll.gp.prbs.KP = FindGPPs(fll.gp.lg.prbs.KP)

      print("fll.gp.prbs.KP")
      print(str(fll.gp.prbs.KP))
      print(head(fll.gp.prbs.KP))
      # print("Correct, n_pairs x n_kinships")
      # 
      # print("fll.SYIPs.lst()[[hst.ind]]")
      # print(str(fll.SYIPs.lst()[[hst.ind]]))
      # print("Correct, n_pairs x 2")
      # 
      # print("knshp.st.bln")
      # print(knshp.st.bln())
      
      # Make objective function
      obj = MakeTMBObj(
        ck.start, "genopair",
        k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
        alpha = alpha(), knshp_st_bool = knshp.st.bln(),
        gp_probs = fll.gp.prbs.KP, 
        smp_yr_ind_prs = fll.SYIPs.lst()[[hst.ind]]
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
  
  # Stop R sessions on other nodes
  stopCluster(cl())
  
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
        alpha = alpha(), knshp_st_bool = knshp.st.bln(),
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

# Nullify objects when new datasets simulated
observeEvent(input$simulate, {
  knshp.st(NULL)
  fit.lst(NULL)
  SYIs.lst(NULL)
  fll.SYIPs.lst(NULL)
  offst.SYIPs.lst(NULL)
})

# When fit models button clicked
observeEvent(input$fit, {
  # Update reactive values for model and kinship sets
  mdl.st(input$mdl.st)
  knshp.st(input$knshp.st)
  
  # If fitting a genopair model for first time find sample-year indices
  if (
    is.null(SYIs.lst()) &&
    any(c("Full genopair", "Offset genopair") %in% mdl.st())
  ) {
    # Make empty list
    SYIs.lst.tmp = vector("list", n.sims())
    
    # Loop over histories
    withProgress({
      for (hst.ind in 1:n.sims()) {
        # Get simulated family and capture histories of population of animals
        # over time
        pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
        
        # Sample history matrix, n_individuals x n_surveys, rows ordered by
        # individual ID
        smp.hsts = as.matrix(pop.cap.hist[, 4:(3 + k())])
        
        # Sample-year indices, columns in sample history matrix, representing
        # the survey that each sample came from, ordered by survey-year then
        # individual ID.  Counting from zero as will be passed to TMB objective
        # function
        SYIs.lst.tmp[[hst.ind]] = col(smp.hsts)[as.logical(smp.hsts)] - 1
        
        incProgress(1/n.sims())
      }
    }, value = 0, message = "Finding sample-year indices")
    
    SYIs.lst(SYIs.lst.tmp)
  }
  
  # If fitting full genopair model for first time find survey-year index pairs
  if (
    is.null(fll.SYIPs.lst()) &&
    "Full genopair" %in% mdl.st()
  ) {
    fll.SYIPs.lst.tmp = vector("list", n.sims())
    
    # Loop over histories
    withProgress({
      for (hst.ind in 1:n.sims()) {
        # Sample-year indices, columns in sample history matrix, representing
        # the survey that each sample came from, ordered by survey-year then
        # individual ID.  Counting from zero as will be passed to TMB objective
        # function
        smp.yr.inds = SYIs.lst()[[hst.ind]]
        
        # Sample index pairs, 2 x n_pairs, representing indices of samples in
        # all pairs
        smp.ind.prs.fll = combn(length(smp.yr.inds), 2)

        # Sample-year index pairs, n_pairs x 2, representing survey-year of each
        # sample in each pair, counting from zero as passing into TMB C++
        # objective function
        fll.SYIPs.lst.tmp[[hst.ind]] = 
          matrix(smp.yr.inds[as.vector(t(smp.ind.prs.fll))], ncol = 2)

        incProgress(1/n.sims())
      }
    }, value = 0, message = "Finding all sample-year index pairs")
    
    fll.SYIPs.lst(fll.SYIPs.lst.tmp)
  }
  
  # If fitting offset genopair model for first time find offset survey-year
  # index pairs
  if (
    is.null(offst.SYIPs.lst()) &&
    "Offset genopair" %in% mdl.st()
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
        
        # Sample index pairs for just consecutive
        # pairs
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