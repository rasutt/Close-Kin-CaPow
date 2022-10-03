# Objects for checks sub-tabs

# Nullify checks for simulated studies when new studies simulated
observeEvent(input$simulate, {
  checks.lst(NULL)
  N.t.mat(NULL)
  ns.SPs(NULL)
  ns.POPs(NULL)
  ns.SMPs(NULL)
  ns.SFPs(NULL)
  ns.SibPs(NULL)
  prpn.unkn.prnts(NULL)
})

## Reactives returning functions to find numbers of kin-pairs

# In survey-years
find.KPs.wtn = reactive(function(find.func, KP.type) {
  # Need to specify matrix in case just one column or sapply returns a vector
  withProgress({
    matrix(
      sapply(sim.lst()$hists.lst, function(hist) {
        # Increment progress-bar
        incProgress(1/n.sims())
        
        # Find kin-pairs
        find.func(attributes(hist), k())
      }),
      ncol = k(), dimnames = list(NULL, Survey = srvy.yrs()), byrow = T
    )
  }, value = 0, message = paste("Finding", KP.type, "pairs in survey-years"))
})

# Between pairs of survey-years
find.KPs.btn = reactive(function(find.func, KP.type) {
  # Need to specify matrix in case just one column or sapply returns a vector
  withProgress({
    matrix(
      sapply(sim.lst()$hists.lst, function(hist) {
        # Increment progress-bar
        incProgress(1/n.sims())
        
        # Find kin-pairs
        find.func(attributes(hist), k())
      }),
      ncol = n.srvy.prs(), dimnames = list(NULL, Survey = srvy.prs()), 
      byrow = T
    )
  }, 
  value = 0, 
  message = paste("Finding", KP.type, "pairs between survey-years"))
})

# In years up to final year
find.KPs.t = reactive(function(find.func, KP.type) {
  # Need to specify matrix in case just one column or sapply returns a vector
  withProgress({
    matrix(
      sapply(sim.lst()$hists.lst, function(hist) {
        # Increment progress-bar
        incProgress(1/n.sims())
        
        # Find kin-pairs
        find.func(attributes(hist), n.yrs.chk.t())
      }),
      ncol = n.yrs.chk.t(), dimnames = list(NULL, Year = yrs.chk.t()), 
      byrow = T
    )
  }, 
  value = 0, 
  message = paste("Finding", KP.type, "pairs in years up to final year"))
})

# Compute separate checks ----
observeEvent({
  input$check.sub.tabs
  input$nav.tab
}, {
  if (input$nav.tab == "check.tab") {
    # Find population sizes over time
    if (
      input$check.sub.tabs %in% c("populations", "all.pairs") && 
      is.null(N.t.mat())
    ) {
      withProgress({
        N.t.mat(t(sapply(sim.lst()$hists.lst, function(hist) {
          # Increment progress-bar
          incProgress(1/n.sims())
          
          # Record population curve
          attributes(hist)$N.t.vec
        })))
      }, value = 0, message = "Finding population sizes")
    }
    
    # Find numbers of self-pairs between survey-years in simulated populations
    if (
      input$check.sub.tabs == "SPs.tab" && 
      is.null(ns.SPs())
    ) {
      ns.SPs(find.KPs.btn()(find.SPs, "self-pairs"))
    }

    # Find proportions of individuals with unknown parents
    if (
      input$check.sub.tabs %in% 
      c("POPs.tab", "SMPs.tab", "SFPs.tab", "SibPs.tab", "bias.tab") && 
      is.null(prpn.unkn.prnts())
    ) {
      # Update reactive value
      prpn.unkn.prnts(list(
        # In survey-years
        pns.UPs.wtn = find.KPs.wtn()(find.pns.UPs.wtn, "unknown parents"),
        
        # In pairs of survey-years
        pns.UPs.btn = find.KPs.btn()(find.pns.UPs.btn, "unknown parents")
      ))
    }
    
    # Find numbers of parent-offspring pairs
    if (
      input$check.sub.tabs == "POPs.tab" && 
      is.null(ns.POPs())
    ) {
      # Update reactive value
      ns.POPs(list(
        # In survey-years
        ns.POPs.wtn = find.KPs.wtn()(find.POPs.wtn, "parent-offspring"),
        
        # Between pairs of survey-years
        ns.POPs.btn = find.KPs.btn()(find.POPs.btn, "parent-offspring")
      ))
    }
    
    # Find numbers of same-mother pairs
    if (
      input$check.sub.tabs %in% c("SMPs.tab", "SibPs.tab") && 
      is.null(ns.SMPs.wtn())
    ) {
      # In survey-years
      ns.SMPs.wtn(find.KPs.wtn()(find.SMPs.wtn, "same-mother"))
    }
    if (
      input$check.sub.tabs == "SMPs.tab" && 
      is.null(ns.SMPs())
    ) {
      # Update reactive value
      ns.SMPs(list(
        # In survey-years
        ns.SMPs.wtn = ns.SMPs.wtn(),
        
        # Between pairs of survey-years
        ns.SMPs.btn = find.KPs.btn()(find.SMPs.btn, "same-mother"),
        
        # Same-mother pairs with ages known, in final year
        ns.SMPs.age.knwn = find.KPs.t()(find.SMPs.age.knwn, "same-mother")
      ))
    }
    
    # Find numbers of same-father pairs
    if (
      input$check.sub.tabs %in% c("SFPs.tab", "SibPs.tab") && 
      is.null(ns.SFPs.wtn())
    ) {
      # In survey-years
      ns.SFPs.wtn(find.KPs.wtn()(find.SFPs.wtn, "same-father"))
    }
    if (
      input$check.sub.tabs == "SFPs.tab" && 
      is.null(ns.SFPs())
    ) {
      # Update reactive value
      ns.SFPs(list(
        ns.SFPs.wtn = ns.SFPs.wtn(), 
        
        # # Between pairs of survey-years
        # ns.SFPs.btn = find.KPs.btn()(find.SFPs.btn, "same-father"),
        
        # Same-father pairs with ages known, in final year
        ns.SFPs.age.knwn = find.KPs.t()(find.SFPs.age.knwn, "same-father"),
        
        # Same-father pairs with same ages known, in final year
        ns.SFPs.same.age = find.KPs.t()(find.SFPs.same.age, "same-father")
      ))
    }
    
    # Find numbers of full and half-sibling pairs
    if (
      input$check.sub.tabs == "SibPs.tab" && 
      is.null(ns.SibPs())
    ) {
      # Full-sibling pairs within survey-years
      ns.FSPs.wtn = find.KPs.wtn()(find.FSPs.wtn, "full-sibling")
      
      # Update reactive value
      ns.SibPs(list(
        ns.FSPs.wtn = ns.FSPs.wtn,
        
        # Half-sibling pairs within survey years
        ns.HSPs.wtn = ns.SMPs.wtn() + ns.SFPs.wtn() - 2 * ns.FSPs.wtn
      ))
    }
  }
})

# Checks based on population sizes
N.s.yrs = reactive(N.t.mat()[, s.yr.inds()])
ns.APs.wtn.pop = reactive({
  ns = choose(N.s.yrs(), 2)
  colnames(ns) = srvy.yrs()
  ns
})
ns.APs.btn.pop = reactive({
  # If there is only one survey-pair apply returns a vector so we have to make
  # it a matrix explicitly 
  ns = t(matrix(
    apply(N.s.yrs(), 1, combn, 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2]),
    nrow = n.srvy.prs()
  ))
  colnames(ns) = srvy.prs()
  ns
})

# Compute combined checks ----
observeEvent({
  input$check.sub.tabs
  input$nav.tab
}, {
  if (
    input$nav.tab == "check.tab" &&
    !(input$check.sub.tabs %in% 
      c("populations", "all.pairs", "SPs.tab", "POPs.tab", "SMPs.tab", 
        "SFPs.tab", "SibPs.tab")) && 
    is.null(checks.lst())
  ) {
    # Objects to store results
    # ns.caps.mat = ns.clvng.caps.mat = ns.clvng.mat = 
    #   matrix(nrow = n.sims(), ncol = k())
    prpn.prnts.unkn.vec = Ns.vec = numeric(n.sims())
    ns.kps.pop.wtn.arr = array(
      dim = c(n.sims(), k(), n.kp.tps.pop.wtn),
      dimnames = list(NULL, Survey = srvy.yrs(), kp.type = kp.tps.pop.wtn)
    )
    ns.kps.pop.btn.arr = array(
      dim = c(n.sims(), n.srvy.prs(), n.kp.tps.pop.btn),
      dimnames = list(NULL, Survey_pair = srvy.prs(), kp.type = kp.tps.pop.btn)
    )
    # ns.kps.cap.wtn.arr = exp.ns.kps.cap.wtn.arr = array(
    #   dim = c(n.sims(), k(), n.kp.tps.cap.wtn),
    #   dimnames = list(NULL, Survey = srvy.yrs(), kp.type = kp.tps.cap.wtn)
    # )
    # ns.kps.cap.btn.arr = exp.ns.kps.cap.btn.arr = array(
    #   dim = c(n.sims(), n.srvy.prs(), n.kp.tps.cap.btn),
    #   dimnames = list(NULL, Survey_pair = srvy.prs(), kp.type = kp.tps.cap.btn)
    # )
    ns.kps.t.arr = array(
      dim = c(n.sims(), n.yrs.chk.t(), n.kp.tps.t),
      dimnames = list(NULL, Year = yrs.chk.t(), kp.type = kp.tps.t)
    )
    
    # Loop over histories
    withProgress({
      for (hist.ind in 1:n.sims()) {
        # Get simulated family and capture histories of population of animals
        # over time
        pop.cap.hist = sim.lst()$hists.lst[[hist.ind]]
        
        # Insert population sizes in survey years for comparison with numbers of
        # kin pairs
        ns.kps.pop.wtn.arr[hist.ind, , 1] = 
          attributes(pop.cap.hist)$N.t.vec[s.yr.inds()]
        
        # Record super-population size
        Ns.vec[hist.ind] = attributes(pop.cap.hist)$Ns
        
        # Get numbers captured and calving in each survey
        ns.caps = attributes(pop.cap.hist)$ns.caps
        # ns.caps.mat[hist.ind, ] = ns.caps
        # ns.clvng.mat[hist.ind, ] = attributes(pop.cap.hist)$ns.clvng
        # ns.clvng.caps.mat[hist.ind, ] = colSums(
        #   pop.cap.hist[, 4:(3 + k())] * pop.cap.hist[, (4 + k()):(3 + 2 * k())]
        # )
        
        # Find proportion captured with unknown parents
        prpn.prnts.unkn.vec[hist.ind] = mean(is.na(pop.cap.hist$mum))
        
        # Numbers of kin-pairs in whole population
        ns.kps.pop.lst = FindNsKinPairsPop(pop.cap.hist, s.yr.inds(), k())
        ns.kps.pop.wtn.arr[hist.ind, , -1] = t(ns.kps.pop.lst$wtn)
        ns.kps.pop.btn.arr[hist.ind, , ] = t(ns.kps.pop.lst$btn)
        
        # # Find numbers of kin pairs among samples
        # ns.kps.cap.lst = FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
        # ns.kps.cap.wtn.arr[hist.ind, , ] = t(ns.kps.cap.lst$wtn)
        # ns.kps.cap.btn.arr[hist.ind, , ] = t(ns.kps.cap.lst$btn)
        # 
        # # Find expected numbers of kin pairs among samples
        # exp.ns.kps.cap.lst = FindExpNsKPs(
        #   k(), n.srvy.prs(), exp.N.fin(), lambda(), fnl.year(), srvy.yrs(), 
        #   phi(), rho(), ns.caps, alpha()
        # )
        # exp.ns.kps.cap.wtn.arr[hist.ind, , ] = t(exp.ns.kps.cap.lst$wtn)
        # exp.ns.kps.cap.btn.arr[hist.ind, , ] = t(exp.ns.kps.cap.lst$btn)
        
        # Find numbers of same-mother/father pairs in the population including
        # animals born in each year in the population history
        ns.kps.t.arr[hist.ind, , ] = 
          FindNsKPsT(pop.cap.hist, hist.len(), n.yrs.chk.t())
        
        # Increment progress-bar
        incProgress(1/n.sims())
      }
    }, value = 0, message = "Checking simulations")
    
    checks.lst(list(
      prpn.prnts.unkn.vec = prpn.prnts.unkn.vec,
      # ns.caps.mat = ns.caps.mat,
      ns.kps.pop.wtn.arr = ns.kps.pop.wtn.arr,
      ns.kps.pop.btn.arr = ns.kps.pop.btn.arr,
      # ns.kps.cap.wtn.arr = ns.kps.cap.wtn.arr,
      # ns.kps.cap.btn.arr = ns.kps.cap.btn.arr,
      # exp.ns.kps.cap.wtn.arr = exp.ns.kps.cap.wtn.arr,
      # exp.ns.kps.cap.btn.arr = exp.ns.kps.cap.btn.arr,
      kps.t.arr = ns.kps.t.arr
    ))
  }
})

# Estimated numbers of kin-pairs for whole population
est.ns.kps.pop.lst = reactive({
  FindEstNsKPsPop(
    exp.N.t(), s.yr.inds(), phi(), lambda(), alpha(), srvy.yrs(), k()
  )
})

# Temporally estimated numbers of kin-pairs
est.ns.kps.t = reactive({
  FindEstNsKPsT(
    exp.N.fin(), phi(), lambda(), alpha(), hist.len(), exp.N.t(), n.yrs.chk.t()
  )
})

