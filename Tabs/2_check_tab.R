# Objects for checks sub-tabs

## Reactives returning functions to find numbers of kin-pairs ----

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
      ncol = n.srvy.prs(), dimnames = list(NULL, Survey_pair = srvy.prs()), 
      byrow = T
    )
  }, 
  value = 0, 
  message = paste("Finding", KP.type, "pairs between survey-years"))
})

# Between pairs of survey-years
find.KPs.age.known.btn = reactive(
  function(find.func, KP.type) {
  # Need to specify matrix in case just one column or sapply returns a vector
  withProgress({
    matrix(
      sapply(sim.lst()$hists.lst, function(hist) {
        # Increment progress-bar
        incProgress(1/n.sims())
        
        # Find kin-pairs
        find.func(attributes(hist), k(), fnl.year(), srvy.yrs())
      }),
      ncol = n.srvy.prs(), dimnames = list(NULL, Survey_pair = srvy.prs()), 
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

# Nullify things depending on last simulations ----
bindEvent(observe({
  checks.lst(NULL)
  N.t.mat(NULL)
  avg.phi.obs(NULL)
  ns.SPs(NULL)
  ns.POPs(NULL)
  ns.SMPs(NULL)
  ns.SFPs(NULL)
  ns.SMPs.t(NULL)
  ns.SFPs.t(NULL)
  ns.SibPs(NULL)
  pns.UPs(NULL)
}), input$simulate)

# Compute separate checks ----
observeEvent({
  input$check.sub.tabs
  input$nav.tab
}, {
  if (input$nav.tab == "check.tab") {
    # Find population sizes over time
    if (
      input$check.sub.tabs %in% c("populations", "all.pairs", "bias.tab") && 
      is.null(N.t.mat())
    ) {
      withProgress({
        N.t.mat.new = t(sapply(sim.lst()$hists.lst, function(hist) {
          # Increment progress-bar
          incProgress(1/n.sims())
          
          # Record population curve
          attributes(hist)$N.t.vec
        }))}, value = 0, message = "Finding population sizes")
      dimnames(N.t.mat.new) = list(NULL, Year = fst.year():fnl.year())
      N.t.mat(N.t.mat.new)
    }
    
    # Find observed survival rates
    if (
      input$check.sub.tabs %in% c("demo.tab", "bias.tab") && 
      is.null(avg.phi.obs())
    ) {
      withProgress({
        avg.phi.obs(sapply(sim.lst()$hists.lst, function(hist) {
          # Increment progress-bar
          incProgress(1/n.sims())
          
          # Record population curve
          attributes(hist)$avg.phi.obs
        }))
      }, value = 0, message = "Finding population sizes")
    }
    
    # Find numbers of self-pairs between survey-years in simulated populations
    if (
      is.null(ns.SPs()) &&
      input$check.sub.tabs %in% 
      c("SPs.tab", "SMPs.tab", "SibPs.tab", "bias.tab")
    ) {
      ns.SPs(list(find.KPs.btn()(find.SPs, "self-pairs")))
    }

    # Find proportions of individuals with unknown parents
    if (
      input$check.sub.tabs %in% 
      c("POPs.tab", "SMPs.tab", "SFPs.tab", "SibPs.tab", "bias.tab") && 
      is.null(pns.UPs())
    ) {
      # Update reactive value
      pns.UPs(list(
        # In survey-years
        pns.UPs.wtn = find.KPs.wtn()(find.pns.UPs.wtn, "unknown parents"),
        
        # In pairs of survey-years
        pns.UPs.btn = find.KPs.btn()(find.pns.UPs.btn, "unknown parents")
      ))
    }
    
    # Find numbers of parent-offspring pairs
    if (
      input$check.sub.tabs %in% c("POPs.tab", "bias.tab") && 
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
      is.null(ns.SMPs()) &&
      input$check.sub.tabs %in% c("SMPs.tab", "SibPs.tab", "bias.tab")
    ) {
      # In survey-years
      ns.SMPs(list(
        wtn = find.KPs.wtn()(find.SMPs.wtn, "same-mother"),
        btn = find.KPs.btn()(find.SMSPs.btn, "same-mother and self-pairs") - 
          ns.SPs()
      ))
    }
    if (
      is.null(ns.SMPs.t()) &&
      input$check.sub.tabs %in% c("SMPs.tab", "bias.tab")
    ) {
      # Update reactive value
      ns.SMPs.t(list(
        # Same-mother pairs with ages known, in final year
        ns.SMPs.age.knwn = find.KPs.t()(find.SMPs.age.knwn, "same-mother"),
        
        # In survey-years
        ns.SMPs.wtn = ns.SMPs()$wtn,
        
        # Between pairs of survey-years
        ns.SMPs.btn = ns.SMPs()$btn,
        
        # Between survey-years, ages known
        ns.SMPs.age.known.btn = 
          find.KPs.age.known.btn()(find.SMPs.age.known.btn, "same-mother")
      ))
    }
    
    # Find numbers of same-father pairs
    if (
      is.null(ns.SFPs()) &&
      input$check.sub.tabs %in% c("SFPs.tab", "SibPs.tab", "bias.tab")
    ) {
      # In survey-years
      ns.SFPs(list(
        wtn = find.KPs.wtn()(find.SFPs.wtn, "same-father"),
        btn = find.KPs.btn()(find.SFSPs.btn, "same-father and self-pairs") - 
          ns.SPs()
      ))
    }
    if (
      is.null(ns.SFPs.t()) &&
      input$check.sub.tabs %in% c("SFPs.tab", "bias.tab") 
    ) {
      # Update reactive value
      ns.SFPs.t(list(
        # Same-father pairs with ages known, in final year
        ns.SFPs.age.knwn = find.KPs.t()(find.SFPs.age.knwn, "same-father"),
        
        # Same-father pairs with same ages known, in final year
        ns.SFPs.same.age = find.KPs.t()(find.SFPs.same.age, "same-father"),
        
        # In survey-years
        ns.SFPs.wtn = ns.SFPs()$wtn,
        
        # Between pairs of survey-years
        ns.SFPs.btn = ns.SFPs()$btn
      ))
    }
    
    # Find numbers of full and half-sibling pairs
    if (
      is.null(ns.SibPs()) &&
      input$check.sub.tabs %in% c("SibPs.tab", "bias.tab")
    ) {
      # Full-sibling pairs
      ns.FSPs.wtn = find.KPs.wtn()(find.FSPs.wtn, "full-sibling")
      ns.FSPs.btn = find.KPs.btn()(find.FSPs.btn, "full-sibling")

      # Update reactive value
      ns.SibPs(list(
        ns.FSPs.wtn = ns.FSPs.wtn,
        ns.FSPs.btn = ns.FSPs.btn,
        
        # Half-sibling pairs
        ns.HSPs.wtn = ns.SMPs()$wtn + ns.SFPs()$wtn - 2 * ns.FSPs.wtn,
        ns.HSPs.btn = ns.SMPs()$btn + ns.SFPs()$btn - 2 * ns.FSPs.btn
      ))
    }
  }
})

# Population sizes in survey-years
N.s.yrs = reactive({
  N.s.yrs.new = N.t.mat()[, s.yr.inds()]

  # Change label for column-names from Year to Survey_year
  dimnames(N.s.yrs.new) = list(NULL, Survey_year = srvy.yrs())
  N.s.yrs.new
})

# Total numbers of pairs within and between survey-years
ns.APs = reactive(list(
  wtn = choose(N.s.yrs(), 2),
  
  # If there is only one survey-pair apply returns a vector so we have to make
  # it a matrix explicitly 
  btn = t(matrix(
    apply(N.s.yrs(), 1, combn, 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2]),
    nrow = n.srvy.prs(), dimnames = list(Survey_pair = srvy.prs(), NULL)
  ))
))

# Full and half-sibling pairs
ns.FSPs = reactive(ns.SibPs()[1:2])
ns.HSPs = reactive(ns.SibPs()[3:4])

# Compute combined checks ----
observeEvent({
  input$check.sub.tabs
  input$nav.tab
}, {
  # if (
  #   input$nav.tab == "check.tab" &&
  #   !(input$check.sub.tabs %in%
  #     c("populations", "all.pairs", "SPs.tab", "POPs.tab", "SMPs.tab",
  #       "SFPs.tab", "SibPs.tab")) &&
  #   is.null(checks.lst())
  # ) {
  #   # Objects to store results
  #   # ns.caps.mat = ns.clvng.caps.mat = ns.clvng.mat =
  #   #   matrix(nrow = n.sims(), ncol = k())
  #   # Ns.vec = numeric(n.sims())
  #   # ns.kps.cap.wtn.arr = exp.ns.kps.cap.wtn.arr = array(
  #   #   dim = c(n.sims(), k(), n.kp.tps.cap.wtn),
  #   #   dimnames = list(NULL, Survey = srvy.yrs(), kp.type = kp.tps.cap.wtn)
  #   # )
  #   # ns.kps.cap.btn.arr = exp.ns.kps.cap.btn.arr = array(
  #   #   dim = c(n.sims(), n.srvy.prs(), n.kp.tps.cap.btn),
  #   #   dimnames = list(NULL, Survey_pair = srvy.prs(), kp.type = kp.tps.cap.btn)
  #   # )
  # 
  #   # Loop over histories
  #   withProgress({
  #     for (hist.ind in 1:n.sims()) {
  #       # Get simulated family and capture histories of population of animals
  #       # over time
  #       pop.cap.hist = sim.lst()$hists.lst[[hist.ind]]
  # 
  #       # Insert population sizes in survey years for comparison with numbers of
  #       # kin pairs
  #       # ns.kps.pop.wtn.arr[hist.ind, , 1] =
  #       #   attributes(pop.cap.hist)$N.t.vec[s.yr.inds()]
  # 
  #       # Record super-population size
  #       # Ns.vec[hist.ind] = attributes(pop.cap.hist)$Ns
  # 
  #       # Get numbers captured and calving in each survey
  #       # ns.caps = attributes(pop.cap.hist)$ns.caps
  #       # ns.caps.mat[hist.ind, ] = ns.caps
  #       # ns.clvng.mat[hist.ind, ] = attributes(pop.cap.hist)$ns.clvng
  #       # ns.clvng.caps.mat[hist.ind, ] = colSums(
  #       #   pop.cap.hist[, 4:(3 + k())] * pop.cap.hist[, (4 + k()):(3 + 2 * k())]
  #       # )
  # 
  #       # # Find numbers of kin pairs among samples
  #       # ns.kps.cap.lst = FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
  #       # ns.kps.cap.wtn.arr[hist.ind, , ] = t(ns.kps.cap.lst$wtn)
  #       # ns.kps.cap.btn.arr[hist.ind, , ] = t(ns.kps.cap.lst$btn)
  #       #
  #       # # Find expected numbers of kin pairs among samples
  #       # exp.ns.kps.cap.lst = FindExpNsKPs(
  #       #   k(), n.srvy.prs(), exp.N.fin(), lambda(), fnl.year(), srvy.yrs(),
  #       #   phi(), rho(), ns.caps, alpha()
  #       # )
  #       # exp.ns.kps.cap.wtn.arr[hist.ind, , ] = t(exp.ns.kps.cap.lst$wtn)
  #       # exp.ns.kps.cap.btn.arr[hist.ind, , ] = t(exp.ns.kps.cap.lst$btn)
  # 
  #       # Increment progress-bar
  #       incProgress(1/n.sims())
  #     }
  #   }, value = 0, message = "Checking simulations")
  # 
  #   checks.lst(list(
  #     # ns.caps.mat = ns.caps.mat,
  #     # ns.kps.cap.wtn.arr = ns.kps.cap.wtn.arr,
  #     # ns.kps.cap.btn.arr = ns.kps.cap.btn.arr,
  #     # exp.ns.kps.cap.wtn.arr = exp.ns.kps.cap.wtn.arr,
  #     # exp.ns.kps.cap.btn.arr = exp.ns.kps.cap.btn.arr,
  #   ))
  # }
})

