# Objects for checks sub-tabs

# Nullify checks for simulated studies when new studies simulated
observeEvent(input$simulate, {
  checks.lst(NULL)
  N.t.mat(NULL)
  ns.SPs(NULL)
  ns.POPs(NULL)
  ns.SMPs(NULL)
  ns.SFPs(NULL)
})

# Compute separate checks ----
observeEvent({
  input$check.sub.tabs
  input$nav.tab
}, {
  # Find population sizes over time
  if (
    input$check.sub.tabs %in% c("populations", "all.pairs") && 
    is.null(N.t.mat())
  ) {
    # Object to store results
    N.t.mat.new = matrix(
      nrow = n.sims(), ncol = hist.len(), 
      dimnames = list(NULL, fst.year():fnl.year())
    )
    
    # Loop over histories
    withProgress({
      for (hist.ind in 1:n.sims()) {
        # Record population curve
        N.t.mat.new[hist.ind, ] = 
          attributes(sim.lst()$hists.lst[[hist.ind]])$N.t.vec
        
        # Increment progress-bar
        incProgress(1/n.sims())
      }
    }, value = 0, message = "Finding population sizes")
    
    # Update reactive value
    N.t.mat(N.t.mat.new)
  }
  
  # Find numbers of self-pairs between survey-years in simulated populations
  if (
    input$check.sub.tabs == "self.pairs" && 
    is.null(ns.SPs())
  ) {
    # Object to store results
    ns.SPs.new = matrix(
      nrow = n.sims(), ncol = n.srvy.prs(), 
      dimnames = list(NULL, srvy.prs())
    )
    
    # Loop over histories
    withProgress({
      for (hist.ind in 1:n.sims()) {
        # Which animals alive in survey years 
        alv.s.yrs = attributes(sim.lst()$hists.lst[[hist.ind]])$alv.s.yrs
        
        # Self-pairs between survey years
        ID = attributes(sim.lst()$hists.lst[[hist.ind]])$ID
        ns.SPs.new[hist.ind, ] = as.vector(combn(1:k(), 2, function(s.inds) {
          sum(ID[alv.s.yrs[, s.inds[1]]] %in% ID[alv.s.yrs[, s.inds[2]]])
        }))
        
        # Increment progress-bar
        incProgress(1/n.sims())
      }
    }, value = 0, message = "Finding self-pairs")
    
    # Update reactive value
    ns.SPs(ns.SPs.new)
  }
  
  # Find numbers of parent-offspring pairs between survey-years in simulated
  # populations
  if (
    input$check.sub.tabs == "POPs.tab" && 
    is.null(ns.POPs())
  ) {
    # Objects to store results
    ns.POPs.wtn = matrix(
      nrow = n.sims(), ncol = k(), 
      dimnames = list(NULL, srvy.yrs())
    )
    ns.POPs.btn = matrix(
      nrow = n.sims(), ncol = n.srvy.prs(), 
      dimnames = list(NULL, srvy.prs())
    )
    
    # Loop over histories
    withProgress({
      for (hist.ind in 1:n.sims()) {
        # Which animals alive in survey years 
        alv.s.yrs = attributes(sim.lst()$hists.lst[[hist.ind]])$alv.s.yrs
        
        ID = attributes(sim.lst()$hists.lst[[hist.ind]])$ID
        mum = attributes(sim.lst()$hists.lst[[hist.ind]])$mum
        dad = attributes(sim.lst()$hists.lst[[hist.ind]])$dad
        
        # Parent-offspring pairs within survey years
        ns.POPs.wtn[hist.ind, ] = sapply(1:k(), function(s.ind) {
          sum(
            mum[alv.s.yrs[, s.ind]] %in% ID[alv.s.yrs[, s.ind]], 
            dad[alv.s.yrs[, s.ind]] %in% ID[alv.s.yrs[, s.ind]]
          )
        })
        
        # Parent-offspring pairs between survey years
        ns.POPs.btn[hist.ind, ] = as.vector(combn(1:k(), 2, function(s.inds) {
          sum(
            mum[alv.s.yrs[, s.inds[1]]] %in% ID[alv.s.yrs[, s.inds[2]]], 
            dad[alv.s.yrs[, s.inds[1]]] %in% ID[alv.s.yrs[, s.inds[2]]],
            mum[alv.s.yrs[, s.inds[2]]] %in% ID[alv.s.yrs[, s.inds[1]]], 
            dad[alv.s.yrs[, s.inds[2]]] %in% ID[alv.s.yrs[, s.inds[1]]]
          )
        }))
        
        # Increment progress-bar
        incProgress(1/n.sims())
      }
    }, value = 0, message = "Finding self-pairs")

    # Update reactive value
    ns.POPs(list(ns.POPs.wtn = ns.POPs.wtn, ns.POPs.btn = ns.POPs.btn))
  }
  
  # Find numbers of same-mother pairs between survey-years in simulated
  # populations
  if (
    input$check.sub.tabs == "SMPs.tab" && 
    is.null(ns.SMPs())
  ) {
    # Objects to store results
    ns.SMPs.wtn = matrix(
      nrow = n.sims(), ncol = k(), 
      dimnames = list(NULL, srvy.yrs())
    )
    ns.SMPs.btn = matrix(
      nrow = n.sims(), ncol = n.srvy.prs(), 
      dimnames = list(NULL, srvy.prs())
    )
    
    # Loop over histories
    withProgress({
      for (hist.ind in 1:n.sims()) {
        # Which animals alive in survey years 
        alv.s.yrs = attributes(sim.lst()$hists.lst[[hist.ind]])$alv.s.yrs
        
        # List of frequency tables of mums in each survey year
        mum = attributes(sim.lst()$hists.lst[[hist.ind]])$mum
        max.mum = max(mum, na.rm = T)
        mum.tab.lst = lapply(1:k(), function(s.ind) {
          tabulate(mum[alv.s.yrs[, s.ind]], max.mum)
        })
        
        # Same-mother pairs within survey years
        ns.SMPs.wtn[hist.ind, ] = sapply(1:k(), function(s.ind) {
          sum(choose(mum.tab.lst[[s.ind]], 2))
        })
        
        # Same-mother pairs between survey years
        ns.SMPs.btn[hist.ind, ] = as.vector(combn(1:k(), 2, function(s.inds) {
          mum.tab.lst[[s.inds[1]]] %*% mum.tab.lst[[s.inds[2]]]
        }))
        
        # Increment progress-bar
        incProgress(1/n.sims())
      }
    }, value = 0, message = "Finding self-pairs")
    
    # Update reactive value
    ns.SMPs(list(ns.SMPs.wtn = ns.SMPs.wtn, ns.SMPs.btn = ns.SMPs.btn))
  }
  
  # Find numbers of same-father pairs between survey-years in simulated
  # populations
  if (
    input$check.sub.tabs == "SFPs.tab" && 
    is.null(ns.SFPs())
  ) {
    # Objects to store results
    ns.SFPs.wtn = matrix(
      nrow = n.sims(), ncol = k(), 
      dimnames = list(NULL, srvy.yrs())
    )
    ns.SFPs.btn = matrix(
      nrow = n.sims(), ncol = n.srvy.prs(), 
      dimnames = list(NULL, srvy.prs())
    )
    
    # Loop over histories
    withProgress({
      for (hist.ind in 1:n.sims()) {
        # Which animals alive in survey years 
        alv.s.yrs = attributes(sim.lst()$hists.lst[[hist.ind]])$alv.s.yrs
        
        # List of frequency tables of dads in each survey year
        dad = attributes(sim.lst()$hists.lst[[hist.ind]])$dad
        max.dad = max(dad, na.rm = T)
        dad.tab.lst = lapply(1:k(), function(s.ind) {
          tabulate(dad[alv.s.yrs[, s.ind]], max.dad)
        })
        
        # Same-father pairs within survey years
        ns.SFPs.wtn[hist.ind, ] = sapply(1:k(), function(s.ind) {
          sum(choose(dad.tab.lst[[s.ind]], 2))
        })
        
        # Same-father pairs between survey years
        # ns.SFPs.btn[hist.ind, ] = as.vector(combn(1:k(), 2, function(s.inds) {
        #   dad.tab.lst[[s.inds[1]]] %*% dad.tab.lst[[s.inds[2]]]
        # }))
        
        # Increment progress-bar
        incProgress(1/n.sims())
      }
    }, value = 0, message = "Finding self-pairs")
    
    # Update reactive value
    ns.SFPs(list(ns.SFPs.wtn = ns.SFPs.wtn, ns.SFPs.btn = ns.SFPs.btn))
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
    !(input$check.sub.tabs %in% 
      c("populations", "all.pairs", "self.pairs", "POPs.tab", "SMPs.tab", 
        "SFPs.tab")) && 
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
    ns.kps.cap.wtn.arr = exp.ns.kps.cap.wtn.arr = array(
      dim = c(n.sims(), k(), n.kp.tps.cap.wtn),
      dimnames = list(NULL, Survey = srvy.yrs(), kp.type = kp.tps.cap.wtn)
    )
    ns.kps.cap.btn.arr = exp.ns.kps.cap.btn.arr = array(
      dim = c(n.sims(), n.srvy.prs(), n.kp.tps.cap.btn),
      dimnames = list(NULL, Survey_pair = srvy.prs(), kp.type = kp.tps.cap.btn)
    )
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

