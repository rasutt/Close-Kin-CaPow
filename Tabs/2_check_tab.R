# Objects for checks sub-tabs

# Calculate checks for simulated studies
observeEvent(input$simulate, {
  # Objects to store results
  N.t.mat = matrix(nrow = n.sims(), ncol = hist.len())
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
  # ns.kps.t.arr = array(dim = c(n.sims(), hist.len() - 2, n.kp.tps.t))
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
      
      # Record population curve
      N.t.mat[hist.ind, ] = attributes(pop.cap.hist)$N.t.vec
      
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
        FindNsKPsT(pop.cap.hist, hist.len(), n.kp.tps.t, n.yrs.chk.t())
      
      # Increment progress-bar
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Checking simulations")
  
  # Insert population sizes in survey years for comparison with numbers of kin
  # pairs
  ns.kps.pop.wtn.arr[, , 1] = N.t.mat[, s.yr.inds()]
  
  checks.lst(list(
    N.t.mat = N.t.mat, 
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

# Function to find estimate errors as proportions of estimates
find.est.errs = function(vals, ests, samp = F) {
  # For sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) ests = rep(ests, each = n.sims())
  vals / ests - 1
}

# Estimate errors as proportions of estimates
ns.kps.pop.wtn.est.errs = reactive({
  find.est.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn)
})
ns.kps.pop.btn.est.errs = reactive({
  find.est.errs(checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn)
})
ns.kps.t.est.errs = reactive({
  find.est.errs(checks.lst()$kps.t.arr, est.ns.kps.t())
})
ns.kps.prb.wtn.est.errs = reactive({
  # Remove population sizes and total numbers of pairs then divide by the latter
  find.est.errs(
    checks.lst()$ns.kps.pop.wtn.arr[, , -1:-2] / 
      array(
        rep(checks.lst()$ns.kps.pop.wtn.arr[, , 2], n.kp.tps.prb.wtn), 
        c(n.sims(), k(), n.kp.tps.prb.wtn)
      ), 
    est.ns.kps.pop.lst()$wtn[, -1:-2] / est.ns.kps.pop.lst()$wtn[, 2]
  )
})
ns.kps.prb.btn.est.errs = reactive({
  # Remove population sizes and total numbers of pairs then divide by the latter
  # (drop = F retains array dimensions when only one survey pair)
  find.est.errs(
    checks.lst()$ns.kps.pop.btn.arr[, , -1, drop = F] / 
      array(
        rep(checks.lst()$ns.kps.pop.btn.arr[, , 1], n.kp.tps.prb.btn), 
        c(n.sims(), n.srvy.prs(), n.kp.tps.prb.btn)
      ), 
    est.ns.kps.pop.lst()$btn[, -1] / est.ns.kps.pop.lst()$btn[, 1]
  )
})
# ns.kps.cap.wtn.est.errs = reactive({
#   find.est.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn, T)
# })
# ns.kps.cap.btn.est.errs = reactive({
#   find.est.errs(checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn, T)
# })
