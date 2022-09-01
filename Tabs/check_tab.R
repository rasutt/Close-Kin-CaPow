# Global variables and list of checks objects for checks sub-tabs

# Types of kin-pairs to be displayed
kp.tps = c(
  "Population sizes", "All-pairs", "Self-pairs", "Parent-offspring pairs", 
  "Same-mother pairs", "Same-father pairs", "Full-sibling pairs", 
  "Half-sibling pairs"
)
kp.tps.pop.wtn = kp.tps[c(1:2, 5:8)]
kp.tps.pop.btn = kp.tps[c(2:3, 5)]
kp.tps.prbs.wtn = kp.tps[5:8]
kp.tps.prbs.btn = kp.tps[c(3, 5)]
kp.tps.cap.wtn = kp.tps[c(4:5, 8)]
kp.tps.cap.btn = kp.tps[3:4]

kp.tps.t = c(
  "SMP{t,f.yr,f.yr}", "SFP{t,f.yr,f.yr}", "SFP{t,t,f.yr}", 
  "SMP{t,f.yr,tm2,tm1}", "SMP{t,f.yr,tm1,f.yrm1}",
  "SMP{t,f.yr,tm1,btwn.t.f.yr}", "SMP{t,f.yr,fst.yr,btwn}"
)

# Numbers of types of kin-pairs
n.kp.tps.pop.wtn = length(kp.tps.pop.wtn)
n.kp.tps.pop.btn = length(kp.tps.pop.btn)
n.kp.tps.cap.wtn = length(kp.tps.cap.wtn)
n.kp.tps.cap.btn = length(kp.tps.cap.btn)
n.kp.tps.prb.wtn = length(kp.tps.prbs.wtn)
n.kp.tps.prb.btn = length(kp.tps.prbs.btn)
n.kp.tps.t = length(kp.tps.t)

# Calculate checks for simulated studies
checks.lst = reactiveVal(checks)
observeEvent(input$simulate, {
  # Objects to store results
  N.t.mat = matrix(nrow = n_sims(), ncol = hist.len())
  ns.kps.t.arr = array(dim = c(n_sims(), hist.len() - 2, n.kp.tps.t))
  # ns.caps.mat = ns.clvng.caps.mat = ns.clvng.mat = 
  #   matrix(nrow = n_sims(), ncol = k())
  ns.kps.pop.wtn.arr = array(dim = c(n_sims(), k(), n.kp.tps.pop.wtn))
  ns.kps.pop.btn.arr = array(dim = c(n_sims(), n.srvy.prs(), n.kp.tps.pop.btn))
  ns.kps.cap.wtn.arr = exp.ns.kps.cap.wtn.arr = 
    array(dim = c(n_sims(), k(), n.kp.tps.cap.wtn))
  ns.kps.cap.btn.arr = exp.ns.kps.cap.btn.arr = 
    array(dim = c(n_sims(), n.srvy.prs(), n.kp.tps.cap.btn))
  prpn.prnts.unkn.vec = Ns.vec = numeric(n_sims())
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n_sims()) {
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
      #   k(), n.srvy.prs(), exp.N.fin(), lambda(), f.year(), srvy.yrs(), phi(), 
      #   rho(), ns.caps, alpha()
      # )
      # exp.ns.kps.cap.wtn.arr[hist.ind, , ] = t(exp.ns.kps.cap.lst$wtn)
      # exp.ns.kps.cap.btn.arr[hist.ind, , ] = t(exp.ns.kps.cap.lst$btn)
      
      # Find numbers of same-mother/father pairs in the population including
      # animals born in each year in the population history
      # ns.kps.t.arr[hist.ind, , ] = t(FindNsKPsT(pop.cap.hist, hist.len()))
      
      # Increment progress-bar
      incProgress(1/n_sims())
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
    ns.kps.pop.btn.arr = ns.kps.pop.btn.arr
    # ns.kps.cap.wtn.arr = ns.kps.cap.wtn.arr,
    # ns.kps.cap.btn.arr = ns.kps.cap.btn.arr,
    # exp.ns.kps.cap.wtn.arr = exp.ns.kps.cap.wtn.arr,
    # exp.ns.kps.cap.btn.arr = exp.ns.kps.cap.btn.arr,
    # kps.t.arr = ns.kps.t.arr
  ))
})

