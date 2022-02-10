# Display parameter values
output$checkParVals <- renderTable({
  par_vals_df(sim.par.vals(), sim.par.names())
}, digits = 3)
# Display simulation values
output$checkSimVals = renderTable(sim.vals())

# Calculate checks for simulated studies
checks.lst = reactive({
  # Objects to store results
  N.t.mat <- matrix(nrow = n_sims(), ncol = hist.len())
  SMFPs.t.arr = array(dim = c(n_sims(), hist.len() - 2, 7))
  # ns.caps.mat <- ns.clvng.caps.mat <- ns.clvng.mat <- 
  #   matrix(nrow = n_sims(), ncol = k())
  ns.KPs.wtn.pop.arr <- array(dim = c(n_sims(), k(), 6))
  ns.KPs.btn.pop.arr <- array(dim = c(n_sims(), n.srvy.prs(), 4))
  ns.KPs.wtn.cap.arr <- exp.ns.KPs.wtn.cap.arr <- 
    array(dim = c(n_sims(), k(), 3))
  ns.KPs.btn.cap.arr <- exp.ns.KPs.btn.cap.arr <- 
    array(dim = c(n_sims(), n.srvy.prs(), 2))
  prpn.prnts.unkn.vec <- Ns.vec <- numeric(n_sims())
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n_sims()) {
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hist.ind]]
      
      # Record population curve
      N.t.mat[hist.ind, ] <- attributes(pop.cap.hist)$N.t.vec
      
      # Record super-population size
      Ns.vec[hist.ind] <- attributes(pop.cap.hist)$Ns
      
      # Get numbers captured and calving in each survey
      ns.caps <- attributes(pop.cap.hist)$ns.caps
      # ns.caps.mat[hist.ind, ] <- ns.caps
      # ns.clvng.mat[hist.ind, ] <- attributes(pop.cap.hist)$ns.clvng
      # ns.clvng.caps.mat[hist.ind, ] <- colSums(
      #   pop.cap.hist[, 4:(3 + k())] * pop.cap.hist[, (4 + k()):(3 + 2 * k())]
      # )
      
      # Find proportion captured with unknown parents
      prpn.prnts.unkn.vec[hist.ind] <- mean(is.na(pop.cap.hist$mum))
      
      # Numbers of kin-pairs in whole population
      ns.kps.pop.lst = FindNsKinPairsPop(pop.cap.hist, s.yr.inds(), k())
      ns.KPs.wtn.pop.arr[hist.ind, , -1] = t(ns.kps.pop.lst$wtn)
      ns.KPs.btn.pop.arr[hist.ind, , ] = t(ns.kps.pop.lst$btn)
      
      # Find numbers of kin pairs among samples
      ns.kps.cap.lst <- FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
      ns.KPs.wtn.cap.arr[hist.ind, , ] = t(ns.kps.cap.lst$wtn)
      ns.KPs.btn.cap.arr[hist.ind, , ] = t(ns.kps.cap.lst$btn)
      
      # Find expected numbers of kin pairs among samples
      exp.ns.kps.cap.lst <- FindExpNsKPs(
        k(), n.srvy.prs(), exp.N.fin(), lambda(), f.year(), srvy.yrs(), phi(), 
        rho(), ns.caps, alpha()
      )
      exp.ns.KPs.wtn.cap.arr[hist.ind, , ] = t(exp.ns.kps.cap.lst$wtn)
      exp.ns.KPs.btn.cap.arr[hist.ind, , ] = t(exp.ns.kps.cap.lst$btn)

      # Find numbers of same-mother/father pairs in the population including
      # animals born in each year in the population history
      SMFPs.t.arr[hist.ind, , ] = t(FindSMFPsT(pop.cap.hist, hist.len()))
      
      # Increment progress-bar
      incProgress(1/n_sims())
    }
  }, value = 0, message = "Checking simulations")

  # Insert population sizes in survey years for comparison with numbers of kin
  # pairs
  ns.KPs.wtn.pop.arr[, , 1] = N.t.mat[, s.yr.inds()]

  list(
    N.t.mat = N.t.mat, 
    prpn.prnts.unkn.vec = prpn.prnts.unkn.vec,
    # ns.caps.mat = ns.caps.mat,
    ns.KPs.wtn.pop.arr = ns.KPs.wtn.pop.arr,
    ns.KPs.btn.pop.arr = ns.KPs.btn.pop.arr,
    ns.KPs.wtn.cap.arr = ns.KPs.wtn.cap.arr,
    ns.KPs.btn.cap.arr = ns.KPs.btn.cap.arr,
    exp.ns.KPs.wtn.cap.arr = exp.ns.KPs.wtn.cap.arr,
    exp.ns.KPs.btn.cap.arr = exp.ns.KPs.btn.cap.arr,
    SMFPs.t.arr = SMFPs.t.arr
  )
})

# Show check results

# Plot population sizes over time
output$popPlot <- renderPlot({
  # Plot population trajectories
  matplot(
    (f.year() - hist.len() + 1):f.year(), t(checks.lst()$N.t.mat), type = 'l',
    col = rgb(0, 0, 0, alpha = 0.1), lty = 1, 
    ylim = c(0, max(checks.lst()$N.t.mat)),
    xlab = 'Year', ylab = 'Nt', main = "Population sizes over time"
  )
  # Expected value
  lines((f.year() - hist.len() + 1):f.year(), exp.N.t(), col = 2, lwd = 2)
  # Average value
  lines(
    (f.year() - hist.len() + 1):f.year(), colMeans(checks.lst()$N.t.mat), 
    col = 4, lwd = 2
  )
  # Surveys
  abline(v = srvy.yrs(), lty = 2)
  # Add legend
  legend(
    "topleft", legend = c("Expected", "Simulated", "Mean", "Survey years"),
    col = c(2, rgb(0, 0, 0, alpha = 0.1), 4, 1), lwd = c(2, 1, 2, 1), 
    lty = c(1, 1, 1, 2)
  )
})

# Table simulated versus expected self-pair probability
output$slfPrPrb = renderTable({
  spr = mean(checks.lst()$ns.KPs.pop.lst$ns.SPs.btn.pop.mat /
               checks.lst()$ns.KPs.pop.lst$ns.APs.btn.pop.mat)
  exp_spr = mean(exp.ns.KPs.pop.lst()$exp.ns.SPs.btn /
                   exp.ns.KPs.pop.lst()$exp.ns.APs.btn)
  df = data.frame(c(spr, exp_spr))
  names(df) = c("Self-pair rate")
  rownames(df) = c("Average simulated", "Average Expected")
  df
}, rownames = T, digits = 6)

# Checks for numbers of kin-pairs

# Show percentage of animals captured for which the parents are unknown
output$percUnknPrnts <- renderText({
  paste0(
    "Percentage of captured animals with unknown parents: ", 
    round(mean(checks.lst()$prpn.prnts.unkn.vec) * 100, 1), "%"
  )
})

# Types of kin-pairs
kn.pr.wtn.pop.tps = c(
  "Population sizes", "All-pairs", "Same-mother pairs", "Same-father pairs",
  "Full-sibling pairs", "Half-sibling pairs"
)
kn.pr.btn.pop.tps = c("All-pairs", "Self-pairs", "Same-mother pairs")
kn.pr.prb.tps = c("Self-pair", "Same-mother pair")
kn.pr.wtn.cap.tps = 
  c("Parent-offspring pairs", "Same-mother pairs", "Half-sibling pairs")
kn.pr.btn.cap.tps = c("Self-pairs", "Parent-offspring pairs")

# Numbers of types of kin-pairs
n.kn.pr.wtn.pop.tps = length(kn.pr.wtn.pop.tps)
n.kn.pr.btn.pop.tps = length(kn.pr.btn.pop.tps)
n.kn.pr.wtn.cap.tps = length(kn.pr.wtn.cap.tps)
n.kn.pr.btn.cap.tps = length(kn.pr.btn.cap.tps)
n.kn.pr.prb.tps = length(kn.pr.prb.tps)

# Indices for types of kin-pairs
KP_pop_inds = c(
  "Survey", "Survey", "Survey-pair", "Survey-pair", "Survey", "Survey", 
  "Survey", "Survey", "Survey-pair"
)
KP_prob_inds = c("Survey-pair", "Survey")
KP_inds = c("Survey-pair", "Survey", "Survey-pair", "Survey", "Survey")

# Expected numbers of kin-pairs for whole population
exp.ns.KPs.pop.lst = reactive({
  # Expected population sizes in survey years
  exp.N.srvy.yrs = exp.N.t()[s.yr.inds()]
  # Approximation for expected number of pairs of animals in survey years
  exp.ns.APs.wtn = choose(exp.N.srvy.yrs, 2)
  # All pairs within surveys and between pairs of surveys
  exp.ns.APs.btn = as.vector(combn(exp.N.srvy.yrs, 2, function(x) x[1] * x[2]))
  
  # Find kin-pair probabilities
  kp.prbs.lst = FindKPProbs(
    exp.N.srvy.yrs, exp.ns.APs.wtn, phi(), lambda(), alpha(), srvy.yrs(), k()
  )
  exp.ns.kps.wtn = t(t(kp.prbs.lst$wtn) * exp.ns.APs.wtn)
  exp.ns.kps.btn = t(t(kp.prbs.lst$btn) * exp.ns.APs.btn)
  
  # Return as list
  list(
    exp.N.srvy.yrs = exp.N.srvy.yrs,
    exp.ns.APs.wtn = exp.ns.APs.wtn,
    exp.ns.APs.btn = exp.ns.APs.btn,
    exp.ns.SPs.btn = exp.ns.kps.btn[1, ],
    exp.ns.SMPs.wtn = exp.ns.kps.wtn[1, ],
    exp.ns.SFPs.wtn = exp.ns.kps.wtn[2, ],
    exp.ns.FSPs.wtn = exp.ns.kps.wtn[3, ],
    exp.ns.HSPs.wtn = exp.ns.kps.wtn[4, ],
    exp.ns.SMPs.btn = exp.ns.kps.btn[2, ]
  )
})

# Kin-pair probabilities
KPs.prob.lst = reactive({
  list(
    SPs.btn = exp.ns.KPs.pop.lst()$exp.ns.SPs.btn / 
      exp.ns.KPs.pop.lst()$exp.ns.APs.btn,
    SMPs.wtn = exp.ns.KPs.pop.lst()$exp.ns.SMPs.wtn / 
      exp.ns.KPs.pop.lst()$exp.ns.APs.wtn
  )
})

# Kin-pair rates observed
KPs.rate.lst = reactive({
  list(
    SPs.btn = checks.lst()$ns.KPs.pop.lst$ns.SPs.btn.pop.mat / 
      checks.lst()$ns.KPs.pop.lst$ns.APs.btn.pop.mat,
    SMPs.wtn = checks.lst()$ns.KPs.pop.lst$ns.SMPs.wtn.pop.mat / 
      checks.lst()$ns.KPs.pop.lst$ns.APs.wtn.pop.mat
  )
})

# Function to table average percentage differences between simulated and
# expected numbers of kin-pairs
KPstab = function(n_types, ns, exp.ns, type_names, pop) {
  diffs = numeric(n_types)
  for (i in 1:n_types) {
    if (pop) {
      diffs[i] = perc(mean(t((t(ns[[i]]) - exp.ns[[i]]) / exp.ns[[i]])))
    } else {
      diffs[i] = perc(mean((ns[[i]] - exp.ns[[i]]) / exp.ns[[i]]))
    }
  }
  df = data.frame(matrix(diffs, nrow = 1))
  names(df) = type_names
  df
}

# Function to table average percentage differences between simulated and
# expected numbers of kin-pairs
kn.prs.tbl.arr = function(n_types, ns, exp.ns, type_names, pop) {
  diffs = numeric(n_types)
  for (i in 1:n_types) {
    if (pop) {
      diffs[i] = perc(mean(t(t(ns[, , i]) / exp.ns[[i]])) - 1)
    } else {
      diffs[i] = perc(mean(ns[, , i] / exp.ns[, , i]) - 1)
    }
  }
  df = data.frame(matrix(diffs, nrow = 1))
  names(df) = type_names
  df
}

# Table average percentage differences for whole population
output$nsKPsPopWtn = renderTable({
  kn.prs.tbl.arr(
    n.kn.pr.wtn.pop.tps, checks.lst()$ns.KPs.wtn.pop.arr, 
    exp.ns.KPs.pop.lst()[c(1:2, 5:8)], kn.pr.wtn.pop.tps, T
  )
})
output$nsKPsPopBtn = renderTable({
  kn.prs.tbl.arr(
    n.kn.pr.btn.pop.tps, checks.lst()$ns.KPs.btn.pop.arr, 
    exp.ns.KPs.pop.lst()[c(3:4, 9)], kn.pr.btn.pop.tps, T
  )
})

# Table average percentage differences for whole population
output$nsKPsProb = renderTable({
  KPstab(n.kn.pr.prb.tps, KPs.rate.lst(), KPs.prob.lst(), kn.pr.prb.tps, T)
})

# Table average percentage differences for captured animals
output$nsKPsCapWtn = renderTable({
  kn.prs.tbl.arr(
    n.kn.pr.wtn.cap.tps, checks.lst()$ns.KPs.wtn.cap.arr, 
    checks.lst()$exp.ns.KPs.wtn.cap.arr, kn.pr.wtn.cap.tps, F
  )
})
output$nsKPsCapBtn = renderTable({
  kn.prs.tbl.arr(
    n.kn.pr.btn.cap.tps, checks.lst()$ns.KPs.btn.cap.arr, 
    checks.lst()$exp.ns.KPs.btn.cap.arr, kn.pr.btn.cap.tps, F
  )
})

# Numbers of same-mother/father pairs in the population including animals born
# in each year in the population history
output$nsSMFPsT = renderTable({
  # Find average values simulated
  mean.SMFPs.t = colMeans(checks.lst()$SMFPs.t.arr)
  
  # Expected numbers
  
  # Same-mother pairs
  exp.ns.SMPs.f.yr = 2 * exp.N.t()[hist.len()] * (1 - phi() / lambda())^2 *
    (lambda() / phi())^alpha() * (phi()^2 / lambda())^((hist.len() - 2):1)
  # Same-father pairs with different ages (one always age zero)
  exp.ns.SFPs.f.yr = exp.ns.SMPs.f.yr * phi()
  # Same-father pairs with same ages (for all ages up to length of population
  # history)
  exp.ns.SFPs.t.t = phi()^(2 * (hist.len() - 2):1) * 
    (exp.N.t()[hist.len()] / lambda()^(alpha() + (hist.len() - 2):1) - 1) *
    (lambda()^2 / phi())^alpha() * lambda() * (1 - phi() / lambda())^2
  # Same-mother pairs between each year and final year, born in two years
  # preceding each year
  exp.ns.SMPs.tm2.tm1.t.f.yr = 4 * exp.N.t()[2:(hist.len() - 1)] * 
    (1 - phi() / lambda())^2 * (lambda() / phi())^alpha() * 
    (phi()^2 / lambda())^2 * phi()^((hist.len() - 2):1)
  # Same-mother pairs between each year and final year, born in years
  # preceding both
  exp.ns.SMPs.tm1.f.yrm1.t.f.yr = 2 * exp.N.t()[2:(hist.len() - 1)] * 
    (1 - phi() / lambda())^2 * (lambda() / phi())^alpha() * 
    (phi()^2 / lambda()) * phi()^((hist.len() - 2):1)
  # Same-mother pairs between each year and final year, born year before each
  # year, and all years in between each year and final year
  exp.ns.SMPs.tm1.btwn.t.f.yr.t.f.yr = 2 * exp.N.t()[2:(hist.len() - 1)] * 
    (1 - phi() / lambda())^2 * (lambda() / phi())^alpha() * 
    (phi()^2 / lambda()) * phi()^((hist.len() - 2):1) * (hist.len() - 2):1
    
  # Same-mother pairs between each year and final year, born in first year of
  # history, and all years in between first of history and each year
  exp.ns.SMPs.fst.yr.btwn.t.f.yr = 4 * exp.N.t()[2:(hist.len() - 1)] * 
    (1 - phi() / lambda())^2 * (lambda() / phi())^alpha() * 
    (phi()^2 / lambda())^(1:(hist.len() - 2)) * phi()^((hist.len() - 2):1) *
    1:(hist.len() - 2)
  
  # Combine for output
  df = rbind(
    mean.SMFPs.t[, 1], exp.ns.SMPs.f.yr, 
    mean.SMFPs.t[, 2], exp.ns.SFPs.f.yr,
    mean.SMFPs.t[, 3], exp.ns.SFPs.t.t,
    mean.SMFPs.t[, 4], exp.ns.SMPs.tm2.tm1.t.f.yr,
    mean.SMFPs.t[, 5], exp.ns.SMPs.tm1.f.yrm1.t.f.yr,
    mean.SMFPs.t[, 6], exp.ns.SMPs.tm1.btwn.t.f.yr.t.f.yr,
    mean.SMFPs.t[, 7], exp.ns.SMPs.fst.yr.btwn.t.f.yr
  )
  rownames(df) = c(
    "Avg(SMP{t,f.yr,f.yr})", "E(SMP{t,f.yr,f.yr})",
    "Avg(SFP{t,f.yr,f.yr})", "E(SFP{t,f.yr,f.yr})",
    "Avg(SFP{t,t,f.yr})", "E(SFP{t,t,f.yr})",
    "Avg(SMP{t,f.yr,tm2,tm1})", "E(SMP{t,f.yr,tm2,tm1})",
    "Avg(SMP{t,f.yr,tm1,f.yrm1})", "E(SMP{t,f.yr,tm1,f.yrm1})",
    "Avg(SMP{t,f.yr,tm1,btwn.t.f.yr})", "E(SMP{t,f.yr,tm1,btwn.t.f.yr})",
    "Avg(SMP{t,f.yr,fst.yr,btwn})", "E(SMP{t,f.yr,fst.yr,btwn})"
  )
  colnames(df) = (f.year() - hist.len() + 2):(f.year() - 1)
  df
}, rownames = T, digits = 1)

# Function to plot simulated versus expected numbers of kin-pairs for one type
# of kin-pair
nsKPsPlot = function(i, pop = F, prob = F) {
  if (pop) {
    diffs = 
      t(t(checks.lst()$ns.KPs.pop.lst[[i]]) / exp.ns.KPs.pop.lst()[[i]] - 1)
    xlab = KP_pop_inds[i]
    main = KP_pop_names[i]
  } else {
    diffs = checks.lst()$ns.KPs.lst[[i]] / checks.lst()$exp.ns.KPs.lst[[i]] - 1 
    xlab = KP_inds[i]
    main = KP_names[i]
  }
  if (prob) {
    diffs = 
      t(t(KPs.rate.lst()[[i]]) / KPs.prob.lst()[[i]] - 1)
    xlab = KP_prob_inds[i]
    main = KP_prob_names[i]
  }
  if (xlab == "Survey") colnames(diffs) = srvy.yrs()
  else colnames(diffs) = apply(combn(srvy.yrs(), 2), 2, paste, collapse = "-")
  {
    boxplot(
      diffs, main = main, xlab = xlab, 
      ylab = "(Observed - expected numbers) / expected numbers"
    )
    abline(h = 0, col = 'red')
    abline(h = mean(diffs), col = 'blue')
    legend(
      "topleft", col = c(2, 4), lty = 1,
      legend = c(
        "Expected difference over all surveys (zero)", 
        "Average difference over all surveys"
      ),
    )
  }
}

# Apply plot function to each type of kin-pair

# # All pairs within survey years for whole population
# output$nsAPsWtnPop <- renderPlot(nsKPsPlot(1, T))
# # All pairs between survey years for whole population
# output$nsAPsBtnPop <- renderPlot(nsKPsPlot(2, T))
# # Self-pairs between survey years for whole population
# output$nsSPsBtnPop <- renderPlot(nsKPsPlot(3, T))
# # Same-mother pairs within survey years for whole population
# output$nsSMPsWtnPop <- renderPlot(nsKPsPlot(4, T))
# # Same-father pairs within survey years for whole population
# output$nsSFPsWtnPop <- renderPlot(nsKPsPlot(5, T))
# # Same-mother pairs between survey years for whole population
# output$nsSMPsBtnPop <- renderPlot(nsKPsPlot(9, T))

# # Self-pair probabilities between survey years for whole population
# output$probSPsBtnPop <- renderPlot(nsKPsPlot(1, prob = T))
# # Same-mother pair probabilities within survey years for whole population
# output$probSMPsWtnPop <- renderPlot(nsKPsPlot(2, prob = T))

# # Self-pairs between samples
# output$nsSPsBtn <- renderPlot(nsKPsPlot(1))
# # Parent-offspring pairs within samples
# output$nsPOPsWtn <- renderPlot(nsKPsPlot(2))
# # Parent-offspring pairs between samples
# output$nsPOPsBtn <- renderPlot(nsKPsPlot(3))
# # Same-mother pairs within samples
# output$nsSMPsWtn <- renderPlot(nsKPsPlot(4))
# # Half-sibling pairs within samples
# output$nsHSPsWtn <- renderPlot(nsKPsPlot(5))

# First life histories from first study
output$alive = renderTable({
  alv_mat = attributes(sim.lst()$hists.lst[[1]])$alv.mat
  mode(alv_mat) = "integer"
  df = data.frame(head(alv_mat))
  names(df) = (f.year() - hist.len() + 1):f.year()
  df
}, rownames = T)

# Print head of first study
output$dataHead <- renderTable(head(data.frame(sim.lst()$hists.lst[[1]])))