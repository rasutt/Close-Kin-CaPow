# Display parameter values
output$checkParVals = renderTable({
  par_vals_df(sim.par.vals(), sim.par.names())
}, digits = 3)

# Display simulation values
output$checkSimVals = renderTable(sim.vals())

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
checks.lst = reactive({
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

  list(
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
  )
})

# Show check results

# Plot population sizes over time
output$popPlot = renderPlot({
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

# Checks for numbers of kin-pairs

# Show percentage of animals captured for which the parents are unknown
output$percUnknPrnts = renderText({
  paste0(
    "Percentage of captured animals with unknown parents: ", 
    round(mean(checks.lst()$prpn.prnts.unkn.vec) * 100, 1), "%"
  )
})

# Estimated numbers of kin-pairs for whole population
est.ns.kps.pop.lst = reactive({
  FindEstNsKPsPop(
    exp.N.t(), s.yr.inds(), phi(), lambda(), alpha(), srvy.yrs(), k()
  )
})

### Kin-pair estimator biases (tables of average percentage differences)

# Function to find estimate errors as proportions of estimates
find.est.errs = function(vals, ests, samp = F) {
  # For sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) ests = rep(ests, each = n_sims())
  vals / ests - 1
}

# Errors in estimates for numbers of kin-pairs in whole population within
# surveys
ns.kps.pop.wtn.est.errs = reactive({
  find.est.errs(checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn)
})

## Function to find them from true and expected values
est.bias = function(vals, exp.vals, type_names, samp = F) {
  # For the sampled animals the expected values are different for each study but
  # otherwise they are repeated
  if (!samp) exp.vals = rep(exp.vals, each = n_sims())
  df = data.frame(
    matrix(perc(colMeans(vals / exp.vals, dims = 2) - 1), nrow = 1)
  )
  names(df) = type_names
  df
}

## Numbers in whole population

# Within surveys
output$biasNsKPsPopWtn = renderTable({
  est.bias(
    checks.lst()$ns.kps.pop.wtn.arr, est.ns.kps.pop.lst()$wtn, kp.tps.pop.wtn
  )
})

# Between surveys
output$biasNsKPsPopBtn = renderTable({
  est.bias(
    checks.lst()$ns.kps.pop.btn.arr, est.ns.kps.pop.lst()$btn, kp.tps.pop.btn
  )
})

## Probabilities (numbers divided by total numbers of pairs)

# Within surveys
output$biasProbsKPsWtn = renderTable({
  # Remove population sizes and total numbers of pairs then divide by the latter
  est.bias(
    checks.lst()$ns.kps.pop.wtn.arr[, , -1:-2] / 
      array(
        rep(checks.lst()$ns.kps.pop.wtn.arr[, , 2], n.kp.tps.prb.wtn), 
        c(n_sims(), k(), n.kp.tps.prb.wtn)
      ), 
    est.ns.kps.pop.lst()$wtn[, -1:-2] / est.ns.kps.pop.lst()$wtn[, 2], 
    kp.tps.prbs.wtn
  )
})

# Between surveys
output$biasProbsKPsBtn = renderTable({
  est.bias(
    # Remove total numbers of pairs then divide by them
    checks.lst()$ns.kps.pop.btn.arr[, , -1] / 
      array(
        rep(checks.lst()$ns.kps.pop.btn.arr[, , 1], n.kp.tps.prb.btn), 
        c(n_sims(), n.srvy.prs(), n.kp.tps.prb.btn)
      ), 
    est.ns.kps.pop.lst()$btn[, -1] / est.ns.kps.pop.lst()$btn[, 1], 
    kp.tps.prbs.btn
  )
})

## Numbers among sampled animals

# Within surveys
output$biasNsKPsCapWtn = renderTable({
  est.bias(
    checks.lst()$ns.kps.cap.wtn.arr,
    checks.lst()$exp.ns.kps.cap.wtn.arr, kp.tps.cap.wtn,
    cap = T
  )
})

# Between surveys
output$biasNsKPsCapBtn = renderTable({
  est.bias(
    checks.lst()$ns.kps.cap.btn.arr,
    checks.lst()$exp.ns.kps.cap.btn.arr, kp.tps.cap.btn,
    cap = T
  )
})

## Temporal estimates vs observed averages (mainly for debugging)

# Numbers of same-mother/father pairs in the population including animals born
# in each year in the population history
output$nsKPsTemp = renderTable({
  # Find average values simulated
  mean.kps.t = t(apply(checks.lst()$kps.t.arr, 2:3, mean))

  # Find estimated values
  exp.kps.t = FindExpNsKPsT(
    exp.N.fin(), phi(), lambda(), alpha(), hist.len(), exp.N.t()
  )

  # Combine for output
  df = rbind(mean.kps.t, exp.kps.t)
  rownames(df) = paste0(
    rep(c("Avg", "Exp"), each = n.kp.t.tps),
    rep(kp.t.tps, 2)
  )
  colnames(df) = (f.year() - hist.len() + 2):(f.year() - 1)
  df[rep(seq(1:n.kp.t.tps), each = 2) + c(0, n.kp.t.tps), ]
}, rownames = T, digits = 1)

# Kin-pair probability estimates
exp.probs.KPs.lst = reactive({
  list(
    SPs.btn = est.ns.kps.pop.lst()$btn[2, ] / est.ns.kps.pop.lst()$btn[1, ],
    SMPs.wtn = est.ns.kps.pop.lst()$btn[3, ] / est.ns.kps.pop.lst()$btn[1, ]
  )
})

# Kin-pair probabilities observed
probs.KPs.lst = reactive({
  list(
    SPs.btn = checks.lst()$ns.kps.pop.btn.arr[, , 2] / 
      checks.lst()$ns.kps.pop.btn.arr[, , 1],
    SMPs.wtn = checks.lst()$ns.kps.pop.btn.arr[, , 3] / 
      checks.lst()$ns.kps.pop.btn.arr[, , 1]
  )
})

# Indices for types of kin-pairs
KP_pop_inds = c(
  "Survey", "Survey", "Survey-pair", "Survey-pair", "Survey", "Survey", 
  "Survey", "Survey", "Survey-pair"
)
KP_prob_inds = c("Survey-pair", "Survey")
KP_inds = c("Survey-pair", "Survey", "Survey-pair", "Survey", "Survey")

# Function to plot simulated versus expected numbers of kin-pairs for one type
# of kin-pair
nsKPsPlot = function(i, pop = F, prob = F) {
  if (pop) {
    diffs = ns.kps.pop.wtn.est.errs()[, , i]
    xlab = KP_pop_inds[i]
    main = kp.tps.pop.wtn[i]
  } else {
    diffs = checks.lst()$ns.KPs.lst[[i]] / checks.lst()$exp.ns.KPs.lst[[i]] - 1 
    xlab = KP_inds[i]
    main = KP_names[i]
  }
  if (prob) {
    diffs = 
      t(t(probs.KPs.lst()[[i]]) / exp.probs.KPs.lst()[[i]] - 1)
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

# All pairs within survey years for whole population
output$nsAPsWtnPop = renderPlot(nsKPsPlot(1, T))
# # All pairs between survey years for whole population
# output$nsAPsBtnPop = renderPlot(nsKPsPlot(2, T))
# # Self-pairs between survey years for whole population
# output$nsSPsBtnPop = renderPlot(nsKPsPlot(3, T))
# # Same-mother pairs within survey years for whole population
# output$nsSMPsWtnPop = renderPlot(nsKPsPlot(4, T))
# # Same-father pairs within survey years for whole population
# output$nsSFPsWtnPop = renderPlot(nsKPsPlot(5, T))
# # Same-mother pairs between survey years for whole population
# output$nsSMPsBtnPop = renderPlot(nsKPsPlot(9, T))

# # Self-pair probabilities between survey years for whole population
# output$probSPsBtnPop = renderPlot(nsKPsPlot(1, prob = T))
# # Same-mother pair probabilities within survey years for whole population
# output$probSMPsWtnPop = renderPlot(nsKPsPlot(2, prob = T))

# # Self-pairs between samples
# output$nsSPsBtn = renderPlot(nsKPsPlot(1))
# # Parent-offspring pairs within samples
# output$nsPOPsWtn = renderPlot(nsKPsPlot(2))
# # Parent-offspring pairs between samples
# output$nsPOPsBtn = renderPlot(nsKPsPlot(3))
# # Same-mother pairs within samples
# output$nsSMPsWtn = renderPlot(nsKPsPlot(4))
# # Half-sibling pairs within samples
# output$nsHSPsWtn = renderPlot(nsKPsPlot(5))

### First study simulated
## First life-histories
output$firstLifeHists = renderTable({
  alv_mat = attributes(sim.lst()$hists.lst[[1]])$alv.mat
  mode(alv_mat) = "integer"
  df = data.frame(head(alv_mat))
  names(df) = (f.year() - hist.len() + 1):f.year()
  df
}, rownames = T)

## First sample-histories
output$firstSampHists = renderTable(head(data.frame(sim.lst()$hists.lst[[1]])))

## Numbers of kin-pairs
# Within surveys
output$firstNsKPsPopWtn = renderTable({
  kps.mat = t(checks.lst()$ns.kps.pop.wtn.arr[1, , ])
  mode(kps.mat) = 'integer'
  df = data.frame(kps.mat, row.names = kp.tps.pop.wtn)
  names(df) = srvy.yrs()
  df
}, rownames = T)

# Between surveys
output$firstNsKPsPopBtn = renderTable({
  kps.mat = t(checks.lst()$ns.kps.pop.btn.arr[1, , ])
  mode(kps.mat) = 'integer'
  df = data.frame(kps.mat, row.names = kp.tps.pop.btn)
  names(df) = apply(combn(srvy.yrs(), 2), 2, paste, collapse = "-")
  df
}, rownames = T)