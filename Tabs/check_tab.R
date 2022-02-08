# Display parameter values
output$checkParVals <- renderTable({
  par_vals_df(sim.par.vals(), sim.par.names())
}, digits = 3)
# Display simulation values
output$checkSimVals = renderTable(sim.vals())

# Calculate checks for simulated studies
checks.lst = reactive({
  # Create matrices for population trajectories (as columns), numbers of mature
  # animals, numbers of animals of ages up to the history length that survived
  # to the final year, numbers of same-mother pairs with one born in final year,
  # numbers of births and deaths, numbers captured, and expected and observed
  # numbers of kin pairs
  N.t.mat <- matrix(nrow = n_sims(), ncol = hist.len())
  ns.SMPs.t.f.yr.pop.mat <- ns.SFPs.t.f.yr.pop.mat <- ns.SFPs.t.t.pop.mat <-
    matrix(nrow = n_sims(), ncol = hist.len() - 1)
  ns.caps.mat <- ns.clvng.caps.mat <- ns.clvng.mat <- ns.APs.wtn.pop.mat <-
    ns.SMPs.wtn.pop.mat <- ns.SFPs.wtn.pop.mat <- ns.FSPs.wtn.pop.mat <-
    ns.HSPs.wtn.pop.mat <-
    ns.POPs.wtn.mat <- ns.SMPs.wtn.mat <- ns.HSPs.wtn.mat <- 
    exp.ns.HSPs.wtn.mat <- exp.ns.POPs.wtn.mat <- exp.ns.SMPs.wtn.mat <-
    matrix(nrow = n_sims(), ncol = k())
  ns.APs.btn.pop.mat <- ns.SPs.btn.pop.mat <- ns.SMPs.btn.pop.mat <-
    ns.SFPs.btn.pop.mat <- ns.POPs.btn.mat <- 
    ns.SPs.btn.mat <- exp.ns.POPs.btn.mat <- exp.ns.SPs.btn.mat <- 
    matrix(nrow = n_sims(), ncol = n.srvy.prs())
  
  # Create vectors for proportions with unknown parents, super-population
  # sizes, and number of same-mother pairs in final year
  prpn.prnts.unkn.vec <- Ns.vec <- numeric(n_sims())
  
  # Display progress
  cat("Checking study: ")
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n_sims()) {
      # Display progress
      if (hist.ind %% 100 == 1) cat(hist.ind, "")
      
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hist.ind]]
      
      # Record population curve
      N.t.mat[hist.ind, ] <- attributes(pop.cap.hist)$N.t.vec
      
      # Record super-population size
      Ns.vec[hist.ind] <- attributes(pop.cap.hist)$Ns
      
      # Get numbers captured and calving in each survey
      ns.caps <- attributes(pop.cap.hist)$ns.caps
      ns.caps.mat[hist.ind, ] <- ns.caps
      ns.clvng.mat[hist.ind, ] <- attributes(pop.cap.hist)$ns.clvng
      ns.clvng.caps.mat[hist.ind, ] <- colSums(
        pop.cap.hist[, 4:(3 + k())] * pop.cap.hist[, (4 + k()):(3 + 2 * k())]
      )
      
      # Find proportion captured with unknown parents
      prpn.prnts.unkn.vec[hist.ind] <- mean(is.na(pop.cap.hist$mum))
      
      # Numbers of kin-pairs in whole population
      ns.kps.pop.lst = FindNsKinPairsPop(
        attributes(pop.cap.hist)$N.t.vec[s.yr.inds()], 
        attributes(pop.cap.hist)$alv.mat[, s.yr.inds()], 
        attributes(pop.cap.hist)$ID, 
        attributes(pop.cap.hist)$mum, 
        attributes(pop.cap.hist)$dad, 
        k()
      )
      ns.APs.wtn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.APs.wtn.pop
      ns.APs.btn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.APs.btn.pop
      ns.SPs.btn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.SPs.btn.pop
      ns.SMPs.wtn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.SMPs.wtn.pop
      ns.SFPs.wtn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.SFPs.wtn.pop
      ns.FSPs.wtn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.FSPs.wtn.pop
      ns.HSPs.wtn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.HSPs.wtn.pop
      ns.SMPs.btn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.SMPs.btn.pop
      ns.SFPs.btn.pop.mat[hist.ind, ] = ns.kps.pop.lst$ns.SFPs.btn.pop
      
      
      # Find numbers of known kin pairs
      ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
      ns.POPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.wtn
      ns.HSPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.HSPs.wtn
      ns.SMPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.SMPs.wtn
      ns.POPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.btn
      ns.SPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.SPs.btn
      
      # Find expected numbers of kin pairs
      exp.ns.kps.lst <- FindExpNsKPs(
        k(), n.srvy.prs(), exp.N.fin(), lambda(), f.year(), srvy.yrs(), phi(), 
        rho(), ns.caps, alpha()
      )
      exp.ns.POPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.wtn
      exp.ns.HSPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.HSPs.wtn
      exp.ns.SMPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.SMPs.wtn
      exp.ns.POPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.btn
      exp.ns.SPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.SPs.btn
      
      # Same-mother and father pairs whole population with one born in final
      # year
      age.0 = attributes(pop.cap.hist)$f.age == 0
      mums.of.age.0 = attributes(pop.cap.hist)$mum[age.0]
      dads.of.age.0 = attributes(pop.cap.hist)$dad[age.0]
      alv.f.yr = attributes(pop.cap.hist)$alv.mat[, hist.len()] == 1
      for (t in (hist.len() - 2):1) {
        age.t = attributes(pop.cap.hist)$f.age == t & alv.f.yr
        mums.of.age.t = attributes(pop.cap.hist)$mum[age.t]
        dads.of.age.t = attributes(pop.cap.hist)$dad[age.t]
        ns.SMPs.t.f.yr.pop.mat[hist.ind, hist.len() - 1 - t] = 
          sum(mums.of.age.t %in% mums.of.age.0)
        
        # There may be more than one animal born in the final year with the same
        # father, and there is one same-father pair for each of them, for each
        # animal born in the current year with that father.  This code
        # multiplies the corresponding counts for each father and takes the sum
        # (computes the dot product).
        max.dad.id = max(dads.of.age.t, dads.of.age.0)
        ns.SFPs.t.f.yr.pop.mat[hist.ind, hist.len() - 1 - t] = 
          tabulate(dads.of.age.t, max.dad.id) %*% 
          tabulate(dads.of.age.0, max.dad.id)
        
        # Same-father pairs born in the current year
        ns.SFPs.t.t.pop.mat[hist.ind, hist.len() - 1 - t] =
          sum(choose(tabulate(dads.of.age.t), 2))
      }
      
      incProgress(1/n_sims())
    }
  }, value = 0, message = "Checking simulations")
  
  list(
    N.t.mat = N.t.mat, 
    prpn.prnts.unkn.vec = prpn.prnts.unkn.vec,
    ns.caps.mat = ns.caps.mat,
    ns.KPs.pop.lst = list(
      N.wtn.pop.mat = N.t.mat[, s.yr.inds()],
      ns.APs.wtn.pop.mat = ns.APs.wtn.pop.mat,
      ns.APs.btn.pop.mat = ns.APs.btn.pop.mat,
      ns.SPs.btn.pop.mat = ns.SPs.btn.pop.mat,
      ns.SMPs.wtn.pop.mat = ns.SMPs.wtn.pop.mat,
      ns.SFPs.wtn.pop.mat = ns.SFPs.wtn.pop.mat,
      ns.FSPs.wtn.pop.mat = ns.FSPs.wtn.pop.mat,
      ns.HSPs.wtn.pop.mat = ns.HSPs.wtn.pop.mat,
      ns.SMPs.btn.pop.mat = ns.SMPs.btn.pop.mat
    ),
    ns.KPs.lst = list(
      ns.SPs.btn.mat = ns.SPs.btn.mat,
      ns.POPs.wtn.mat = ns.POPs.wtn.mat,
      ns.POPs.btn.mat = ns.POPs.btn.mat,
      ns.SMPs.wtn.mat = ns.SMPs.wtn.mat,
      ns.HSPs.wtn.mat = ns.HSPs.wtn.mat
    ),
    exp.ns.KPs.lst = list(
      exp.ns.SPs.btn.mat = exp.ns.SPs.btn.mat,
      exp.ns.POPs.wtn.mat = exp.ns.POPs.wtn.mat,
      exp.ns.POPs.btn.mat = exp.ns.POPs.btn.mat,
      exp.ns.SMPs.wtn.mat = exp.ns.SMPs.wtn.mat,
      exp.ns.HSPs.wtn.mat = exp.ns.HSPs.wtn.mat
    ),
    ns.SMPs.t.f.yr.pop.mat = ns.SMPs.t.f.yr.pop.mat,
    ns.SFPs.t.f.yr.pop.mat = ns.SFPs.t.f.yr.pop.mat,
    ns.SFPs.t.t.pop.mat = ns.SFPs.t.t.pop.mat
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
output$obsParVals = renderTable({
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

# Names of types of kin-pairs
KP_pop_names = c(
  "Population sizes in survey years",
  "All-pairs within surveys", "All-pairs between surveys",
  "Self-pairs between surveys", "Same-mother pairs within surveys",
  "Same-father pairs within surveys", "Full-sibling pairs within surveys",
  "Half-sibling pairs within surveys", "Same-mother pairs between surveys"
)
KP_prob_names = c(
  "Self-pair probabilities between surveys", 
  "Same-mother pair probabilities within surveys"
)
KP_names = c(
  "Self-pairs between surveys", "Parent-offspring pairs within surveys", 
  "Parent-offspring pairs between surveys", "Same-mother pairs within surveys",
  "Half-sibling pairs within surveys"
)
# Number of types of kin-pairs
n_KP_pop_types = length(KP_pop_names)
n_KP_prob_types = length(KP_prob_names)
n_KP_types = length(KP_names)
# Indices for types of kin-pairs
KP_pop_inds = c(
  "Survey", "Survey", "Survey-pair", "Survey-pair", "Survey", "Survey", 
  "Survey", "Survey", "Survey-pair"
)
KP_prob_inds = c("Survey-pair", "Survey")
KP_inds = c("Survey-pair", "Survey", "Survey-pair", "Survey", "Survey")

# Expected numbers of kin-pairs for whole population
exp.ns.KPs.pop.lst = reactive({
  # All pairs within surveys and between pairs of surveys
  exp.N.srvy.yrs = exp.N.fin() / lambda()^(f.year() - srvy.yrs())
  exp.ns.APs.wtn = choose(exp.N.srvy.yrs, 2)
  exp.ns.APs.btn = as.vector(combn(exp.N.srvy.yrs, 2, function(x) x[1] * x[2]))
  
  # Self-pairs between pairs of surveys
  exp.ns.SPs.btn = as.vector(combn(srvy.yrs(), 2, function(srvy.pr) {
    exp.N.fin() / lambda()^(f.year() - srvy.pr[1]) *
      phi()^(srvy.pr[2] - srvy.pr[1])
  }))
  
  # Same-mother pairs within surveys
  exp.ns.SMPs.wtn = 2 * exp.N.srvy.yrs * (1 - phi() / lambda())^2 *
    (lambda() / phi())^alpha() * lambda() * phi()^2 / (lambda() - phi()^2)^2
  
  # Same-father pairs within surveys
  exp.ns.SFPs.wtn = (1 - phi() / lambda())^2 *
    (lambda() / phi())^alpha() * lambda() * 
    (exp.N.srvy.yrs * 
       (2 * phi()^3 / (lambda() - phi()^2)^2 + 
          lambda() / (lambda() - phi()^2)) -
       lambda()^alpha() / (1 - phi()^2))
  
  # Full-sibling pairs within surveys
  exp.ns.FSPs.wtn = 4 *  (1 - phi() / lambda())^2 *
    (lambda() / phi())^(2 * alpha()) * 
    phi()^4 / ((lambda() - phi()^3) * (1 - phi()^2))
  
  # Half-sibling pairs within surveys
  exp.ns.HSPs.wtn = exp.ns.SMPs.wtn + exp.ns.SFPs.wtn - 2 * exp.ns.FSPs.wtn
  
  # Same-mother pairs between surveys
  exp.ns.SMPs.btn = as.vector(combn(1:k(), 2, function(s.inds) {
    exp.ns.SMPs.wtn[s.inds[1]] * phi()^(s.inds[2] - s.inds[1]) *
      ((lambda() - phi()^2) / phi()^2 * (s.inds[2] - s.inds[1]) + 2)
  }))
  
  # Return as list
  list(
    exp.N.srvy.yrs = exp.N.srvy.yrs,
    exp.ns.APs.wtn = exp.ns.APs.wtn,
    exp.ns.APs.btn = exp.ns.APs.btn,
    exp.ns.SPs.btn = exp.ns.SPs.btn,
    exp.ns.SMPs.wtn = exp.ns.SMPs.wtn,
    exp.ns.SFPs.wtn = exp.ns.SFPs.wtn,
    exp.ns.FSPs.wtn = exp.ns.FSPs.wtn,
    exp.ns.HSPs.wtn = exp.ns.HSPs.wtn,
    exp.ns.SMPs.btn = exp.ns.SMPs.btn
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

# Table average percentage differences for whole population
output$nsKPsPop = renderTable({
  KPstab(
    n_KP_pop_types, checks.lst()$ns.KPs.pop.lst, 
    exp.ns.KPs.pop.lst(), KP_pop_names, T
  )
})

# Table average percentage differences for whole population
output$nsKPsProb = renderTable({
  KPstab(
    n_KP_prob_types, KPs.rate.lst(), 
    KPs.prob.lst(), KP_prob_names, T
  )
})

# Table average percentage differences for captured animals
output$nsKPsCap = renderTable({
  KPstab(
    n_KP_types, checks.lst()$ns.KPs.lst, 
    checks.lst()$exp.ns.KPs.lst, KP_names, F
  )
})

# Same-mother pairs in final year between animals born in final year and each
# other year in population history
output$SMPsFYear = renderTable({
  mean.ns.SMPs.f.yr = colMeans(checks.lst()$ns.SMPs.t.f.yr.pop.mat)
  mean.ns.SFPs.f.yr = colMeans(checks.lst()$ns.SFPs.t.f.yr.pop.mat)
  mean.ns.SFPs.t.t = colMeans(checks.lst()$ns.SFPs.t.t.pop.mat)
  
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
  
  df = rbind(
    mean.ns.SMPs.f.yr, c(exp.ns.SMPs.f.yr, NA), 
    mean.ns.SFPs.f.yr, c(exp.ns.SFPs.f.yr, NA),
    mean.ns.SFPs.t.t, c(exp.ns.SFPs.t.t, NA)
  )
  rownames(df) = c(
    "Avg(SMP{t,f.yr,f.yr})", "E(SMP{t,f.yr,f.yr})",
    "Avg(SFP{t,f.yr,f.yr})", "E(SFP{t,f.yr,f.yr})",
    "Avg(SFP{t,t,f.yr})", "E(SFP{t,t,f.yr})"
  )
  colnames(df) = (f.year() - hist.len() + 2):f.year()
  df
}, rownames = T, digits = 3)

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

# All pairs within survey years for whole population
output$nsAPsWtnPop <- renderPlot(nsKPsPlot(1, T))
# All pairs between survey years for whole population
output$nsAPsBtnPop <- renderPlot(nsKPsPlot(2, T))
# Self-pairs between survey years for whole population
output$nsSPsBtnPop <- renderPlot(nsKPsPlot(3, T))
# Same-mother pairs within survey years for whole population
output$nsSMPsWtnPop <- renderPlot(nsKPsPlot(4, T))
# Same-father pairs within survey years for whole population
output$nsSFPsWtnPop <- renderPlot(nsKPsPlot(5, T))
# Same-mother pairs between survey years for whole population
output$nsSMPsBtnPop <- renderPlot(nsKPsPlot(9, T))

# Self-pair probabilities between survey years for whole population
output$probSPsBtnPop <- renderPlot(nsKPsPlot(1, prob = T))
# Same-mother pair probabilities within survey years for whole population
output$probSMPsWtnPop <- renderPlot(nsKPsPlot(2, prob = T))

# Self-pairs between samples
output$nsSPsBtn <- renderPlot(nsKPsPlot(1))
# Parent-offspring pairs within samples
output$nsPOPsWtn <- renderPlot(nsKPsPlot(2))
# Parent-offspring pairs between samples
output$nsPOPsBtn <- renderPlot(nsKPsPlot(3))
# Same-mother pairs within samples
output$nsSMPsWtn <- renderPlot(nsKPsPlot(4))
# Half-sibling pairs within samples
output$nsHSPsWtn <- renderPlot(nsKPsPlot(5))

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