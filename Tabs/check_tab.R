# Display parameter values
output$checkParVals <- renderTable({
  par_vals_df(sim.par.vals(), sim.par.names())
}, digits = 3)
# Display simulation values
output$checkSimVals = renderTable(sim.vals())

# Check simulated studies
checks.lst = reactive({
  # Expected final population size
  exp.N.fin <- exp.N.t()[hist.len()]
  
  # Create matrices for population trajectories (as columns), numbers of mature
  # animals, numbers of animals of ages up to the history length that survived
  # to the final year, numbers of same-mother pairs with one born in final year,
  # numbers of births and deaths, numbers captured, and expected and observed
  # numbers of kin pairs
  N.t.mat <- ns.mtr.mat <- matrix(nrow = n_sims(), ncol = hist.len())
  ns_bths <- ns_dths <- ns.srvvd.f.yr.mat <- ns.SMPs.f.yr.mat <-
    matrix(nrow = n_sims(), ncol = hist.len() - 1)
  ns.caps.mat <- ns.clvng.caps.mat <- ns.clvng.mat <- ns.POPs.wtn.mat <- 
    ns.SMPs.wtn.mat <- ns.HSPs.wtn.mat <- 
    exp.ns.HSPs.wtn.mat <- exp.ns.POPs.wtn.mat <- 
    exp.ns.SMPs.wtn.mat <-
    matrix(nrow = n_sims(), ncol = k())
  ns.POPs.btn.mat <- ns.SPs.btn.mat <- exp.ns.POPs.btn.mat <- 
    exp.ns.SPs.btn.mat <- matrix(nrow = n_sims(), ncol = n.srvy.prs())
  
  # Create vectors for proportions with unknown parents, and super-population
  # sizes
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
      
      # Find birth and survival rates
      f_alv_mat = attributes(pop.cap.hist)$alv.mat
      bs_n_ds = f_alv_mat[, -1] != f_alv_mat[, -hist.len()]
      ns_bths[hist.ind, ] = colSums(bs_n_ds & f_alv_mat[, -hist.len()] == 0)
      ns_dths[hist.ind, ] = colSums(bs_n_ds & f_alv_mat[, -hist.len()] == 1)
      
      # Record numbers of mature animals
      ns.mtr.mat[hist.ind, ] <- attributes(pop.cap.hist)$ns.mtr
      
      # Record number animals of ages up to the history length that survived to
      # the final year
      ns.srvvd.f.yr.mat[hist.ind, ] <- table(
        factor(
          attributes(pop.cap.hist)$f.age[f_alv_mat[, hist.len()] == 1],
          levels = (hist.len() - 2):0
        )
      )      
      
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
      
      # Find numbers of known kin pairs
      ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), pop.cap.hist)
      
      # Record in matrices
      ns.POPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.wtn
      ns.HSPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.HSPs.wtn
      ns.SMPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.SMPs.wtn
      ns.POPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.btn
      ns.SPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.SPs.btn
      
      # Same-mother pairs with one born in final year
      age.0 = attributes(pop.cap.hist)$f.age == 0
      mums.of.age.0 = attributes(pop.cap.hist)$mum[age.0]
      alv.f.yr = f_alv_mat[, hist.len()] == 1
      for (t in (hist.len() - 2):1) {
        age.t = attributes(pop.cap.hist)$f.age == t & alv.f.yr
        mums.of.age.t = attributes(pop.cap.hist)$mum[age.t]
        ns.SMPs.f.yr.mat[hist.ind, hist.len() - 1 - t] = 
          sum(mums.of.age.t %in% mums.of.age.0)
      }

      # Find expected numbers of kin pairs
      exp.ns.kps.lst <- FindExpNsKPs(
        k(), n.srvy.prs(), exp.N.fin, lambda(), f.year(), srvy.yrs(), phi(), 
        rho(), ns.caps, alpha()
      )
      
      # Record in matrices
      exp.ns.POPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.wtn
      exp.ns.HSPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.HSPs.wtn
      exp.ns.SMPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.SMPs.wtn
      exp.ns.POPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.btn
      exp.ns.SPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.SPs.btn
      
      incProgress(1/n_sims())
    }
  }, value = 0, message = "Checking simulations")
  
  list(
    N.t.mat = N.t.mat, 
    mean.N.t = colMeans(N.t.mat),
    ns_bths = ns_bths,
    ns_dths = ns_dths,
    ns.mtr.mat = ns.mtr.mat,
    ns.srvvd.f.yr.mat = ns.srvvd.f.yr.mat,
    ns.KPs.lst = list(
      ns.SPs.btn.mat = ns.SPs.btn.mat,
      ns.POPs.wtn.mat = ns.POPs.wtn.mat,
      ns.POPs.btn.mat = ns.POPs.btn.mat,
      ns.SMPs.wtn.mat = ns.SMPs.wtn.mat,
      ns.HSPs.wtn.mat = ns.HSPs.wtn.mat
    ),
    ns.SMPs.f.yr.mat = ns.SMPs.f.yr.mat,
    exp.ns.KPs.lst = list(
      exp.ns.SPs.btn.mat = exp.ns.SPs.btn.mat,
      exp.ns.POPs.wtn.mat = exp.ns.POPs.wtn.mat,
      exp.ns.POPs.btn.mat = exp.ns.POPs.btn.mat,
      exp.ns.SMPs.wtn.mat = exp.ns.SMPs.wtn.mat,
      exp.ns.HSPs.wtn.mat = exp.ns.HSPs.wtn.mat
    ),
    prpn.prnts.unkn.vec = prpn.prnts.unkn.vec,
    ns.caps.mat = ns.caps.mat
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
    (f.year() - hist.len() + 1):f.year(), checks.lst()$mean.N.t, 
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

# Observed parameter values
output$obsParVals = renderTable({
  br = mean(checks.lst()$ns_bths / checks.lst()$N.t.mat[, -hist.len()])
  sr = 1 - mean(checks.lst()$ns_dths / checks.lst()$N.t.mat[, -hist.len()])
  gr = mean(checks.lst()$mean.N.t[-1] / checks.lst()$mean.N.t[-hist.len()])
  brmf = 2 * mean(checks.lst()$ns_bths / checks.lst()$ns.mtr.mat[, -1])
  # Mature females the year before
  # exp_brmf = 2 * (lambda() / phi() - 1) * (lambda() / phi())^alpha()
  exp_brmf = 2 * (1 - phi() / lambda()) * (lambda() / phi())^alpha()
  df = data.frame(br, gr, sr, brmf, exp_brmf)
  names(df) = c(
    "Birth rate", "Growth rate", "Survival rate", 
    "Birth rate among mature females", 
    "Expected birth rate among mature females"
  )
  df
}, digits = 3)

# Checks for numbers of kin-pairs

# Names of types of kin-pairs
KP_names = c(
  "Self-pairs between surveys", "Parent-offspring pairs within surveys", 
  "Parent-offspring pairs between surveys", "Same-mother pairs within surveys",
  "Half-sibling pairs within surveys"
)
# Number of types of kin-pairs
n_KP_types = length(KP_names)
# Indices for types of kin-pairs
KP_inds = c("Survey-pair", "Survey", "Survey-pair", "Survey", "Survey")

# Show average percentage differences between simulated and expected numbers
# of kin-pairs
output$nsKPs = renderTable({
  diffs = numeric(n_KP_types)
  for (i in 1:n_KP_types) {
    diffs[i] = perc(mean(
      (checks.lst()$ns.KPs.lst[[i]] - checks.lst()$exp.ns.KPs.lst[[i]]) /
      checks.lst()$exp.ns.KPs.lst[[i]]
    ))
  }
  df = data.frame(matrix(diffs, nrow = 1))
  names(df) = KP_names
  df
})

# Show percentage of animals captured for which the parents are unknown
output$percUnknPrnts <- renderText({
  paste0(
    "Percentage of captured animals with unknown parents (1DP): ", 
    round(mean(checks.lst()$prpn.prnts.unkn.vec) * 100, 1), "%"
  )
})

# Function to plot simulated versus expected numbers of kin-pairs for one type
# of kin-pair
nsKPsPlot = function(i) {
  diffs = checks.lst()$ns.KPs.lst[[i]] - checks.lst()$exp.ns.KPs.lst[[i]]
  xlab = KP_inds[i]
  if (xlab == "Survey") colnames(diffs) = srvy.yrs()
  else colnames(diffs) = apply(combn(srvy.yrs(), 2), 2, paste, collapse = "-")
  {
    boxplot(
      diffs, main = KP_names[i], xlab = xlab, 
      ylab = "Observed - expected numbers"
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

# Apply function to each type of kin-pair

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
  df = data.frame(head(attributes(sim.lst()$hists.lst[[1]])$alv.mat))
  names(df) = (f.year() - hist.len() + 1):f.year()
  df
}, rownames = T)

# Effect of dependency between births and numbers mature on probability of
# breeding
output$bthsNMtr = renderTable({
  mean_bths = colMeans(checks.lst()$ns_bths)
  exp_bths = exp.N.t()[-1] * (1 - phi() / lambda())
  mean.ns.mtr = colMeans(checks.lst()$ns.mtr.mat[, -1])
  exp.ns.mtr = exp.N.t()[-1] * (phi() / lambda())^alpha()
  mean.ns.srvvd.f.yr = colMeans(checks.lst()$ns.srvvd.f.yr.mat)
  exp.ns.srvvd.f.yr = (1 - phi() / lambda()) * 
    (phi() / lambda())^((hist.len() - 2):0) * exp.N.t()[hist.len()]
  mean.ns.SMPs.f.yr = colMeans(checks.lst()$ns.SMPs.f.yr.mat)
  exp.ns.SMPs.f.yr = 2 * exp.N.t()[hist.len()] * (1 - phi() / lambda())^2 *
    (lambda() / phi())^alpha() * (phi()^2 / lambda())^((hist.len() - 2):1)
  
  df = rbind(
    mean_bths, exp_bths, mean.ns.mtr, exp.ns.mtr, mean.ns.srvvd.f.yr,
    exp.ns.srvvd.f.yr, mean.ns.SMPs.f.yr, c(exp.ns.SMPs.f.yr, NA)
  )
  rownames(df) = c(
    "Average numbers born", "Expected numbers born", "Average numbers mature",
    "Expected numbers mature", 
    "Average numbers born and survived to final year",
    "Expected numbers born and survived to final year",
    "Average numbers same-mother pairs with one born in final year",
    "Expected numbers same-mother pairs with one born in final year"
  )
  colnames(df) = (f.year() - hist.len() + 2):f.year()
  df
}, rownames = T, digits = 1)

# Print head of first study
output$dataHead <- renderTable(head(data.frame(sim.lst()$hists.lst[[1]])))