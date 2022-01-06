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
  
  # Create matrices for population trajectories (as columns), numbers of births
  # and deaths, numbers captured, and expected and observed numbers of kin pairs
  N.t.mat <- ns.mtr.mat <- matrix(nrow = n_sims(), ncol = hist.len())
  ns_bths <- ns_dths <- matrix(nrow = n_sims(), ncol = hist.len() - 1)
  ns.caps.mat <- ns.clvng.caps.mat <- ns.clvng.mat <- ns.POPs.wtn.mat <- 
    ns.HSPs.wtn.mat <- exp.ns.HSPs.wtn.mat <- exp.ns.POPs.wtn.mat <- 
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
      ns.POPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.btn
      ns.SPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.SPs.btn
      
      # Find expected numbers of kin pairs
      exp.ns.kps.lst <- FindExpNsKPs(
        k(), n.srvy.prs(), exp.N.fin, lambda(), f.year(), srvy.yrs(), phi(), 
        rho(), ns.caps, alpha()
      )
      
      # Record in matrices
      exp.ns.POPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.wtn
      exp.ns.HSPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.HSPs.wtn
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
    ns.POPs.wtn.mat = ns.POPs.wtn.mat,
    ns.HSPs.wtn.mat = ns.HSPs.wtn.mat,
    ns.POPs.btn.mat = ns.POPs.btn.mat,
    ns.SPs.btn.mat = ns.SPs.btn.mat,
    exp.ns.POPs.wtn.mat = exp.ns.POPs.wtn.mat,
    exp.ns.HSPs.wtn.mat = exp.ns.HSPs.wtn.mat,
    exp.ns.POPs.btn.mat = exp.ns.POPs.btn.mat,
    exp.ns.SPs.btn.mat = exp.ns.SPs.btn.mat,
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
  df = data.frame(br, gr, sr)
  names(df) = c("Birth rate", "Growth rate", "Survival rate")
  df
}, digits = 3)

# Effect of dependency between births and numbers mature on probability of
# breeding
output$bthsNMtr = renderTable({
  mean_bths = colMeans(checks.lst()$ns_bths)
  mean.ns.mtr = colMeans(checks.lst()$ns.mtr.mat[, -hist.len()])
  ratio.mns = mean_bths / mean.ns.mtr
  mn.ratio = 
    colMeans(checks.lst()$ns_bths / checks.lst()$ns.mtr.mat[, -hist.len()])
  df = rbind(mean_bths, mean.ns.mtr, mn.ratio, ratio.mns)
  rownames(df) = c(
    "Average numbers born", "Average numbers mature", "Average ratio",
    "Ratio of averages"
  )
  colnames(df) = (f.year() - hist.len() + 2):f.year()
  df
}, rownames = T, digits = 3)

# First life histories from first study
output$alive = renderTable({
  df = data.frame(head(attributes(sim.lst()$hists.lst[[1]])$alv.mat))
  names(df) = (f.year() - hist.len() + 1):f.year()
  df
}, rownames = T)

# Function to plot simulated versus expected numbers of kin-pairs
nsKPsPlot = function(n_obs, n_exp, x, xlab, kp_type) {
  diffs = n_obs - n_exp
  colnames(diffs) = x
  boxplot(
    diffs, main = kp_type, xlab = xlab, ylab = "Observed - expected numbers"
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

# Parent-offspring pairs within samples
output$nsPOPsWtn <- renderPlot({
  nsKPsPlot(
    checks.lst()$ns.POPs.wtn.mat, checks.lst()$exp.ns.POPs.wtn.mat, 
    srvy.yrs(), "Survey", "Parent-offspring pairs within samples"
  )
})

# Half-sibling pairs within samples
output$nsHSPsWtn <- renderPlot({
  nsKPsPlot(
    checks.lst()$ns.HSPs.wtn.mat, checks.lst()$exp.ns.HSPs.wtn.mat, 
    srvy.yrs(), "Survey", "Half-sibling pairs within samples"
  )
})

# Parent-offspring pairs between samples
output$nsPOPsBtn <- renderPlot({
  nsKPsPlot(
    checks.lst()$ns.POPs.btn.mat, checks.lst()$exp.ns.POPs.btn.mat, 
    apply(combn(srvy.yrs(), 2), 2, paste, collapse = "-"), "Survey-pair", 
    "Parent-offspring pairs between samples"
  )
})

# Self-pairs between samples
output$nsSPsBtn <- renderPlot({
  nsKPsPlot(
    checks.lst()$ns.SPs.btn.mat, checks.lst()$exp.ns.SPs.btn.mat, 
    apply(combn(srvy.yrs(), 2), 2, paste, collapse = "-"), "Survey-pair", 
    "Self-pairs between samples"
  )
})

# Display percentage of animals captured for which the parents are unknown
output$percUnknPrnts <- renderText({
  # nPOPs = checks.lst()$ns.POPs.wtn.mat
  # 
  # # Pairs are lost quadratically with animals
  # prpns.pairs.lost = 1 - (1 - checks.lst()$prpn.prnts.unkn.vec)^2
  # 
  # paste(
  #   "Expected number of parent-offspring pairs within samples
  #   lost due to unknown parents:",
  #   signif(mean(prpns.pairs.lost * checks.lst()$ns.POPs.wtn.mat), 3),
  #   "\n"
  # )
  paste0(
    "Percentage of captured animals with unknown parents (1DP): ", 
    round(mean(checks.lst()$prpn.prnts.unkn.vec) * 100, 1), "%"
  )
})

# Print head of first study
output$dataHead <- renderTable(head(data.frame(sim.lst()$hists.lst[[1]])))