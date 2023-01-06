# Find allele frequencies from genotypes in 2 x L x n_individuals arrays,
# representing two binary SNPs at each locus for each individual, excluding
# repeated samples of the same individual in different surveys.  Frequencies are
# returned as 2 x L matrices representing the frequencies of 0 and 1-coded SNP
# alleles at each locus
frst.ale.frqs = reactive(FindAleFrqs(FS.atts()$unq.smp.gts))

# Probabilities of possible genopairs at each locus as 3 x 3 x L arrays, where
# rows represent the first genotypes, and columns the second, ordered as 00, 01,
# and 11, for binary SNPs.
frst.pss.gp.prbs.KPs = reactive(FindPssGPPsKPs(frst.ale.frqs(), L()))

# Sample histories from first study, (n_animals x n_surveys)
frst.smp.hsts = reactive(as.matrix(fst.std()[, 4:(3 + k())]))

# Indices of samples from first study, (n_samples)
frst.smp.inds = reactive(row(frst.smp.hsts())[as.logical(frst.smp.hsts())])

# Sample genotypes, (2 x L x n_samples), indexed from unique sampled genotypes,
# (2 x L x n_animals), representing two binary SNPs at each locus for
# each individual sampled at least once. They are expanded to arrays by indexing
# sample genotypes.
frst.smp.gts = reactive(FS.atts()$unq.smp.gts[, , frst.smp.inds()])

# Indices of survey-years for each sample, starting from zero, as used in TMB
# C++ objective function
frst.smp.yr.inds = reactive({
  col(frst.smp.hsts())[as.logical(frst.smp.hsts())] - 1
})

# Sample index pairs, 2 x n_pairs matrix of indices of samples in each pair to
# include in likelihood, possibly all pairs or just consecutive pairs
frst.SIPs.fll = reactive(combn(length(frst.smp.inds()), 2))
frst.SIPs.offst = reactive(FindSIPsOffset(k(), frst.smp.yr.inds()))

# Sample-year index pairs
FindSYIPs = reactive(function(smp.ind.prs) {
  # n_pairs x 2 matrix of survey-years for each sample in each pair
  matrix(
    frst.smp.yr.inds()[as.vector(t(smp.ind.prs))], ncol = 2,
    dimnames = list(pair = NULL, sample_year_index = 1:2)
  )
})

# Indices of survey-years for each sample in each pair, starting at zero for
# C++ template, and ordered by survey-year of first sample
frst.SYIPs.fll = reactive(FindSYIPs()(frst.SIPs.fll()))
frst.SYIPs.offst = reactive(FindSYIPs()(frst.SIPs.offst()))

# Nullify genopair log-probabilities when new datasets simulated
observeEvent(input$simulate, frst.LGPPs.KP.fll(NULL))

# Find genopair log-probabilities as n.pairs x n.kp.tps matrix with rows for
# pairs of individuals, and columns for types of kinships considered
observeEvent({
    input$nav.tab
    input$check.sub.tabs
}, {
  if (
    input$nav.tab == "check.tab" && 
    input$check.sub.tabs %in% c("frst.gts.tb", "frst.ests.tb") &&
    is.null(frst.LGPPs.KP.fll())
  ) {
    frst.LGPPs.KP.fll(FindLogGPProbsKP(
      frst.pss.gp.prbs.KPs(), frst.smp.gts(), frst.SIPs.fll(), L()
    ))
  }
})
frst.LGPPs.KP.offst = reactive({
  FindLogGPProbsKP(
    frst.pss.gp.prbs.KPs(), frst.smp.gts(), frst.SIPs.offst(), L()
  )
})

# Expected values of HSP vs UP PLODs given kinships
exp.plod.KP = reactive({
  # Possible values of HSP vs UP PLODs, leaving division by number of loci to
  # next step after summation
  pss.plods = 
    log(frst.pss.gp.prbs.KPs()[, , , 2] / frst.pss.gp.prbs.KPs()[, , , 1])
  
  # Unrelated, parent-offspring, and self-pairs
  prbs.plods = frst.pss.gp.prbs.KPs() * rep(pss.plods, 4)
  exp.plod.base = colSums(prbs.plods, dims = 3) / L()
  
  # First-cousin, avuncular, and half-sibling pairs
  exp.plod.extd = (c(7, 3) * exp.plod.base[1] + exp.plod.base[3]) / c(8, 4)
  
  # Combine in order from furthest to closest kinship, add names, and return.
  vec = c(exp.plod.base[1], exp.plod.extd, exp.plod.base[2:4])
  names(vec) = c(
    "Unrelated", "First cousin", "Avuncular", "Half-sibling",
    "Parent-offspring", "Self"
  )
  vec
})

# Find half-sibling vs unrelated pairs PLODs from log genopair probabilities
first.plods = reactive({
  (frst.LGPPs.KP.fll()[, 2] - frst.LGPPs.KP.fll()[, 1]) / L()
})

## Outputs

# Table of genotypes of first few individuals captured (can show kin-pairs
# later)
output$firstGTs = renderTable({
  df = data.frame(cbind(
    rep(fst.std()$ID[frst.smp.inds()[1:3]], each = 2), 
    rep(paste0(c("m", "p"), "aternal"), 3),
    matrix(frst.smp.gts()[, , 1:3], nrow = 6)
  ))
  names(df) = c("ID", "Allele", paste0("L", 1:L()))
  df
})

# Table showing allele frequencies
output$firstAFs = renderTable({
  df = data.frame(frst.ale.frqs())
  # df = data.frame(asNumeric(ale.frqs()))
  names(df) = paste0("L", 1:L())
  row.names(df) = 0:1
  df
}, rownames = T)

# Function to format genopair probabilities for display
frmt.gpps = function(gpps) {
  df = data.frame(gpps)
  # df = data.frame(asNumeric(gpps))
  names(df) = row.names(df) = c("00", "01", "11")
  df
}

# Table showing GPPs given unrelated
output$firstGPPsUP = renderTable({
  frmt.gpps(frst.pss.gp.prbs.KPs()[, , 1, 1])
}, rownames = T)

# Table showing GPPs given half-siblings
output$firstGPPsHSP = renderTable({
  frmt.gpps(frst.pss.gp.prbs.KPs()[, , 1, 2])
}, rownames = T)

# Table showing GPPs given parent-offspring
output$firstGPPsPOP = renderTable({
  frmt.gpps(frst.pss.gp.prbs.KPs()[, , 1, 3])
}, rownames = T)

# Table showing GPPs given self-resample
output$firstGPPsSP = renderTable({
  frmt.gpps(frst.pss.gp.prbs.KPs()[, , 1, 4])
}, rownames = T)

# Table of genopair probabilities of first few pairs captured (can show
# kin-pairs later)
output$firstFewLGPPs = renderTable({
  df = data.frame(cbind(
    matrix(fst.std()[frst.smp.inds()[t(frst.SIPs.fll()[, 1:3])], 1], ncol = 2),
    frst.SYIPs.fll()[1:3, ],
    frst.LGPPs.KP.fll()[1:3, ]
  ))
  df[, 1:4] = as.integer(as.matrix(df[, 1:4]))
  names(df) = c("ID1", "ID2", "Survey index 1", "Survey index 2", gp.prb.KP.tps)
  df
})

# Histograms of genopair log-probabilities given basic kinships
output$firstLGPPs = renderPlot({
  par(mfrow = c(1, 4))
  lapply(1:4, function(i) {
    # Plotting may fail if all log-probabilities negative infinity
    if(any(is.finite(frst.LGPPs.KP.fll()[, i]))) {
      hist(
        frst.LGPPs.KP.fll()[, i], main = gp.prb.KP.tps[i],
        xlab = "Log-probability", br = 50
      )
    } else plot.new()
  })
})

# Histograms of genopair log-probabilities given basic kinships
output$frstGpPs = renderPlot({
  par(mfrow = c(1, 4))
  frst.Gp.Ps = FindGPPs(frst.LGPPs.KP.fll())
  lapply(1:4, function(i) {
    hst.dt = hist(frst.Gp.Ps[, i], br = 50, plot = F)
    hst.dt$counts = log(hst.dt$counts + 1)
    plot(
      hst.dt, main = gp.prb.KP.tps[i], xlab = "Probability", 
      ylab = "log (frequency + 1)"
    )
  })
})

# Plot PLODs
output$firstPLODs = renderPlot({
  # Plot plods
  hist(
    first.plods(), main = "HSP vs UP PLODs for all samples",
    sub = paste(L(), "loci"), xlab = "PLOD", breaks = 200,
  )
  
  # Plot expected values
  abline(v = exp.plod.KP()[1], col = 2, lwd = 2)
  
  # Add legend
  legend(
    "topright", col = 1:7, lty = c(0, rep(1, 6)), cex = 1, lwd = 2,
    legend = c(
      "Expected value given kinship", "Unrelated", "First cousin",
      "Avuncular", "Half-sibling", "Parent-offspring", "Self"
    )
  )
})

# Plot PLODs and zoom in on rare values representing likely close-kin
output$firstPLODsRare = renderPlot({
  # Plot uncommon values
  hst.data = hist(first.plods(), breaks = 200, plot = F) 
  hst.data$counts = log(hst.data$counts + 1)
  plot(
    hst.data, main = "Frequency on log scale",
    xlab = "PLOD", ylab = "Log (frequency + 1)"
  )
  
  # Plot expected values
  abline(v = exp.plod.KP(), col = 2:7, lwd = 2)
})

