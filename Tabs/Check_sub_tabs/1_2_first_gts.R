# Genotypes, 2 x L x n_individuals, indexed from unique sampled genotypes,
# (2 x L x n_animals), representing two binary SNPs at each locus for
# each individual sampled at least once. They are expanded to arrays by indexing
# sample genotypes.
frst.gts = reactive(FS.atts()$unq.smp.gts)

# Allele frequencies for first study, 2 x n_loci, representing 0 and 1-coded SNP
# alleles at each locus
frst.ale.frqs = reactive(FindAleFrqs(frst.gts()))

# Possible genopair probabilities given kinships at each locus for first study,
# 3 x 3 x n_loci x n_kinships, for possible genotypes ordered 00, 01, 11
frst.pgpsgks = reactive(FindPssGPPsKPs(frst.ale.frqs(), L(), knshp.chcs))

# Sample histories from first study, n_animals x n_surveys
frst.smp.hsts = reactive(as.matrix(frst.std()[, 4:(3 + k())]))

# Sample-individual indices from first study, n_samples
frst.siis = reactive(row(frst.smp.hsts())[as.logical(frst.smp.hsts())])

# Sample-year indices from first study, n_samples, starting from zero as used in
# TMB C++ objective function
frst.syis = reactive(col(frst.smp.hsts())[as.logical(frst.smp.hsts())] - 1)

# Offset sample index pairs for first study
frst.osips = reactive(FindSIPsOffset(k(), frst.syis()))

# Function to find sample-individual index pairs, n_pairs x 2
FindSIIPs = function(siis, sips) {
  matrix(siis[as.vector(sips)], ncol = 2)
}

# Full and offset sample-individual index pairs, n_pairs x 2
frst.fsiips = reactive(t(combn(frst.siis(), 2)))
frst.osiips = reactive(FindSIIPs(frst.siis(), frst.osips()))

# Function to find sample-year index pairs, n_pairs x 2
FindSYIPs = function(syis, sips) {
  matrix(syis[as.vector(sips)], ncol = 2)
}

# Indices of survey-years for each sample in each pair, starting at zero for
# C++ template, and ordered by survey-year of first sample
frst.fsyips = reactive(t(combn(frst.syis(), 2)))
frst.osyips = reactive(FindSYIPs(frst.syis(), frst.osips()))

# Nullify genopair log-probabilities when new datasets simulated
observeEvent(input$simulate, frst.fglps(NULL))

# If tab changed in app
observeEvent({
    input$nav.tab
    input$check.sub.tabs
}, {
  # If new tab requires genopair probabilities and they have not been computed
  if (
    input$nav.tab == "check.tab" && 
    input$check.sub.tabs %in% c(
      "frst.gts.tb", "frst.lklhds.tb", "frst.ests.tb"
    ) &&
    is.null(frst.fglps())
  ) {
    # Find full set of genopair log-probabilities for first study, n.pairs x
    # n.kinships
    frst.fglps(FindGLPs(frst.pgpsgks(), frst.gts(), frst.fsiips(), L()))
  }
})

# Offset genopair log-probabilities for first study
frst.oglps = reactive(FindGLPs(frst.pgpsgks(), frst.gts(), frst.osiips(), L()))

# Genopair probabilities for optimization
frst.fgps = reactive(FindGPPs(frst.fglps()))
frst.ogps = reactive(FindGPPs(frst.oglps()))

# Expected values of HSP vs UP PLODs given various kinships
exp.plods = reactive({
  # Possible  genopair HSP vs UP log-probability ratios at each locus,
  # 3 x 3 x n_loci
  pglprs = log(frst.pgpsgks()[, , , 2] / frst.pgpsgks()[, , , 1])
  
  # Expected plods for kinship basis, unrelated, half-sibling, parent-offspring,
  # and self-pairs
  exp.plds.bss = colSums(frst.pgpsgks() * rep(pglprs, 4), dims = 3) / L()
  
  # Expected plods for extended kinships, first-cousin and avuncular pairs
  exp.plds.extd = (c(7, 3) * exp.plds.bss[1] + exp.plds.bss[3]) / c(8, 4)
  
  # Combine in order from furthest to closest kinship, add names, and return
  vec = c(exp.plds.bss[1], exp.plds.extd, exp.plds.bss[2:4])
  names(vec) = c(
    "Unrelated", "First cousin", "Avuncular", "Half-sibling",
    "Parent-offspring", "Self"
  )
  vec
})

# Find half-sibling vs unrelated pairs PLODs from log genopair probabilities
first.plods = reactive({
  (frst.fglps()[, 2] - frst.fglps()[, 1]) / L()
})

## Outputs

# Table of genotypes of first few individuals captured
output$firstGTs = renderTable({
  df = data.frame(cbind(
    rep(frst.std()$ID[frst.siis()[1:3]], each = 2), 
    rep(paste0(c("m", "p"), "aternal"), 3),
    matrix(frst.gts()[, , 1:3], nrow = 6)
  ))
  names(df) = c("ID", "Allele", paste0("L", 1:L()))
  df
})

# Table showing allele frequencies
output$firstAFs = renderTable({
  df = data.frame(frst.ale.frqs())
  names(df) = paste0("L", 1:L())
  row.names(df) = 0:1
  df
}, rownames = T)

# Table showing possible genotype probabilities
output$firstGtPrbs = renderTable({
  df = data.frame(FindPssGtPrbs(frst.ale.frqs()))
  names(df) = paste0("L", 1:L())
  row.names(df) = pss.gt.lbls
  df
}, rownames = T)

# Function to format genopair probabilities for display
FrmtGPs = function(gpps) {
  df = data.frame(gpps)
  names(df) = row.names(df) = pss.gt.lbls
  df
}

# Table showing genopair probabilities given unrelated
output$firstGPsUP = renderTable({
  FrmtGPs(frst.pgpsgks()[, , 1, 1])
}, rownames = T)

# Table showing genopair probabilities given half-siblings
output$firstGPsHSP = renderTable({
  FrmtGPs(frst.pgpsgks()[, , 1, 2])
}, rownames = T)

# Table showing genopair probabilities given parent-offspring
output$firstGPsPOP = renderTable({
  FrmtGPs(frst.pgpsgks()[, , 1, 3])
}, rownames = T)

# Table showing genopair probabilities given self-resample
output$firstGPsSP = renderTable({
  FrmtGPs(frst.pgpsgks()[, , 1, 4])
}, rownames = T)

# Table of genopair log-probabilities of first few pairs captured
output$firstFewLGPPs = renderTable({
  df = data.frame(cbind(
    matrix(frst.std()$ID[as.vector(frst.fsiips()[1:3, ])], nrow = 3),
    frst.fsyips()[1:3, ],
    frst.fglps()[1:3, ]
  ))
  df[, 1:4] = as.integer(as.matrix(df[, 1:4]))
  names(df) = c("ID 1", "ID 2", "Survey index 1", "Survey index 2", gpkts)
  df
})

# Histograms of genopair log-probabilities given kinships
output$firstGLPs = renderPlot({
  # Set up four plots
  par(mfrow = c(1, 4))
  
  # Try to plot for each kinship
  lapply(1:4, function(i) {
    # If any log-probabilities greater than negative infinity
    if(any(is.finite(frst.fglps()[, i]))) {
      # Plot them
      hist(
        frst.fglps()[, i], main = gpkts[i],
        xlab = "Log-probability", br = 50
      )
    } 
    
    # Otherwise leave a blank space
    else plot.new()
  })
})

# Plot PLODs on regular scale
output$firstPLODs = renderPlot({
  # Plot plods
  hist(
    first.plods(), main = "HSP vs UP PLODs for all samples",
    sub = paste(L(), "loci"), xlab = "PLOD", breaks = 200,
  )
  
  # Plot expected values
  abline(v = exp.plods()[1], col = 2, lwd = 2)
  
  # Add legend
  legend(
    "topright", col = 1:7, lty = c(0, rep(1, 6)), cex = 1, lwd = 2,
    legend = c(
      "Expected value given kinship", "Unrelated", "First cousin",
      "Avuncular", "Half-sibling", "Parent-offspring", "Self"
    )
  )
})

# Plot PLODs on log-scale
output$firstPLODsRare = renderPlot({
  # Get counts and take logs
  hst.data = hist(first.plods(), breaks = 200, plot = F) 
  hst.data$counts = log(hst.data$counts + 1)
  
  # Plot them
  plot(
    hst.data, main = "Frequency on log scale",
    xlab = "PLOD", ylab = "Log (frequency + 1)"
  )
  
  # Add expected values
  abline(v = exp.plods(), col = 2:7, lwd = 2)
})

# Table of genopair probabilities of first few pairs captured
output$firstFewGPPsFll = renderTable({
  df = data.frame(cbind(
    matrix(frst.std()[frst.siis()[frst.fsiips()[1:3, ]], 1], ncol = 2),
    frst.fsyips()[1:3, ],
    format(head(frst.fgps(), 3), digits = 3, scientific = T)
  ))
  df[, 1:4] = as.integer(as.matrix(df[, 1:4]))
  names(df) = c("ID1", "ID2", "Survey index 1", "Survey index 2", gpkts)
  df
}, digits = 6)

# Show table of genopair probabilities for first few offset pairs (order is
# random to avoid bias due to age representation when individuals repeated to
# include pairs between surveys with different numbers of samples)
output$firstFewGPPsOffst = renderTable({
  df = data.frame(cbind(
    matrix(
      frst.std()[frst.siis()[frst.osiips()[1:3, ]], 1], 
      ncol = 2
    ),
    frst.osyips()[1:3, ],
    format(head(frst.ogps(), 3), digits = 3, scientific = T)
  ))
  df[, 1:4] = as.integer(as.matrix(df[, 1:4]))
  names(df) = c("ID1", "ID2", "Survey index 1", "Survey index 2", gpkts)
  df
}, digits = 6)

# Table of survey-year index pair counts for full genopair model
output$frstSYIPCntsFll = renderTable({
  df = data.frame(matrix(table(data.frame(frst.fsyips())), nrow = k()))
  dimnames(df) = list(Index_1 = 1:k() - 1, Index_2 = 1:k() - 1)
  df
}, rownames = T)

# Table of numbers of pairs with corresponding survey indices for first and
# second samples
output$frstSYIPCntsOffst = renderTable({
  df = data.frame(matrix(table(data.frame(frst.osyips())), nrow = k()))
  dimnames(df) = list(Index_1 = 1:k() - 1, Index_2 = 1:k() - 1)
  df
}, rownames = T)

# Histograms of genopair log-probabilities given basic kinships
output$frstGpPs = renderPlot({
  par(mfrow = c(1, 4))
  frst.Gp.Ps = FindGPPs(frst.fglps())
  lapply(1:4, function(i) {
    hst.dt = hist(frst.Gp.Ps[, i], br = 50, plot = F)
    hst.dt$counts = log(hst.dt$counts + 1)
    plot(
      hst.dt, main = gpkts[i], xlab = "Probability", 
      ylab = "log (frequency + 1)"
    )
  })
})


