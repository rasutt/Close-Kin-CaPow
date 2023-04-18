# Genotypes, 2 x L x n_individuals, indexed from unique sampled genotypes,
# (2 x L x n_animals), representing two binary SNPs at each locus for
# each individual sampled at least once. They are expanded to arrays by indexing
# sample genotypes.
frst.gts = reactive(FS.atts()$ind.gts)

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
frst.osips = reactive(FindOSIPs(k(), frst.syis()))

# Function to find sample-individual or sample-year index pairs, n_pairs x 2
FindSISYIPs = function(sisyis, sips) {
  matrix(sisyis[as.vector(sips)], ncol = 2)
}

# Full and offset sample-individual index pairs, n_pairs x 2
frst.fsiips = reactive(t(combn(frst.siis(), 2)))
frst.osiips = reactive(FindSISYIPs(frst.siis(), frst.osips()))

# Indices of survey-years for each sample in each pair, starting at zero for
# C++ template, and ordered by survey-year of first sample
frst.fsyips = reactive(t(combn(frst.syis(), 2)))
frst.osyips = reactive(FindSISYIPs(frst.syis(), frst.osips()))

# Nullify genopair log-probabilities when new datasets simulated
observeEvent(input$simulate, frst.fglps(NULL))

# Full genopair log-probabilities

# If tab changed in app
observeEvent({
  input$nav.tab
  input$check.sub.tabs
  input$FS.sub.tab
}, {
  # If new tab requires full genopair probabilities and they have not been
  # computed
  if (
    input$nav.tab == "check.tab" && 
    input$check.sub.tabs == "frst.stdy.tb" && 
    input$FS.sub.tab %in% c(
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
frst.fgps = reactive(FindGPs(frst.fglps()))
frst.ogps = reactive(FindGPs(frst.oglps()))

# Expected values of HSP vs UP PLODs given various kinships
exp.plods = reactive({
  # Possible genopair log-probability ratios (given HSP vs UP) at each locus,
  # 3 x 3 x n_loci
  pglprs = log(frst.pgpsgks()[, , , 2] / frst.pgpsgks()[, , , 1])
  
  # Expected plods for kinship basis, unrelated, half-sibling, parent-offspring,
  # and self-pairs
  eplds.bss = colSums(frst.pgpsgks() * rep(pglprs, 4), dims = 3) / L()
  
  # Expected plods for extended kinships, first-cousin and avuncular pairs
  eplds.extd = (c(7, 3) * eplds.bss[1] + eplds.bss[3]) / c(8, 4)
  
  # Combine in order from furthest to closest kinship, add names, and return
  eplds.cmbd = c(eplds.bss[1], eplds.extd, eplds.bss[2:4])
  names(eplds.cmbd) = c(
    "Unrelated", "First cousin", "Avuncular", "Half-sibling",
    "Parent-offspring", "Self"
  )
  eplds.cmbd
})

# Find half-sibling vs unrelated pairs PLODs from log genopair probabilities
first.plods = reactive({
  (frst.fglps()[, 2] - frst.fglps()[, 1]) / L()
})

## Outputs

# Table of genotypes of first few individuals captured
output$firstFSGs = renderTable({
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
output$firstGPs = renderTable({
  df = data.frame(FindPssGtPrbs(frst.ale.frqs()))
  names(df) = paste0("L", 1:L())
  row.names(df) = pss.gt.lbls
  df
}, rownames = T)

# Function to format genopair probabilities and log-probabilities at first locus
# for display
FrmtGDL1 = function(gnpr.dt) {
  df = data.frame(gnpr.dt)
  names(df) = row.names(df) = pss.gt.lbls
  df
}

# Tables showing genopair probabilities given each kinship
output$firstGPsUPL1 = renderTable({
  FrmtGDL1(frst.pgpsgks()[, , 1, 1])
}, rownames = T)
output$firstGPsHSPL1 = renderTable({
  FrmtGDL1(frst.pgpsgks()[, , 1, 2])
}, rownames = T)
output$firstGPsPOPL1 = renderTable({
  FrmtGDL1(frst.pgpsgks()[, , 1, 3])
}, rownames = T)
output$firstGPsSPL1 = renderTable({
  FrmtGDL1(frst.pgpsgks()[, , 1, 4])
}, rownames = T)

# Table showing genopair log-probabilities given each kinship
output$firstGLPsUPL1 = renderTable({
  FrmtGDL1(log(frst.pgpsgks()[, , 1, 1]))
}, rownames = T)
output$firstGLPsHSPL1 = renderTable({
  FrmtGDL1(log(frst.pgpsgks()[, , 1, 2]))
}, rownames = T)
output$firstGLPsPOPL1 = renderTable({
  FrmtGDL1(log(frst.pgpsgks()[, , 1, 3]))
}, rownames = T)
output$firstGLPsSPL1 = renderTable({
  FrmtGDL1(log(frst.pgpsgks()[, , 1, 4]))
}, rownames = T)

# Function to format first sample genopair probabilities and log-probabilities
# for display
FrmtFSGD = function(siips, gnpr.dt) {
  df = data.frame(
    matrix(frst.std()$ID[as.vector(siips)], nrow = 3),
    format(gnpr.dt, digits = 3, scientific = T)
  )
  names(df) = c("ID 1", "ID 2", gpkts)
  df
}

# Table of genopair log-probabilities of first few pairs captured
output$firstFSGLPs = renderTable({
  FrmtFSGD(frst.fsiips()[1:3, ], frst.fglps()[1:3, ])
})

# Histograms of genopair log-probabilities given kinships
output$firstASGLPs = renderPlot({
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
    
    # Otherwise leave a blank plot
    else plot.new()
  })
})

# Tables of genopair probabilities of first few pairs in full and offset models
output$firstFSGPsF = renderTable({
  FrmtFSGD(frst.fsiips()[1:3, ], frst.fgps()[1:3, ])
})
output$firstFSGPsO = renderTable({
  FrmtFSGD(frst.osiips()[1:3, ], frst.ogps()[1:3, ])
})

# Tables of survey-year index pair counts for full and offset genopair models
FrmtSYIPs = function(syips) {
  df = data.frame(matrix(table(data.frame(syips)), nrow = k()))
  dimnames(df) = list(1:k() - 1, 1:k() - 1)
  df
}
output$firstSYIPsF = renderTable(FrmtSYIPs(frst.fsyips()), rownames = T)
output$firstSYIPsO = renderTable(FrmtSYIPs(frst.osyips()), rownames = T)

# Plot PLODs on log-scale
output$firstPLODs = renderPlot({
  # Get counts and take logs
  hst.data = hist(first.plods(), breaks = 200, plot = F) 
  hst.data$counts = log(hst.data$counts + 1)
  
  # Plot them
  plot(
    hst.data, main = "HSP vs UP PLODs for all samples",
    xlab = "PLOD", ylab = "Log (frequency + 1)"
  )
  
  # Add expected values
  abline(v = exp.plods(), col = 2:7, lwd = 2)
  
  # Add legend
  legend(
    "topright", col = 1:7, lty = c(0, rep(1, 6)), cex = 1, lwd = 2,
    legend = c(
      "Expected value given kinship", "Unrelated", "First cousin",
      "Avuncular", "Half-sibling", "Parent-offspring", "Self"
    )
  )
})


