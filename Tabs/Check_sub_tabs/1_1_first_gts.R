# Find allele frequencies from genotypes in 2 x L x n_individuals arrays,
# representing two binary SNPs at each locus for each individual, excluding
# repeated samples of the same individual in different surveys.  Frequencies are
# returned as 2 x L matrices representing the frequencies of 0 and 1-coded SNP
# alleles at each locus
frst.ale.frqs = reactive({
  # Frequencies of 1-coded SNP alleles are means over both alleles at each locus
  # for each sample
  ale.frqs.1 = apply(attributes(fst.std())$unq.smp.gts, 2, mean)
  # ale.frqs.1 = apply(
  #   mpfrArray(
  #     attributes(fst.std())$unq.smp.gts, 10, 
  #     dim(attributes(fst.std())$unq.smp.gts)
  #   ), 2, mean
  # )
  
  # Combine with frequencies for 0-coded alleles
  rbind(1 - ale.frqs.1, ale.frqs.1)
})

# Probabilities of possible genopairs at each locus as 3 x 3 x L arrays, where
# rows represent the first genotypes, and columns the second, ordered as 00, 01,
# and 11, for binary SNPs.
frst.pss.gp.prbs.KPs = reactive({
  find.pss.gp.prbs.KPs(pss.gts, n.pss.gts, frst.ale.frqs(), L())
})

# Sample histories from first study, (n_animals x n_surveys)
frst.smp.hsts = reactive(as.matrix(fst.std()[, 4:(3 + k())]))

# Indices of samples from first study, (n_samples)
frst.smp.inds = reactive(row(frst.smp.hsts())[as.logical(frst.smp.hsts())])

# Sample genotypes, (2 x L x n_samples), indexed from unique sampled genotypes,
# (2 x L x n_animals), representing two binary SNPs at each locus for
# each individual sampled at least once. They are expanded to arrays by indexing
# sample genotypes.
frst.smp.gts = reactive(attributes(fst.std())$unq.smp.gts[, , frst.smp.inds()])

# Indices of pairs of samples from first study, (n_samples x 2)
frst.smp.ind.prs = reactive(t(combn(frst.smp.inds(), 2)))

# Indices of survey-years for each sample in each pair, starting at zero for
# C++ template, and ordered by survey-year of first sample
frst.smp.yr.ind.prs = reactive({
  t(combn(col(frst.smp.hsts())[as.logical(frst.smp.hsts())] - 1, 2))
})

# Indices of within-survey pairs to reduce search space later
frst.wtn.prs.inds = reactive(which(frst.smp.yr.ind.prs()[, 1] == frst.smp.yr.ind.prs[, 2]()))

# Nullify genopair log-probabilities when new datasets simulated
observeEvent(input$simulate, frst.lg.gp.prbs.KP(NULL))

# Find genopair log-probabilities as n.pairs x n.kp.tps matrix with rows for
# pairs of individuals, and columns for types of kinships considered
observeEvent(
  {
    input$nav.tab
    input$check.sub.tabs
  }, 
  {
    if (
      input$nav.tab == "check.tab" && input$check.sub.tabs == "frst.gts.tb" &&
      is.null(frst.lg.gp.prbs.KP())
    ) {
      frst.lg.gp.prbs.KP(
        FindLogGPProbsKP(frst.smp.gts(), L(), frst.pss.gp.prbs.KPs(), frst.wtn.prs.inds)
      )
    }
  }
)

# Expected values of HSP vs UP PLODs given kinships
exp.plod.KP = reactive({
  # Possible values of HSP vs UP PLODs, leaving division by number of loci to
  # next step after summation
  pss.plods = log(frst.pss.gp.prbs.KPs()[, , , 2] / frst.pss.gp.prbs.KPs()[, , , 1])
  
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
first.plods = reactive((frst.lg.gp.prbs.KP()[, 2] - frst.lg.gp.prbs.KP()[, 1]) / L())

## Outputs

# First sample-histories
output$firstSampHists = renderTable(head(data.frame(fst.std())))

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
  # df = data.frame(gpps)
  df = data.frame(asNumeric(gpps))
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
    fst.std()[frst.smp.ind.prs()[1:3, 1], 1],
    fst.std()[frst.smp.ind.prs()[1:3, 2], 1],
    frst.smp.yr.ind.prs()[1:3, ],
    frst.lg.gp.prbs.KP()[1:3, ]
  ))
  names(df) = c("ID1", "ID2", "Survey index 1", "Survey index 2", gp.prb.KP.tps)
  df
})

# Histograms of genopair log-probabilities given basic kinships
output$firstLGPPs = renderPlot({
  par(mfrow = c(1, 4))
  br = 50
  xlab = "Log-probability"
  lapply(1:4, function(i) {
    hist(frst.lg.gp.prbs.KP()[, i], main = gp.prb.KP.tps[i], xlab = xlab, br = br)
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

# Table of estimates from first study optimising genopair likelihood
output$firstGPEsts = renderTable({
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), attributes(fst.std())$N.t.vec[hist.len()])
  ck.lwr <- c(0, 0.75, attributes(fst.std())$ns.caps[k()])
  ck.upr <- c(0.35, 1, Inf)
  
  # Get genopair probabilities (by excluding probabilities giveb half-sibs for
  # now) and check for pairs where all probabilities underflow to zero
  lg.gpp.slct = frst.lg.gp.prbs.KP()[, -2]
  gpp.slct = exp(lg.gpp.slct)
  all_undrflw = rowSums(gpp.slct) == 0
  
  # If there is underflow adjust log-probabilities by factor giving equal
  # weight to smallest and largest values to avoid both under and overflow
  if (any(all_undrflw)) {
    cat("Proportion of pairs for which probabilities given all kinships 
        underflow to zero:", mean(all_undrflw), "\n")
    
    # Want smallest maximum kinship probability and largest probability to be
    # equally far from one
    adj = mean(c(min(apply(lg.gpp.slct, 1, max)), max(lg.gpp.slct)))
    lg.gpp.adj = lg.gpp.slct - adj
    gpp.adj = exp(lg.gpp.adj)
    
    # Show adjustment and results
    cat("Probabilities adjusted by factor of exp(", adj, ")\n", sep = "")
    print("Adjusted log-probabilities and probabilities of genopairs:")
    print(summary(lg.gpp.adj))
    print(summary(gpp.adj))
    cat("Proportion of pairs for which adjusted probabilities given all 
        kinships underflow to zero:", mean(rowSums(gpp.adj) == 0), "\n")
  }
  
  # Try to fit genopair likelihood model
  print(table(frst.smp.yr.ind.prs()[, 1], frst.smp.yr.ind.prs()[, 2]))
  print(str(gpp.slct))
  gp.tmb = TryGenopairTMB(
    if (any(all_undrflw)) gpp.adj else gpp.slct, 
    frst.smp.yr.ind.prs(),
    # smp.yr.ind.prs()[wtn.prs.inds, ], 
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), ck.start, ck.lwr, ck.upr, alpha()
  )
  
  # # Checked that it give the same results in R
  # gp.r = TryGenopair(
  #   if (any(all_undrflw)) gpp.adj else gpp.slct, 
  #   smp.yr.ind.prs() + 1,
  #   k(), fnl.year(), srvy.yrs(), alpha(), 
  #   ck.start, ck.lwr, ck.upr
  # )
  # print(gp.r)
  
  gp.tmb$est.se.df
})

