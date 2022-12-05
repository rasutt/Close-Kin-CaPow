# Find allele frequencies from genotypes in 2 x L x n_individuals arrays,
# representing two binary SNPs at each locus for each individual, excluding
# repeated samples of the same individual in different surveys.  Frequencies are
# returned as 2 x L matrices representing the frequencies of 0 and 1-coded SNP
# alleles at each locus
ale.frqs = reactive({
  # Frequencies of 1-coded SNP alleles are means over both alleles at each locus
  # for each sample
  ale.frqs.1 = apply(attributes(fst.std())$unq.smp.gts, 2, mean)
  
  # Combine with frequencies for 0-coded alleles
  rbind(1 - ale.frqs.1, ale.frqs.1)
})

# Probabilities of possible genopairs at each locus as 3 x 3 x L arrays, where
# rows represent the first genotypes, and columns the second, ordered as 00, 01,
# and 11, for binary SNPs.
pss.gp.prbs.KP = reactive({
  ales.1.inds = pss.gts[1, ]
  ales.2.inds = pss.gts[2, ]
  
  # First genotypes, found by indexing the 2 x L allele frequencies matrix for
  # each allele of each possible genotype (globally defined for SNP genotypes),
  # and multiplying by 2 possible cases for heterozygous genotypes.  Filled into
  # a 3 x 1 x L array which is indexed 3 times to fill the three columns.
  pss.gt.prbs = matrix(
    ale.frqs()[ales.1.inds, ] * ale.frqs()[ales.2.inds, ] * 
      (1 + (ales.1.inds != ales.2.inds)), 
    nrow = 3
  )
  pss.gt.2.prbs = array(
    rep(pss.gt.prbs, each = 3), 
    c(n.pss.gts, n.pss.gts, L())
  )
  pss.gt.1.prbs = aperm(pss.gt.2.prbs, c(2, 1, 3))
  
  # Conditional probabilities given that the pair are unrelated, the products of
  # the respective genotype probabilities, the second found by permuting the
  # array containing the first.
  pss.gp.prbs.UP = pss.gt.1.prbs * pss.gt.2.prbs
  
  # Conditional probabilities of second genotype given that the pair are parent
  # and offspring (unordered).  0.5 for each allele in first genotype being
  # inherited, multiplied by the probability of the second genotype in each
  # case, as in table 3, pg. 269, Bravingtion et al. (2016) Close-Kin
  # Mark-Recapture. Filled into an L x 3 x 3 array which is permuted to the
  # standard dimensions.
  pss.cnd.gt.2.prbs.POP = aperm(
    array(
      # Note: order data enters array is down columns, not across rows
      c(
        ale.frqs()[1, ], 0.5 * ale.frqs()[1, ], rep(0, L()), 
        ale.frqs()[2, ], 0.5 * colSums(ale.frqs()), ale.frqs()[1, ],
        rep(0, L()), 0.5 * ale.frqs()[2, ], ale.frqs()[2, ]
      ), c(L(), n.pss.gts, n.pss.gts)
    ), c(2, 3, 1)
  )
  
  # Conditional probabilities given that the pair are parent and offspring
  # (unordered). Products of the first genotype probabilities and the
  # conditional probabilities of the second genotypes given that the pair are
  # parent and offspring.
  pss.gp.prbs.POP = pss.gt.1.prbs * pss.cnd.gt.2.prbs.POP
  
  # Conditional probabilities given that the pair are the same individual
  # sampled twice. Genotype probabilities when the genotypes are the same, and
  # zero otherwise.
  pss.gp.prbs.SP = aperm(
    array(
      cbind(
        pss.gt.1.prbs[1, 1, ], 0, 0, 0, pss.gt.1.prbs[2, 2, ], 0, 0, 0, 
        pss.gt.1.prbs[3, 3, ]
      ), c(L(), n.pss.gts, n.pss.gts)
    ), c(2, 3, 1)
  )
  
  # Conditional probabilities given that the pair are half-siblings.  Average of
  # probabilities for unrelated and parent-offspring pairs.
  pss.gp.prbs.HSP = (pss.gp.prbs.UP + pss.gp.prbs.POP) / 2

  array(
    c(pss.gp.prbs.UP, pss.gp.prbs.HSP, pss.gp.prbs.POP, pss.gp.prbs.SP),
    c(n.pss.gts, n.pss.gts, L(), 4)
  )
})

# Sample histories from first study, (n_animals x n_surveys)
fst.smp.hsts = reactive(as.matrix(fst.std()[, 4:(3 + k())]))

# Indices of samples from first study, (n_samples)
smp.inds = reactive(row(fst.smp.hsts())[as.logical(fst.smp.hsts())])

# Sample genotypes, (2 x L x n_samples), indexed from unique sampled genotypes,
# (2 x L x n_animals), representing two binary SNPs at each locus for
# each individual sampled at least once. They are expanded to arrays by indexing
# sample genotypes.
smp.gts = reactive(attributes(fst.std())$unq.smp.gts[, , smp.inds()])

# Indices of pairs of samples from first study, (n_samples x 2)
smp.ind.prs = reactive(t(combn(smp.inds(), 2)))

# Indices of survey-years for each sample in each pair, starting at zero for
# C++ template, and ordered by survey-year of first sample
smp.yr.ind.prs = reactive({
  t(combn(col(fst.smp.hsts())[as.logical(fst.smp.hsts())] - 1, 2))
})

# Indices of within-survey pairs to reduce search space later
wtn.prs.inds = reactive(which(smp.yr.ind.prs()[, 1] == smp.yr.ind.prs[, 2]()))

# Nullify genopair log-probabilities when new datasets simulated
observeEvent(input$simulate, lg.gp.prbs.KP(NULL))

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
      is.null(lg.gp.prbs.KP())
    ) {
      lg.gp.prbs.KP(
        FindLogGPProbsKP(smp.gts(), L(), pss.gp.prbs.KP(), wtn.prs.inds)
      )
    }
  }
)

# First sample-histories
output$firstSampHists = renderTable(head(data.frame(fst.std())))

# Table of genotypes of first few individuals captured (can show kin-pairs
# later)
output$firstGTs = renderTable({
  df = data.frame(cbind(
    rep(fst.std()$ID[smp.inds()[1:3]], each = 2), 
    rep(paste0(c("m", "p"), "aternal"), 3),
    matrix(smp.gts()[, , 1:3], nrow = 6)
  ))
  names(df) = c("ID", "Allele", paste0("L", 1:L()))
  df
})

# Table showing allele frequencies
output$firstAFs = renderTable({
  df = data.frame(ale.frqs())
  names(df) = paste0("L", 1:L())
  row.names(df) = 0:1
  df
}, rownames = T)

# Function to format genopair probabilities for display
frmt.gpps = function(gpps) {
  df = data.frame(gpps)
  names(df) = row.names(df) = c("00", "01", "11")
  df
}

# Table showing GPPs given unrelated
output$firstGPPsUP = renderTable({
  frmt.gpps(pss.gp.prbs.KP()[, , 1, 1])
}, rownames = T)

# Table showing GPPs given half-siblings
output$firstGPPsHSP = renderTable({
  frmt.gpps(pss.gp.prbs.KP()[, , 1, 2])
}, rownames = T)

# Table showing GPPs given parent-offspring
output$firstGPPsPOP = renderTable({
  frmt.gpps(pss.gp.prbs.KP()[, , 1, 3])
}, rownames = T)

# Table showing GPPs given self-resample
output$firstGPPsSP = renderTable({
  frmt.gpps(pss.gp.prbs.KP()[, , 1, 4])
}, rownames = T)

# Table of genopair probabilities of first few pairs captured (can show
# kin-pairs later)
output$firstFewLGPPs = renderTable({
  df = data.frame(cbind(
    fst.std()[smp.ind.prs()[1:3, 1], 1],
    fst.std()[smp.ind.prs()[1:3, 2], 1],
    smp.yr.ind.prs()[1:3, ],
    lg.gp.prbs.KP()[1:3, ]
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
    hist(lg.gp.prbs.KP()[, i], main = gp.prb.KP.tps[i], xlab = xlab, br = br)
  })
})

# Expected values of HSP vs UP PLODs given kinships
exp.plod.KP = reactive({
  # Possible values of HSP vs UP PLODs, leaving division by number of loci to
  # next step after summation
  pss.plods = log(pss.gp.prbs.KP()[, , , 2] / pss.gp.prbs.KP()[, , , 1])
  
  # Unrelated, parent-offspring, and self-pairs
  exp.plod.base = colSums(pss.gp.prbs.KP() * rep(pss.plods, 4), dims = 3) / L()
  
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
first.plods = reactive((lg.gp.prbs.KP()[, 2] - lg.gp.prbs.KP()[, 1]) / L())

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
  lg.gpp.slct = lg.gp.prbs.KP()[, -2]
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
  print(table(smp.yr.ind.prs()[, 1], smp.yr.ind.prs()[, 2]))
  print(str(gpp.slct))
  gp.tmb = TryGenopairTMB(
    if (any(all_undrflw)) gpp.adj else gpp.slct, 
    smp.yr.ind.prs(),
    # smp.yr.ind.prs()[wtn.prs.inds, ], 
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), ck.start, ck.lwr, ck.upr, alpha()
  )
  
  gp.tmb$est.se.df
})

# Plot GPPs - This is hard, but I do think it may be worthwhile as gives a
# feeling for what's happening in multiple dimensions. The log-probabilities are
# no good because most go to negative infinity for PO and SPs, but the GPPs
# should still be worth a shot.

# output$firstObsGPPs = renderPlot({
#   lims = apply(gp.prbs.KP(), 2, max)
#   x = seq(0, lims[1], len = 10)
#   y = seq(0, lims[2], len = 10)
#   z = outer(x, y, function(x, y) pmin(pmax(1 - x - y, 0), lims[4]))
#   
#   # Plot GPPs
#   res = persp(
#     x, y, z, main = "Genopair probabilities for all samples",
#     xlab = "Unrelated", ylab = "Parent-offspring", zlab = "Self",
#     xlim = c(0, lims[1]), ylim = c(0, lims[2]), zlim = c(0, lims[4]),
#     tick = "detailed"
#   )
#   points(
#     trans3d(
#       gp.prbs.KP()[, 1], gp.prbs.KP()[, 3], gp.prbs.KP()[, 4], res
#     ), col = rgb(1, 0, 0, 0.1)
#   )  
#   # # Plot expected values
#   # abline(v = exp.plod.KP()[1], col = 2, lwd = 2)
#   
#   # # Add legend
#   # legend(
#   #   "topright", col = 1:7, lty = c(0, rep(1, 6)), cex = 1, lwd = 2,
#   #   legend = c(
#   #     "Expected value given kinship", "Unrelated", "First cousin",
#   #     "Avuncular", "Half-sibling", "Parent-offspring", "Self"
#   #   )
#   # )
# })
