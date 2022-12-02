# Table of genotypes of first few individuals captured (can show kin-pairs
# later)
output$firstGTs = renderTable({
  gt = attributes(fst.std())$unq.smp.gts
  df = data.frame(cbind(
    rep(fst.std()$ID[1:3], each = 2), 
    rep(paste0(c("m", "p"), "aternal"), 3),
    rbind(gt[, , 1], gt[, , 2], gt[, , 3])
  ))
  names(df) = c("ID", "Allele", paste0("L", 1:L()))
  df
})

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

# First genotypes, found by indexing the 2 x L allele frequencies matrix for
# each allele of each possible genotype (globally defined for SNP genotypes),
# and multiplying by 2 possible cases for heterozygous genotypes.  Filled into a
# 3 x 1 x L array which is indexed 3 times to fill the three columns.
pss.gt.1.prbs = reactive({
  array(
    ale.frqs()[ales.1.inds, ] * ale.frqs()[ales.2.inds, ] * 
      (1 + (ales.1.inds != ales.2.inds)), c(n.pss.gts, 1, L())
  )[, rep(1, 3), ]
})

# Conditional probabilities given that the pair are unrelated, the products of
# the respective genotype probabilities, the second found by permuting the array
# containing the first.
pss.gp.prbs.UP = reactive({
  pss.gt.1.prbs() * aperm(pss.gt.1.prbs(), c(2, 1, 3))
})

# Conditional probabilities of second genotype given that the pair are parent
# and offspring (unordered).  0.5 for each allele in first genotype being
# inherited, multiplied by the probability of the second genotype in each case,
# as in table 3, pg. 269, Bravingtion et al. (2016) Close-Kin Mark-Recapture.
# Filled into an L x 3 x 3 array which is permuted to the standard dimensions.
pss.cnd.gt.2.prbs.POP = reactive({
  aperm(
    array(
      # Note: order data enters array is down columns, not across rows
      c(
        ale.frqs()[1, ], 0.5 * ale.frqs()[1, ], rep(0, L()), 
        ale.frqs()[2, ], 0.5 * colSums(ale.frqs()), ale.frqs()[1, ],
        rep(0, L()), 0.5 * ale.frqs()[2, ], ale.frqs()[2, ]
      ), c(L(), n.pss.gts, n.pss.gts)
    ), c(2, 3, 1)
  )
})

# Conditional probabilities given that the pair are parent and offspring
# (unordered). Products of the first genotype probabilities and the conditional
# probabilities of the second genotypes given that the pair are parent and
# offspring.
pss.gp.prbs.POP = reactive({
  pss.gt.1.prbs() * pss.cnd.gt.2.prbs.POP()
})

# Conditional probabilities given that the pair are the same individual sampled
# twice. Genotype probabilities when the genotypes are the same, and zero
# otherwise.
pss.gp.prbs.SP = reactive({
  aperm(
    array(
      cbind(
        pss.gt.1.prbs()[1, 1, ], 0, 0, 0, pss.gt.1.prbs()[2, 2, ], 0, 0, 0, 
        pss.gt.1.prbs()[3, 3, ]
      ), c(L(), n.pss.gts, n.pss.gts)
    ), c(2, 3, 1)
  )
})

# Conditional probabilities given that the pair are half-siblings.  Average of
# probabilities for unrelated and parent-offspring pairs.
pss.gp.prbs.HSP = reactive({
  (pss.gp.prbs.UP() + pss.gp.prbs.POP()) / 2
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
  df = data.frame(gpps()[, , 1])
  names(df) = row.names(df) = c("00", "01", "11")
  df
}

# Table showing GPPs given unrelated
output$firstGPPsUP = renderTable({
  frmt.gpps(pss.gp.prbs.UP)
}, rownames = T)

# Table showing GPPs given half-siblings
output$firstGPPsHSP = renderTable({
  frmt.gpps(pss.gp.prbs.HSP)
}, rownames = T)

# Table showing GPPs given parent-offspring
output$firstGPPsPOP = renderTable({
  frmt.gpps(pss.gp.prbs.POP)
}, rownames = T)

# Table showing GPPs given self-resample
output$firstGPPsSP = renderTable({
  frmt.gpps(pss.gp.prbs.SP)
}, rownames = T)

# HSP vs UP pseudo log-likelihood ratios at each locus, calculated as on pg. 29
# of Aldridge-Sutton (2019) New Methods ..., but subtracting log(2) later, after
# summation, rather than now for each term separately. 
pss.plods.pls.lg.2 = reactive({
  log(1 + pss.gp.prbs.POP() / pss.gp.prbs.UP())
})

# Expected values of PLODs given kinships
exp.plod.KP = reactive({
  # Unrelated, parent-offspring, and self-pairs
  exp.plod.base = sapply(
    list(pss.gp.prbs.UP(), pss.gp.prbs.POP(), pss.gp.prbs.SP()),
    function(gp.prbs) sum(gp.prbs * pss.plods.pls.lg.2()) / L()
  )
  
  # First-cousin, avuncular, and half-sibling pairs
  exp.plod.extd = 
    (c(7, 3, 1) * exp.plod.base[1] + exp.plod.base[2]) / c(8, 4, 2)
  
  # Combine in order from furthest to closest kinship, add names, and return.
  # P(gts|HSP) = (P(gts|UP) + P(gts|POP)) / 2 at each locus, so subtracting
  # log(2) * n.loci divided by n.loci.
  vec = c(exp.plod.base[1], exp.plod.extd, exp.plod.base[2:3]) - log(2)
  names(vec) = c(
    "Unrelated", "First cousin", "Avuncular", "Half-sibling",
    "Parent-offspring", "Self"
  )
  vec
})

fst.smp.hst = reactive({
  as.matrix(fst.std()[, 4:(3 + k())])
})
smp.inds = reactive({
  row(fst.smp.hst())[as.logical(fst.smp.hst())]
})
smp.ind.prs = reactive({
  t(combn(smp.inds(), 2))
})

# Add genotypes for repeated samples. Unique sampled genotypes are 2 x L x
# n_animals arrays, representing two binary SNPs at each locus for each
# individual sampled at least once. They are expanded to 2 x L x n_samples
# arrays by repeating genotypes by numbers of samples per individual. Also order
# by survey-year for easy look up of kinship probabilities in TMB later.
smp.gts = reactive({
  attributes(fst.std())$unq.smp.gts[, , smp.inds()]
})

# Indices of survey-years for each sample in each pair, starting at zero for
# C++ template, and ordered by survey-year of first sample
smp.yr.ind.prs = reactive({
  t(combn(col(fst.smp.hst())[as.logical(fst.smp.hst())] - 1, 2))
})

# Indices of within-survey pairs
wtn.prs.inds = reactive({
  smp.yr.ind.prs()[, 1] == smp.yr.ind.prs[, 2]()
})

# Nullify genopair log-probabilities when new datasets simulated
observeEvent(input$simulate, lg.gp.prbs.KP(NULL))

# Find genopair log-probabilities as n.pairs x n.kp.tps matrix with rows for
# pairs of individuals, and columns for types of kinships considered
observeEvent({
  input$nav.tab
  input$check.sub.tabs
}, {
  if (
    input$nav.tab == "check.tab" && input$check.sub.tabs == "frst.gts.tb" &&
    is.null(lg.gp.prbs.KP())
  ) {
    lg.gp.prbs.KP(FindLogGPProbsKP(
      smp.gts(), L(), pss.gp.prbs.UP(), pss.gp.prbs.POP(), pss.gp.prbs.SP(),
      wtn.prs.inds
    ))
  }
})

gp.prbs.KP = reactive(exp(lg.gp.prbs.KP()))

# Find half-sibling vs unrelated pairs PLODs from log genopair probabilities
first.plods = reactive({
  (lg.gp.prbs.KP()[, 2] - lg.gp.prbs.KP()[, 1]) / L()
})

# Table of genopair probabilities of first few pairs captured (can show
# kin-pairs later)
output$firstGPPs = renderTable({
  df = data.frame(cbind(
    fst.std()[smp.ind.prs()[1:3, 1], 1],
    fst.std()[smp.ind.prs()[1:3, 2], 1],
    smp.yr.ind.prs()[1:3, ],
    gp.prbs.KP()[1:3, ],
    lg.gp.prbs.KP()[1:3, ]
  ))
  names(df) = c(
    "ID1", "ID2", "Survey index 1", "Survey index 2",
    paste0("P_", colnames(gp.prbs.KP())), 
    paste0("Log_P_", colnames(gp.prbs.KP()))
  )
  df
})

output$firstObsGPPs = renderPlot({
  par(mfrow = c(1, 4))
  br = 50
  xlab = "Log-probability"
  hist(lg.gp.prbs.KP()[, "UP"], main = "Unrelated", xlab = xlab, br = br)
  hist(
    lg.gp.prbs.KP()[, "HSP"], main = "Half-sibling", xlab = xlab, br = br
  )
  hist(
    lg.gp.prbs.KP()[, "POP"], main = "Parent-offspring", xlab = xlab, br = br
  )
  hist(
    lg.gp.prbs.KP()[, "SP"], main = "Self-resample", xlab = xlab, br = br
  )
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
  hist(
    first.plods(), main = "Uncommon values suggest likely close-kin",
    xlab = "PLOD", breaks = 200, ylim = c(0, 100)
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
  
  # Try to git genopair likelihood model
  gp.tmb = TryGenopairTMB(
    if (any(all_undrflw)) gpp.adj else gpp.slct, 
    smp.yr.ind.prs(),
    # smp.yr.ind.prs()[wtn.prs.inds, ], 
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), ck.start, ck.lwr, ck.upr, alpha()
  )
  
  gp.tmb$est.se.df
})
