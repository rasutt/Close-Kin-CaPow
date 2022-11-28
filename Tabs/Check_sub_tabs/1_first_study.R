# Outputs for first study sub-tab of checks tab

# Get first study from list
fst.std = reactive(sim.lst()$hists.lst[[1]])

### First study simulated

# Function to format table of integers
frmt.tbl = function(data, rw.nms, cl.nms) {
  mode(data) = "integer"
  df = data.frame(matrix(data, ncol = length(cl.nms)), row.names = rw.nms)
  names(df) = cl.nms
  df
}

## First sample-histories
output$firstSampHists = renderTable({
  head(data.frame(fst.std()))
})

## Numbers of kin-pairs in whole population
# Within surveys
output$firstNsKPsWtnPop = renderTable({
  pop.atts = attributes(fst.std())
  SMPs = find.SMPs.wtn(pop.atts, k())
  SFPs = find.SFPs.wtn(pop.atts, k())
  FSPs = find.FSPs.wtn(pop.atts, k())
  frmt.tbl(
    rbind(
      pop.atts$N.t.vec[s.yr.inds()],
      choose(pop.atts$N.t.vec[s.yr.inds()], 2),
      find.POPs.wtn(pop.atts, k()),
      SMPs, SFPs, FSPs, SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps.pop.wtn, srvy.yrs()
  )
}, rownames = T)

# Between surveys
output$firstNsKPsBtnPop = renderTable({
  pop.atts = attributes(fst.std())
  SPs.prnts.kwn = find.SPs.prnts.kwn(pop.atts, k())
  SMPs = find.SMSPs.btn(pop.atts, k()) - SPs.prnts.kwn
  SFPs = find.SFSPs.btn(pop.atts, k()) - SPs.prnts.kwn
  FSPs = find.FSSPs.btn(pop.atts, k()) - SPs.prnts.kwn
  frmt.tbl(
    rbind(
      combn(
        pop.atts$N.t.vec[s.yr.inds()], 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2]
      ),
      find.SPs(pop.atts, k()),  SPs.prnts.kwn, find.POPs.btn(pop.atts, k()),
      SMPs, SFPs, FSPs, SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps.pop.btn, srvy.prs()
  )
}, rownames = T)

## Estimated numbers of kin-pairs in whole population
# Within surveys
output$firstEstNsKPsWtnPop = renderTable({
  frmt.tbl(t(est.ns.kps.pop()$wtn), kp.tps.pop.wtn, srvy.yrs())
}, rownames = T)

# Between surveys
output$firstEstNsKPsBtnPop = renderTable({
  frmt.tbl(t(est.ns.kps.pop()$btn), kp.tps.pop.btn, srvy.prs())
}, rownames = T)

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

## HSP vs UP PLODs for pairs of individuals captured

# Add genotypes for repeated samples. Unique sampled genotypes are 2 x L x
# n_animals arrays, representing two binary SNPs at each locus for each
# individual sampled at least once. They are expanded to 2 x L x n_samples
# arrays by repeating genotypes by numbers of samples per individual.
smp.gts = reactive({
  # Numbers of samples per individual are sums of rows of binary sample-history
  # matrix
  ns.smps.pr.ind = rowSums(fst.std()[, 4:(3 + k())])
  
  # Index unique genotype array repeatedly for re-sampled individuals
  attributes(fst.std())$unq.smp.gts[
    , , rep(1:length(ns.smps.pr.ind), times = ns.smps.pr.ind)
  ]
})

# Find allele frequencies from genotypes in 2 x L x n_samples arrays,
# representing two binary SNPs at each locus for each sample.  Frequencies are
# returned as 2 x L matrices representing the frequencies of 0 and 1-coded SNP
# alleles at each locus
ale.frqs = reactive({
  # Frequencies of 1-coded SNP alleles are means over both alleles at each locus
  # for each sample
  ale.frqs.1 = apply(smp.gts(), 2, mean)
  
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
      (1 + (ales.1.inds != ales.2.inds)), 
    c(n.pss.gts, 1, L())
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
      c(
        ale.frqs()[1, ], 0.5 * ale.frqs()[1, ], rep(0, L()), 
        ale.frqs()[2, ], 0.5 * colSums(ale.frqs()), ale.frqs()[1, ],
        rep(0, L()), 0.5 * ale.frqs()[2, ], ale.frqs()[2, ]
      ),
      c(L(), n.pss.gts, n.pss.gts)
    ), 
    c(2, 3, 1)
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
      ), 
      c(L(), n.pss.gts, n.pss.gts)
    ),
    c(2, 3, 1)
  )
})

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

# Find log genopair probabilities as n.pairs x n.kp.tps matrix with rows for
# pairs of individuals, and columns for types of kinships considered
lg.gp.prbs.KP = reactive({
  FindLogGPProbsKP(
    smp.gts(), L(), pss.gp.prbs.POP(), pss.gp.prbs.UP(), pss.gp.prbs.SP()
  )
})

# Find half-sibling vs unrelated pairs PLODs from log genopair probabilities
first.plods = reactive({
  (lg.gp.prbs.KP()[, 2] - lg.gp.prbs.KP()[, 1]) / L()
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
  hist(
    first.plods(), main = "Uncommon values suggest likely close-kin",
    xlab = "PLOD", breaks = 200, ylim = c(0, 100)
  )
  
  # Plot expected values
  abline(v = exp.plod.KP(), col = 2:7, lwd = 2)
})

