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

# Find the three possible genotype probabilities at each locus as 3 x L matrix,
# by indexing 2 x L allele frequencies matrix for each allele of each possible
# genotype (globally defined for SNP genotypes), and multiplying by 2 possible
# cases for heterozygous genotypes
pss.gt.prbs = reactive({
  matrix(
    ale.frqs()[ales.1.inds, ] * ale.frqs()[ales.2.inds, ] * 
      (1 + (ales.1.inds != ales.2.inds)), 
    nrow = n.pss.gts, ncol = L()
  )
})

# Find PLODs
first.plods = reactive(FindPlods(smp.gts(), ale.frqs(), pss.gt.prbs(), L()))

# Plot PLODs
output$firstPLODs = renderPlot({
  # Plot plods
  hist(
    first.plods()$plods, 
    main = "HSP vs UP PLODs for all samples",
    sub = paste(L(), "loci"),
    xlab = "PLOD",
    breaks = 200,
  )
  
  # Plot expected values
  abline(v = first.plods()$ev.up, col = 2, lwd = 2)
  
  # Add legend
  legend(
    "topright", col = 1:7, lty = c(0, rep(1, 6)), cex = 1, lwd = 2,
    legend = c(
      "Expected value given kinship", "Unrelated", "First cousin",
      "Avuncular", "Half-sibling", "Parent-offspring", "Self"
    )
  )
})

# Plot rare values of PLODs
output$firstPLODsRare = renderPlot({
  # Plot uncommon values
  hist(
    first.plods()$plods, 
    main = "Uncommon values suggest likely close-kin",
    xlab = "PLOD",
    breaks = 200, 
    ylim = c(0, 100)
  )
  
  # Plot expected values
  abline(v = first.plods()$evs, col = 2:7, lwd = 2)
})

