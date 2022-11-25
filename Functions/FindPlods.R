# Function to find PLODS for SNP genotypes.  Simplified from PLOD code for Emma.
FindPlods = function(smp.gts, ale.frqs, pss.gt.prbs, L) {
  # Repeat possible genotype probabilities to form 3 x 3 x L array representing
  # possible second genotypes as columns in possible genopairs at each locus
  pss.gt.2.prbs = array(
    rep(pss.gt.prbs, each = n.pss.gts), dim = c(n.pss.gts, n.pss.gts, L)
  )
  
  # Create array for possible half-sibling versus unrelated pair plods and add
  # possible second genotype probabilities
  hsp.up.plods.ary = pss.gt.2.prbs
  
  # Find possible genopairs
  gts.1 = pss.gts[, rep(1:n.pss.gts, n.pss.gts)]
  gts.2 = pss.gts[, rep(1:n.pss.gts, each = n.pss.gts)]
  
  # Find allele equalities among possible genopairs
  cis.eqs = gts.1 == gts.2
  cis.eqs.1.mat = matrix(cis.eqs[1, ], nrow = n.pss.gts)
  cis.eqs.2.mat = matrix(cis.eqs[2, ], nrow = n.pss.gts)
  
  trans.eqs.gts.2.htro = gts.1 == gts.2[c(2, 1), ] & 
    rep(gts.2[1, ] != gts.2[2, ], each = 2)
  trans.eqs.gts.2.htro.1.mat = matrix(trans.eqs.gts.2.htro[1, ], nrow = n.pss.gts)
  trans.eqs.gts.2.htro.2.mat = matrix(trans.eqs.gts.2.htro[2, ], nrow = n.pss.gts)
  
  # Add conditional probabilities of possible second genotypes, given
  # parent-offspring with first, to PLODs, when alleles shared
  hsp.up.plods.ary[rep(cis.eqs.1.mat, L)] = 
    hsp.up.plods.ary[rep(cis.eqs.1.mat, L)] + 0.5 *
    ale.frqs[
      cbind(gts.2[2, cis.eqs[1, ]], rep(1:L, each = sum(cis.eqs.1.mat)))
    ]
  hsp.up.plods.ary[rep(cis.eqs.2.mat, L)] = 
    hsp.up.plods.ary[rep(cis.eqs.2.mat, L)] + 0.5 *
    ale.frqs[
      cbind(gts.2[1, cis.eqs[2, ]], rep(1:L, each = sum(cis.eqs.2.mat)))
    ]
  
  # Don't add twice when second genotype homozygous
  hsp.up.plods.ary[rep(trans.eqs.gts.2.htro.1.mat, L)] = 
    hsp.up.plods.ary[rep(trans.eqs.gts.2.htro.1.mat, L)] + 0.5 *
    ale.frqs[
      cbind(
        gts.2[1, trans.eqs.gts.2.htro[1, ]], 
        rep(1:L, each = sum(trans.eqs.gts.2.htro.1.mat))
      )
    ]
  hsp.up.plods.ary[rep(trans.eqs.gts.2.htro.2.mat, L)] = 
    hsp.up.plods.ary[rep(trans.eqs.gts.2.htro.2.mat, L)] + 0.5 *
    ale.frqs[
      cbind(
        gts.2[2, trans.eqs.gts.2.htro[2, ]],
        rep(1:L, each = sum(trans.eqs.gts.2.htro.2.mat))
      )
    ]
  
  # Find conditional second genotype probabilities given parent-offspring with
  # first
  cnd.gt.2.prbs.pop.ary = hsp.up.plods.ary - pss.gt.2.prbs
  
  # Take logs of HSP vs UP pseudo likelihood ratios but skip dividing by 2 for
  # now
  hsp.up.plods.ary = log(hsp.up.plods.ary) - 
    rep(log(pss.gt.prbs), each = n.pss.gts)
  
  # Make array of first genotype probabilities
  gt.1.prbs.mat = array(pss.gt.prbs[, rep(1:L, each = n.pss.gts)], 
                         dim = c(n.pss.gts, n.pss.gts, L))
  
  # Find expected values and combine in list
  ev.sp = sum(hsp.up.plods.ary[cbind(
    rep(1:n.pss.gts, L), rep(1:n.pss.gts, L), rep(1:L, each = n.pss.gts)
  )] * pss.gt.prbs) / L - log(2)
  ev.pop = 
    sum(gt.1.prbs.mat * cnd.gt.2.prbs.pop.ary * hsp.up.plods.ary) / L - log(2)
  ev.up = sum(gt.1.prbs.mat * pss.gt.2.prbs * hsp.up.plods.ary) / L - log(2)
  ev.hsgpop = (ev.pop + ev.up) / 2
  ev.tcggpop = (ev.pop + 3 * ev.up) / 4
  ev.fcgtcgggpop = (ev.pop + 7 * ev.up) / 8
  evs = c(ev.up, ev.fcgtcgggpop, ev.tcggpop, ev.hsgpop, ev.pop, ev.sp)
  names(evs) = c("Unrelated", "First cousin", "Avuncular", "Half-sibling",
                  "Parent-offspring", "Self")
  
  # Find number of samples
  n.samps = dim(smp.gts)[3]
  
  # Find all pairs of sample indices, and number of pairs of samples
  samp.prs.inds = combn(n.samps, 2)
  samp.inds.1 = samp.prs.inds[1, ]
  samp.inds.2 = samp.prs.inds[2, ]
  n.pairs = choose(n.samps, 2)
  
  # Transform smp.gts to genotype indices
  smp.gts.new = colSums(smp.gts) + 1
  
  # Set size of batches of loci to keep memory usage < 1GB
  # n.pairs * btch.sz * 8 < 1Gb = 1e9 ~ 2^30 => btch.sz ~< 2^27 / n.pairs
  btch.sz = 2^(26 - ceiling(log(n.pairs, 2)))
  
  # Find numbers of batches
  ns.btchs = ceiling(L / btch.sz)
  
  # Make vectors for PLODs and numbers of shared markers
  hsp.up.plods = ns.shrd.mrkrs = numeric(n.pairs)
  
  # Set batch counter to zero
  btch.cnt = 0
  
  # Set timer
  s.time = proc.time()[3]
  
  # Display progress
  cat("Finding", n.pairs, "PLODs, in", sum(ns.btchs), "batches, over", L,
      "loci \n")
  cat("Batch: ")
  
  # Loop over batches of loci
  withProgress({
    for(btch.ind in 1:ns.btchs) {
      # Increment batch counter
      btch.cnt = btch.cnt + 1
      
      # Display progress
      cat(btch.cnt, "")
      
      # Find locus indices
      loci.inds = ((btch.ind - 1) * btch.sz + 1):min(btch.ind * btch.sz, L)
      
      # Lookup plods
      hsp.up.plods.obs.mat = matrix(
        hsp.up.plods.ary[
          cbind(
            as.vector(smp.gts.new[loci.inds, samp.inds.1]), 
            as.vector(smp.gts.new[loci.inds, samp.inds.2]), 
            rep(loci.inds, n.pairs)
          )
        ], 
        nrow = n.pairs, byrow = T
      )
      
      # Find and add numbers of shared markers
      ns.shrd.mrkrs = ns.shrd.mrkrs + 
        (length(loci.inds) - rowSums(is.na(hsp.up.plods.obs.mat)))
      
      # Find and add plods over shared markers
      hsp.up.plods = hsp.up.plods + rowSums(hsp.up.plods.obs.mat, na.rm = T)
      
      incProgress(amount = 1 / ns.btchs)
    }
  }, value = 0, message = "Finding HSP vs UP PLODs")
  
  
  # Divide by numbers of shared markers. Include division by two skipped for hsp
  # probs earlier. P(gts|HSP) = (P(gts|UP) + P(gts|POP)) / 2 at each locus, so
  # subtract log(2) * n.loci and divide by n.loci.
  hsp.up.plods = hsp.up.plods / ns.shrd.mrkrs - log(2)
  
  # Display progress
  cat("Done \n")
  cat("Time taken:", proc.time()[3] - s.time, "seconds \n")
  
  list(plods = hsp.up.plods, ev.up = ev.up, evs = evs)
}