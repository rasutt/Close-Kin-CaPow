# Function to find log-genopair probabilities given that the individuals form
# one of a set of possible kinpairs.  Computes over batches of loci to limit
# memory usage to ~1Gb.
FindLogGPProbsKP = function(
    pss.gp.prbs.KP, smp.gts, smp.ind.prs, L, sngl.knshp = F
) {
  # Indices of first and second samples in each pair
  smp.1.inds = smp.ind.prs[1, ]
  smp.2.inds = smp.ind.prs[2, ]
  
  n.pairs = length(smp.1.inds)
  if (sngl.knshp) {
    n.KP.tps = 1
  } else {
    n.KP.tps = dim(pss.gp.prbs.KP)[4]
  }
  pss.gp.lg.prbs.KP = log(pss.gp.prbs.KP)
  
  # Make matrix for genopair probabilities with columns for each kinship
  # considered
  gp.lg.prbs.KP = matrix(
    0, nrow = n.pairs, ncol = n.KP.tps,
    dimnames = list(pair = 1:n.pairs, Kinship = dimnames(pss.gp.prbs.KP)[[4]])
  )
  
  # Transform genotypes from 2 x L x n_samples arrays of alleles to n_samples x
  # L matrices of indices of genotypes ordered 00, 01, and 11, for binary SNPs
  smp.gt.inds = t(colSums(smp.gts)) + 1
  
  # Set size of batches of loci to keep memory usage < 1GB, n.pairs * btch.sz *
  # 8 bytes < 1Gb = 1e9 bytes ~ 2^30 => btch.sz ~< 2^27 / n.pairs
  btch.sz = 2^(26 - ceiling(log(n.pairs, 2)))
  
  # Find numbers of batches
  n.btchs = ceiling(L / btch.sz)
  
  # Set batch counter to zero
  btch.cnt = 0
  
  # Set timer
  s.time = proc.time()[3]
  
  # Loop over batches of loci 
  for(btch.ind in 1:n.btchs) {
    # Increment batch counter
    btch.cnt = btch.cnt + 1
    
    # Find indices for current batch of loci
    loci.inds = ((btch.ind - 1) * btch.sz + 1):min(btch.ind * btch.sz, L)
    n.loci.btch = length(loci.inds)
    
    # Indices of genotypes for current batch of loci, n_samples x n_loci_batch
    smp.gt.inds.btch = smp.gt.inds[, loci.inds]
    
    # Indices of genotypes for first and second samples in each pair for
    # current batch of loci, n_pairs x n_loci_batch
    smp.1.gt.inds.btch = smp.gt.inds.btch[smp.1.inds, ]
    smp.2.gt.inds.btch = smp.gt.inds.btch[smp.2.inds, ]
    
    # Matrix of indices of genotypes for first and second samples in each pair
    # for current batch of loci, (n_pairs x n_loci_batch) x 3
    ind.mat = cbind(
      as.vector(smp.1.gt.inds.btch),
      as.vector(smp.2.gt.inds.btch),
      rep(loci.inds, each = n.pairs)
    )
    
    # Look up genopair log-probabilities for this batch of loci and add to
    # totals. Possible values are stored in n_gts x n_gts x L x kinships
    # array. Using apply to pick out one kinship at a time to limit memory
    # usage. No significant slow down as just a few kinships.
    if (sngl.knshp) {
      gp.lg.prbs.KP = gp.lg.prbs.KP + 
        rowSums(matrix(pss.gp.lg.prbs.KP[ind.mat], nrow = n.pairs))
    } else {
      gp.lg.prbs.KP = gp.lg.prbs.KP + 
        apply(pss.gp.lg.prbs.KP, 4, function(lg.pss.gp.prbs) {
          rowSums(matrix(lg.pss.gp.prbs[ind.mat], nrow = n.pairs))
        })
    }
  }
  
  # Return log genopair probabilities
  gp.lg.prbs.KP
}
