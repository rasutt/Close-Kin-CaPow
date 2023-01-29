# Function to find genopair log-probabilities given that the individuals have a
# particular kinship.  Computes over batches of loci to limit memory usage to
# ~1Gb.
FindGLPs = function(pgps, gts, siips, L, sk = F) {
  # Sample-individual indices of first and second samples in each pair, n_pairs
  siis.1 = siips[, 1]
  siis.2 = siips[, 2]
  
  # Number of pairs of samples
  n.pairs = length(siis.1)
  
  # Number of kinships for which genopair log-probabilities are being found
  n.knshps = if (sk) 1 else dim(pgps)[4]
  
  # Possible genopair log-probabilities at each locus given kinships,
  # n_genotypes x n_genotypes x n_loci x n_kinships
  pglps = log(pgps)
  
  # Actual genopair log-probabilities given kinships, n_pairs x n_kinships
  glps = matrix(
    0, nrow = n.pairs, ncol = n.knshps,
    dimnames = list(pair = 1:n.pairs, Kinship = dimnames(pgps)[[4]])
  )
  
  # Genotype indices at each locus, n_individuals x n_loci, for genotypes
  # ordered 00, 01, 11, computed from binary alleles at each locus, 2 x n_loci x
  # n_individuals
  gt.inds = t(colSums(gts)) + 1
  
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
    
    # Genotype indices for current batch of loci, n_individuals x n_loci_batch
    gt.inds.btch = gt.inds[, loci.inds]
    
    # Sample genotype indices for current batch of loci for first and second
    # samples in each pair, n_pairs x n_loci_batch
    sgisb.1 = gt.inds.btch[siis.1, ]
    sgisb.2 = gt.inds.btch[siis.2, ]
    
    # Sample genotype and locus-indices for current batch of loci for first and
    # second samples in each pair, (n_pairs x n_loci_batch) x 3
    sglisb = cbind(
      as.vector(sgisb.1),
      as.vector(sgisb.2),
      rep(loci.inds, each = n.pairs)
    )
    
    # Look up genopair log-probabilities for this batch of loci and add to
    # totals. Possible values are stored in n_gts x n_gts x L x kinships
    # array. Using apply to pick out one kinship at a time to limit memory
    # usage. No significant slow down as just a few kinships.
    if (sk) {
      glps = glps + rowSums(matrix(pglps[sglisb], nrow = n.pairs))
    } else {
      glps = glps + 
        apply(pglps, 4, function(pglps.knshp) {
          rowSums(matrix(pglps.knshp[sglisb], nrow = n.pairs))
        })
    }
  }
  
  # Return log genopair probabilities
  glps
}
