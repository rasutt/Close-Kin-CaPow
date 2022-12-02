# Function to find log-genopair probabilities given that the individuals form
# one of a set of possible kinpairs.  Computes over batches of loci to limit
# memory usage to ~1Gb.
FindLogGPProbsKP = function(smp.gts, L, pss.gp.prbs.KP, wtn.prs.inds) {
  # Find number of samples
  n.samps = dim(smp.gts)[3]
  
  # Find all pairs of sample indices, and number of pairs of samples - only
  # within-survey pairs for now
  samp.prs.inds = combn(n.samps, 2)
  # samp.prs.inds = combn(n.samps, 2)[, wtn.prs.inds]
  n.pairs = choose(n.samps, 2)
  # n.pairs = sum(wtn.prs.inds)
  
  n.KP.tps = dim(pss.gp.prbs.KP[4])
  
  # Make matrix for genopair probabilities with columns for each kinship
  # considered
  lg.gp.prbs.KP = matrix(
    0, nrow = n.pairs, ncol = n.KP.tps, 
    dimnames = list(Locus = NULL, Kinship = c("UP", "HSP", "POP", "SP"))
  )
  
  # Transform smp.gts to genotype indices
  smp.gts.new = colSums(smp.gts) + 1
  
  # Set size of batches of loci to keep memory usage < 1GB
  # n.pairs * btch.sz * 8 < 1Gb = 1e9 ~ 2^30 => btch.sz ~< 2^27 / n.pairs
  btch.sz = 2^(26 - ceiling(log(n.pairs, 2)))
  
  # Find numbers of batches
  n.btchs = ceiling(L / btch.sz)
  
  # Set batch counter to zero
  btch.cnt = 0
  
  # Set timer
  s.time = proc.time()[3]
  
  # Display progress
  cat("Finding log genopair probabilities for", n.pairs, "pairs, over",
      n.KP.tps, "kinships, in", n.btchs, "batches, over", L, "loci \n")
  cat("Batch: ")
  
  # Show progress-bar
  withProgress({
    # Loop over batches of loci 
    for(btch.ind in 1:n.btchs) {
      # Increment batch counter
      btch.cnt = btch.cnt + 1
      
      # Display progress
      cat(btch.cnt, "")
      
      # Find locus indices
      loci.inds = ((btch.ind - 1) * btch.sz + 1):min(btch.ind * btch.sz, L)
      
      # Create index matrix for locus-batch
      ind.mat = cbind(
        as.vector(smp.gts.new[loci.inds, samp.prs.inds[1, ]]), 
        as.vector(smp.gts.new[loci.inds, samp.prs.inds[2, ]]), 
        rep(loci.inds, n.pairs)
      )
      
      # Look up genopair probabilities for this batch of loci and add to totals
      lg.gp.prbs.KP = lg.gp.prbs.KP + 
        apply(pss.gp.prbs.KP, 4, function(pss.gp.prbs) {
          rowSums(matrix(log(pss.gp.prbs)[ind.mat], nrow = n.pairs, byrow = T))
        })
      
      # Show progress on progress-bar
      incProgress(amount = 1 / n.btchs)
    }
  }, value = 0, message = "Finding genopair log-probabilities")

  # Display progress
  cat("Done \n")
  cat("Time taken:", proc.time()[3] - s.time, "seconds \n")
  
  # Return log genopair probabilities
  lg.gp.prbs.KP
}