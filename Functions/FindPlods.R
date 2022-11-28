# Function to find PLODS for SNP genotypes.  Simplified from PLOD code for Emma.
FindPlods = function(
    smp.gts, L, pss.plods.pls.lg.2
) {
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
        pss.plods.pls.lg.2[
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
  
  list(plods = hsp.up.plods)
}