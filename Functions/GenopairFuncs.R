# Find allele frequencies from genotypes in 2 x L x n_individuals arrays,
# representing two binary SNPs at each locus for each individual, excluding
# repeated samples of the same individual in different surveys.  Frequencies are
# returned as 2 x L matrices representing the frequencies of 0 and 1-coded SNP
# alleles at each locus
FindAleFrqs <- function(unq.smp.gts) {
  # Frequencies of 1-coded SNP alleles are means over both alleles at each
  # locus for each sample
  ale.frqs.1 = apply(unq.smp.gts, 2, mean)
  
  # Combine with frequencies for 0-coded alleles and return
  rbind(1 - ale.frqs.1, ale.frqs.1)
}

# Find sample genotypes, extracted from matrix of individual genotypes,
# n_samples x n_loci, rows ordered by survey-year then individual ID
FindSmpGts <- function(smp.hsts, unq.smp.gts) {
  # Sample indices, rows in sample history matrix, representing the
  # individual that each sample came from, ordered by survey-year then
  # individual ID
  smp.inds = row(smp.hsts)[as.logical(smp.hsts)]
  
  unq.smp.gts[, , smp.inds]
}

# Sample index pairs, 2 x n_pairs matrix of indices of samples in each pair to
# include in likelihood, possibly all pairs or just consecutive pairs with
# individuals shuffled out of ID/birth order
FindSIPsOffset = function(k, smp.yr.inds) {
  # Empty matrix to add sample index pairs to, 2 x 0
  smp.ind.prs = matrix(0, 2, 0)
  
  # Possible survey-year indices
  pss.s.yr.inds = 0:(k - 1)
  
  # Loop over first survey-year indices
  for (i in pss.s.yr.inds) {
    # Get indices for samples in this year and randomize order
    smp.inds.yr.1 = which(smp.yr.inds == i)
    smp.inds.yr.1 = smp.inds.yr.1[sample(length(smp.inds.yr.1))]
    
    # If at least one sample
    if (length(smp.inds.yr.1) > 0) {
      # Add offset index pairs for samples in this year
      smp.ind.prs = cbind(
        smp.ind.prs, 
        rbind(smp.inds.yr.1, c(smp.inds.yr.1[-1], smp.inds.yr.1[1]))
      )
      
      # Loop over second survey-year indices
      for (j in pss.s.yr.inds[pss.s.yr.inds > i]) {
        # Get indices for samples in this year 
        smp.inds.yr.2 = which(smp.yr.inds == j)
        smp.inds.yr.2 = smp.inds.yr.2[sample(length(smp.inds.yr.2))]
        
        # Check at least one sample
        if (length(smp.inds.yr.2) > 0) {
          # Add offset index pairs for samples from these years, shorter index
          # set recycled with warning
          smp.ind.prs = suppressWarnings(
            cbind(smp.ind.prs, rbind(smp.inds.yr.1, smp.inds.yr.2))
          )
        }
      }
    }
  }
  
  smp.ind.prs
}

FindGPPs = function(LGPPs) {
  # Get genopair probabilities (by excluding probabilities giveb half-sibs for
  # now) and check for pairs where all probabilities underflow to zero
  gpp.slct = exp(LGPPs)
  all_undrflw = rowSums(gpp.slct) == 0
  
  # If there is underflow adjust log-probabilities by factor giving equal
  # weight to smallest and largest values to avoid both under and overflow
  if (any(all_undrflw)) {
    cat("Proportion of pairs for which probabilities given all kinships 
        underflow to zero:", mean(all_undrflw), "\n")
    
    # Want smallest maximum kinship probability and largest probability to be
    # equally far from one
    adj = mean(c(min(apply(LGPPs, 1, max)), max(LGPPs)))
    lg.gpp.adj = LGPPs - adj
    gpp.adj = exp(lg.gpp.adj)
    
    # Show adjustment and results
    cat("Probabilities adjusted by factor of exp(", adj, ")\n", sep = "")
    print("Adjusted log-probabilities and probabilities of genopairs:")
    print(summary(lg.gpp.adj))
    print(summary(gpp.adj))
    cat("Proportion of pairs for which adjusted probabilities given all 
        kinships underflow to zero:", mean(rowSums(gpp.adj) == 0), "\n")
    
    return(gpp.adj)
  } 
  else {
    return(gpp.slct)
  }
}

FindGpMdlInpts <- function(pop.cap.hist, L, k, os.mdl, knshp.st) {
  # Sample history matrix, n_individuals x n_surveys, rows ordered by
  # individual ID
  smp.hsts = as.matrix(pop.cap.hist[, 4:(3 + k)])
  
  # Sample genotypes, extracted from matrix of individual genotypes,
  # n_samples x n_loci, rows ordered by survey-year then individual ID
  smp.gts = FindSmpGts(smp.hsts, attributes(pop.cap.hist)$unq.smp.gts)
  
  # Sample-year indices, columns in sample history matrix, representing the
  # survey that each sample came from, ordered by survey-year then
  # individual ID.  Counting from zero as will be passed to TMB objective
  # function
  smp.yr.inds = col(smp.hsts)[as.logical(smp.hsts)] - 1
  
  # Sample index pairs, 2 x n_pairs, representing indices of samples in each
  # pair to include in likelihood, possibly all pairs or just consecutive
  # pairs
  smp.ind.prs = 
    if (os.mdl) {
      FindSIPsOffset(k, smp.yr.inds)
    } else {
      combn(length(smp.yr.inds), 2)
    }
  
  # Sample-year index pairs, n_pairs x 2, representing survey-year of each
  # sample in each pair, counting from zero as passing into TMB C++
  # objective function
  smp.yr.ind.prs = matrix(smp.yr.inds[as.vector(t(smp.ind.prs))], ncol = 2)
  
  # Allele frequencies, 2 x n_loci matrices, representing relative
  # frequencies of 0 and 1-coded SNP alleles at each locus
  ale.frqs = FindAleFrqs(attributes(pop.cap.hist)$unq.smp.gts)
  
  # Possible genopair probabilities given kinship set, n_possible_genotypes
  # x n_possible_genotypes x n_loci x n_kinships, representing probabilities
  # of each possible pair of genotypes at each locus given each kinship
  # considered
  pss.gp.prbs.KPs = FindPssGPPsKPs(ale.frqs, L, knshp.st)
  
  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  lg.gp.prbs.KPs = 
    FindLogGPProbsKP(pss.gp.prbs.KPs, smp.gts, smp.ind.prs, L)
  
  # Exponentiate genopair log-probabilities given kinship set, checking for
  # underflow and trying to adjust if necessary.  Would be good to raise a
  # proper error if adjustment impossible
  GPPs = FindGPPs(lg.gp.prbs.KPs)
  
  list(GPPs = GPPs, smp.yr.ind.prs = smp.yr.ind.prs)
}