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
  # Sample-individual indices, row numbers in sample history matrix,
  # representing the individual that each sample came from, ordered by
  # survey-year then individual ID
  smp.indvdl.inds = row(smp.hsts)[as.logical(smp.hsts)]
  
  unq.smp.gts[, , smp.indvdl.inds]
}

# Function to find sample index pairs for offset model, n_pairs x 2,
# including consecutive pairs within surveys, and between all pairs of surveys,
# with samples repeated for survey with fewer samples
FindSIPsOffset = function(k, syis) {
  # Sample index pairs, 0 x 2, to concatenate pairs to
  sips = matrix(NA, 0, 2)
  
  # Possible survey-year indices
  psyis = 0:(k - 1)
  
  # Loop over first survey-year indices
  for (i in psyis) {
    # Sample indices for survey year one of pair, in random order
    sissy.1 = which(syis == i)
    sissy.1 = sissy.1[sample(length(sissy.1))]
    
    # If at least one sample
    if (length(sissy.1) > 0) {
      # Add offset index pairs for samples in this year
      sips = rbind(sips, cbind(sissy.1, c(sissy.1[-1], sissy.1[1])))
      
      # Loop over second survey-year indices
      for (j in psyis[psyis > i]) {
        # Sample indices for survey year two of pair, in random order
        sissy.2 = which(syis == j)
        sissy.2 = sissy.2[sample(length(sissy.2))]
        
        # Check at least one sample
        if (length(sissy.2) > 0) {
          # Add offset index pairs for samples from these years, shorter index
          # set recycled with warning
          sips = suppressWarnings(rbind(sips, cbind(sissy.1, sissy.2)))
        }
      }
    }
  }
  
  # Return
  sips
}

FindGPsGvnKs = function(LGPPs) {
  # Get genopair probabilities (by excluding probabilities giveb half-sibs for
  # now) and check for pairs where all probabilities underflow to zero
  gpp.slct = exp(LGPPs)
  colnames(gpp.slct) = colnames(LGPPs)
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
    colnames(gpp.adj) = colnames(LGPPs)
    
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