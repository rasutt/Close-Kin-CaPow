# Find allele frequencies from genotypes in 2 x L x n_individuals arrays,
# representing two binary SNPs at each locus for each individual, excluding
# repeated samples of the same individual in different surveys.  Frequencies are
# returned as 2 x L matrices representing the frequencies of 0 and 1-coded SNP
# alleles at each locus
FindAleFrqs <- function(ind.gts) {
  # Frequencies of 1-coded SNP alleles are means over both alleles at each
  # locus for each sample
  ale.frqs.1 = apply(ind.gts, 2, mean)
  
  # Combine with frequencies for 0-coded alleles and return
  rbind(1 - ale.frqs.1, ale.frqs.1)
}

# Find sample genotypes, extracted from matrix of individual genotypes,
# n_samples x n_loci, rows ordered by survey-year then individual ID
FindSmpGts <- function(smp.hsts, ind.gts) {
  # Sample-individual indices, row numbers in sample history matrix,
  # representing the individual that each sample came from, ordered by
  # survey-year then individual ID
  smp.indvdl.inds = row(smp.hsts)[as.logical(smp.hsts)]
  
  ind.gts[, , smp.indvdl.inds]
}

# Function to find offset sample index pairs, n_pairs x 2, consecutive pairs
# within surveys, and between all pairs of surveys, with samples repeated for
# survey with fewer samples
FindOSIPs = function(k, syis) {
  # Sample index pairs, 0 x 2, to concatenate pairs to
  sips = matrix(NA, 0, 2)
  
  # Possible survey-year indices
  psyis = 0:(k - 1)
  
  # Loop over first survey-year indices
  for (i in psyis) {
    # Sample indices for survey year one of pair, in random order
    sissy.1 = which(syis == i)
    n.sissy.1 = length(sissy.1)
    sissy.1 = sissy.1[sample(n.sissy.1)]
    
    # If at least one sample
    if (n.sissy.1 > 0) {
      # Add offset index pairs for samples in this year
      sips = rbind(sips, cbind(sissy.1[-n.sissy.1], sissy.1[-1]))
      
      # Loop over second survey-year indices
      for (j in psyis[psyis > i]) {
        # Sample indices for survey year two of pair, in random order
        sissy.2 = which(syis == j)
        n.sissy.2 = length(sissy.2)
        sissy.2 = sissy.2[sample(n.sissy.2)]
        
        # Check at least one sample
        if (n.sissy.2 > 0) {
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

# Function to find genopair probabilities by exponentiating log-probabilities,
# checking for underflow, and trying to adjust if necessary
FindGPs = function(glps) {
  # Get genopair probabilities (by excluding probabilities giveb half-sibs for
  # now) and check for pairs where all probabilities underflow to zero
  gps = exp(glps)
  colnames(gps) = colnames(glps)
  all_undrflw = rowSums(gps) == 0
  
  # If there is underflow adjust log-probabilities by factor giving equal
  # weight to smallest and largest values to avoid both under and overflow
  if (any(all_undrflw)) {
    cat("Proportion of pairs for which probabilities given all kinships 
        underflow to zero:", mean(all_undrflw), "\n")
    
    # Want smallest maximum kinship probability and largest probability to be
    # equally far from one
    adj = mean(c(min(apply(glps, 1, max)), max(glps)))
    glps.adj = glps - adj
    gps.adj = exp(glps.adj)
    colnames(gps.adj) = colnames(glps)
    
    # Show adjustment and results
    cat("Probabilities adjusted by factor of exp(", adj, ")\n", sep = "")
    print("Adjusted log-probabilities and probabilities of genopairs:")
    print(summary(glps.adj))
    print(summary(gps.adj))
    cat("Proportion of pairs for which adjusted probabilities given all 
        kinships underflow to zero:", mean(rowSums(gps.adj) == 0), "\n")
    
    return(gps.adj)
  } 
  else {
    return(gps)
  }
}

# Function to find genopair probabilities from genotypes and sample-individual
# index pairs
FindGPsMdl <- function(pop.cap.hist, L, knshp.st, siips) {
  # Get individual genotypes
  gts = attributes(pop.cap.hist)$ind.gts
    
  # Allele frequencies, 2 x n_loci matrices, representing relative
  # frequencies of 0 and 1-coded SNP alleles at each locus
  ale.frqs = FindAleFrqs(gts)
  
  # Possible genopair probabilities given kinship set, n_possible_genotypes
  # x n_possible_genotypes x n_loci x n_kinships, representing probabilities
  # of each possible pair of genotypes at each locus given each kinship
  # considered
  pgps = FindPssGPPsKPs(ale.frqs, L, knshp.st)
  
  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  glps = FindGLPs(pgps, gts, siips, L)
  
  # Exponentiate genopair log-probabilities given kinship set, checking for
  # underflow and trying to adjust if necessary
  FindGPs(glps)
}