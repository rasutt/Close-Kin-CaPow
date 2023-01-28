# Allele frequencies for each study, 2 x n_loci X n_sims
ale.frqs.ary = reactive({
  # Make empty list
  ale.frqs.ary.tmp = array(dim = c(2, L(), n.sims()))
  
  # Loop over histories
  withProgress({
    for (hst.ind in 1:n.sims()) {
      # Get simulated family and capture histories of population of animals
      # over time
      pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
      
      # Allele frequencies, 2 x n_loci matrices, representing relative
      # frequencies of 0 and 1-coded SNP alleles at each locus
      ale.frqs.ary.tmp[, , hst.ind] = 
        FindAleFrqs(attributes(pop.cap.hist)$unq.smp.gts)
      
    }
  }, value = 0, message = "Finding allele frequencies")
  
  ale.frqs.ary.tmp
})

# Possible genotype probabilities for each study, 3 x n_loci x n_sims
pss.gt.prbs.ary = reactive({
  FindPssGtPrbsAry(ale.frqs.ary())
})

# Possible first genotype probabilities for sample pairs for each study, 3 x 3 x
# n_loci x n_sims, first genotypes are rows
pss.gt.1.prbs.ary = reactive({
  FindPssFrstGtPrbsAry(pss.gt.prbs.ary(), L(), n.sims())
})

# Possible genopair probabilities for each kinship for each study, 3 x 3 x
# n_loci x n_sims
pss.gp.prbs.UPs.ary = reactive({
  FindPssGpPsUPsAry(pss.gt.1.prbs.ary())
})
pss.gp.prbs.SPs.ary = reactive({
  FindPssGpPsSPsAry(pss.gt.prbs.ary(), L(), n.sims())
})
pss.gp.prbs.POPs.ary = reactive({
  FindPssGpPsPOPsAry(pss.gt.1.prbs.ary(), ale.frqs.ary(), L(), n.sims())
})
pss.gp.prbs.HSPs.ary = reactive({
  FindPssGpPsHSPsAry(pss.gp.prbs.UPs.ary(), pss.gp.prbs.POPs.ary())
})

# Find full set of genopair log-probabilities for each kinship

# Normal function, taking an array of possible genopair probabilities given a
# certain kinship, over all studies, as an input, and using reactive objects
FindFllLgGpPrbsKP = function(hst.ind, pss.gp.prbs.KP.ary, kp.tp) {
  # Get simulated family and capture histories of population of animals
  # over time
  pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
  
  # Sample history matrix, n_individuals x n_surveys, rows ordered by
  # individual ID
  smp.hsts = as.matrix(pop.cap.hist[, 4:(3 + k())])
  
  # Sample index pairs, 2 x n_pairs, representing indices of samples in each
  # pair
  smp.ind.prs = combn(sum(smp.hsts), 2)
  
  # Sample genotypes, extracted from matrix of individual genotypes,
  # n_samples x n_loci, rows ordered by survey-year then individual ID
  smp.gts = FindSmpGts(smp.hsts, attributes(pop.cap.hist)$unq.smp.gts)
  
  # Possible genopair probabilities over genopairs and loci for this study
  pss.gp.prbs.KP = pss.gp.prbs.KP.ary[, , , hst.ind]
  
  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  FindLogGPProbsKP(
    pss.gp.prbs.KP, smp.gts, smp.ind.prs, L(), sngl.knshp = T, kp.tp
  )
}

# Function using "parallel" library to use multiple CPU cores, giving > 2x
# speedup for datasets >= 100 loci and 100 studies. Sample histories, genotypes,
# and possible genopair probabilities are exported to the new R sessions in
# advance
FindFllLgGpPrbsKPPrll = function(hst.ind) {
  # Get simulated family and capture histories of population of animals
  # over time
  pop.cap.hist <- hst.lst.prll[[hst.ind]]
  
  # Individual genotypes, n_individuals x n_loci, representing genotypes for
  # each individual sampled
  indvdl.gts = attributes(pop.cap.hist)$unq.smp.gts
  
  # Possible genopair probabilities over genopairs and loci for this study
  pss.gp.prbs.KP = pss.gp.prbs.KP.ary[, , , hst.ind]
  
  # Sample-individual index pairs, n_pairs x 2, representing individual that
  # each sample in each pair came from
  SIIPs = SIIPs.lst.prll[[hst.ind]]
  
  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  FindLogGPProbsKP(pss.gp.prbs.KP, indvdl.gts, SIIPs, L.prll, sngl.knshp = T)
}

# Set up multicore processing and find log-probabilities given unrelated pairs 
MakeFLGPsKPRctv = function(kp.tp) {
  reactive({
    # Record time
    s = Sys.time()
    
    pss.gp.prbs.KP.ary = switch(
      kp.tp,
      "unrelated" = pss.gp.prbs.UPs.ary(),
      "self" = pss.gp.prbs.SPs.ary(),
      "parent-offspring" = pss.gp.prbs.POPs.ary(),
      "half-sibling" = pss.gp.prbs.HSPs.ary()
    )
    clusterExport(cl(), list("pss.gp.prbs.KP.ary"), environment())
    
    # Find log-probabilities for different chunks of studies on different nodes
    # and combine results when done
    withProgress(
      {
        lst = parLapply(cl(), 1:n.sims(), function(hst.ind) {
          FindFllLgGpPrbsKPPrll(hst.ind)
        })
      },
      value = 0,
      message = paste(
        "Finding genopair log-probabilities given", kp.tp, "pairs"
      ),
      detail = "Using multiple cores so no progress updates available"
    )
    
    # Show time taken
    cat(
      "Found log probabilities, given", kp.tp, "pairs, over",
      L(), "loci and",
      n.sims(), "studies in",
      Sys.time() - s, "\n\n"
    )
    
    # Return results
    lst
  })
}
FGLPs.UPs.lst = MakeFLGPsKPRctv("unrelated")
FGLPs.SPs.lst = MakeFLGPsKPRctv("self")
FGLPs.POPs.lst = MakeFLGPsKPRctv("parent-offspring")
FGLPs.HSPs.lst = MakeFLGPsKPRctv("half-sibling")
