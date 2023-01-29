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
pgps.UPs.ary = reactive({
  FindPssGpPsUPsAry(pss.gt.1.prbs.ary())
})
pgps.SPs.ary = reactive({
  FindPssGpPsSPsAry(pss.gt.prbs.ary(), L(), n.sims())
})
pgps.POPs.ary = reactive({
  FindPssGpPsPOPsAry(pss.gt.1.prbs.ary(), ale.frqs.ary(), L(), n.sims())
})
pgps.HSPs.ary = reactive({
  FindPssGpPsHSPsAry(pgps.UPs.ary(), pgps.POPs.ary())
})

# Find full set of genopair log-probabilities for each kinship

# Normal function, taking an array of possible genopair probabilities given a
# certain kinship, over all studies, as an input, and using reactive objects
FindFLGPs = function(hst.ind, pgps.ary) {
  # Update progress bar
  incProgress(1/n.sims())
  
  # Get simulated family and capture histories of population of animals
  # over time
  pop.cap.hist <- sim.lst()$hists.lst[[hst.ind]]
  
  # Individual genotypes, n_individuals x n_loci, representing genotypes for
  # each individual sampled
  gts = attributes(pop.cap.hist)$unq.smp.gts
  
  # Possible genopair probabilities for this study, n_genotypes x n_genotypes x
  # n_loci
  pgps = pgps.ary[, , , hst.ind]
  
  # Sample-individual index pairs, n_pairs x 2, representing the individual that
  # each sample in each pair came from
  siips = fsisyips.lst()$siips.lst[[hst.ind]]
  
  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  FindGLPs(pgps, gts, siips, L(), sk = T)
}

# Function using "parallel" library to use multiple CPU cores, giving > 2x
# speedup for datasets >= 100 loci and 100 studies. Sample histories, genotypes,
# and possible genopair probabilities are exported to the new R sessions in
# advance
FindFLGPsPrll = function(hst.ind) {
  # Get simulated family and capture histories of population of animals
  # over time
  pop.cap.hist <- hst.lst.prll[[hst.ind]]
  
  # Individual genotypes, n_individuals x n_loci, representing genotypes for
  # each individual sampled
  gts = attributes(pop.cap.hist)$unq.smp.gts
  
  # Possible genopair probabilities for this study, n_genotypes x n_genotypes x
  # n_loci
  pgps = pgps.ary[, , , hst.ind]
  
  # Sample-individual index pairs, n_pairs x 2, representing individual that
  # each sample in each pair came from
  siips = siips.lst.prll[[hst.ind]]
  
  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  FindGLPs(pgps, gts, siips, L.prll, sk = T)
}

# Make full genopair log-probabilities reactive function given a particular
# kinship
MakeFGLPsRctv = function(knshp) {
  reactive({
    # Possible genopair probabilities given selected kinship, n_genotypes x
    # n_genotypes x n_loci
    pgps.ary = switch(
      knshp,
      "unrelated" = pgps.UPs.ary(),
      "self" = pgps.SPs.ary(),
      "parent-offspring" = pgps.POPs.ary(),
      "half-sibling" = pgps.HSPs.ary()
    )
    
    # If the numbers of studies and/or loci are large use multicore processing
    if (n.sims() * L() > 1e3) {
      # Record time 
      s = Sys.time()
      
      # Send possible genopair probabilities to R sessions on multiple CPU cores
      clusterExport(cl(), list("pgps.ary"), environment())
      
      # Find actual genopair log-probabilities for different chunks of studies
      # on different nodes and combine results when done
      withProgress(
        {
          lst = parLapply(cl(), 1:n.sims(), FindFLGPsPrll)
        },
        value = 0,
        message = paste(
          "Finding genopair log-probabilities given", knshp, "pairs"
        ),
        detail = "Using multiple cores so no progress updates available"
      )
      
      # Show time taken
      cat(
        "Found log probabilities, given", knshp, "pairs, over",
        L(), "loci and",
        n.sims(), "studies in",
        Sys.time() - s, "\n\n"
      )
      
      # Return results
      lst
    } 
    
    # Otherwise use normal code
    else {
      # Find actual genopair log-probabilities
      withProgress(
        {
          lst = lapply(1:n.sims(), FindFLGPs, pgps.ary)
        },
        value = 0,
        message = paste(
          "Finding genopair log-probabilities given", knshp, "pairs"
        )
      )
      
      # Return results
      lst
    }
  })
}
fglps.UPs.lst = MakeFGLPsRctv("unrelated")
fglps.SPs.lst = MakeFGLPsRctv("self")
fglps.POPs.lst = MakeFGLPsRctv("parent-offspring")
fglps.HSPs.lst = MakeFGLPsRctv("half-sibling")
