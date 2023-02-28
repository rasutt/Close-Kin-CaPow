# Code to find full set of genopair log-probabilities

# Sample-individual and sample-year indices, length n_samples, starting from
# zero for sample-years for TMB functions, for all studies
sisyis = reactive({
  # Make lists
  siis = syis = vector("list", n.sims())
  
  # Show progress-bar
  withProgress(
    # Loop over histories
    for (hst.ind in 1:n.sims()) {
      # Get simulated family and capture histories of population of
      # animals over time
      stdy <- sim.lst()$hists.lst[[hst.ind]]
      
      # Sample history matrix, n_individuals x n_surveys, rows ordered by
      # individual ID, using == 1 as simpler for data frame input
      smp.hsts.bln = stdy[, 4:(3 + k())] == 1
      
      # Sample-individual indices
      siis[[hst.ind]] = row(smp.hsts.bln)[smp.hsts.bln]
      
      # Sample-year indices, starting from zero for TMB functions
      syis[[hst.ind]] = col(smp.hsts.bln)[smp.hsts.bln] - 1
      
      # Update progress-bar
      incProgress(1/n.sims())
    }, 
    value = 0, 
    message = "Finding all sample-individual and sample-year indices"
  )
  
  # Return list
  list(siis = siis, syis = syis)
})

# Full set of sample-individual and sample-year index pairs for genopair
# models, 2 x n_pairs x 2, representing individual and year that each
# sample came from. Sample-years start from zero for TMB C++ objective
# function
fsisyips = reactive({
  # Make lists
  siips = syips = vector("list", n.sims())
  
  # Show progress-bar
  withProgress(
    # Loop over histories
    for (hst.ind in 1:n.sims()) {
      # Sample-individual and sample-year index pairs
      siips[[hst.ind]] = t(combn(sisyis()$siis[[hst.ind]], 2))
      syips[[hst.ind]] = t(combn(sisyis()$syis[[hst.ind]], 2))
      
      # Update progress-bar
      incProgress(1/n.sims())
    }, 
    value = 0, 
    message = "Finding all sample-individual and sample-year index pairs"
  )
  
  # Return lists
  list(siips = siips, syips = syips)
})

# Allele frequencies for each study, 2 x n_loci X n_sims
ale.frqs.ary = reactive({
  # Loop over histories
  withProgress({
    sapply(
      sim.lst()$hists.lst, 
      function(stdy) FindAleFrqs(attributes(stdy)$ind.gts),
      simplify = "array"
    )
  }, value = 0, message = "Finding allele frequencies")
})

# Possible genotype probabilities for each study, 3 x n_loci x n_sims
pgtps.ary = reactive(FindPssGtPrbsAry(ale.frqs.ary()))

# Possible first genotype probabilities for genopairs in each study, 3 x 3 x
# n_loci x n_sims, first genotypes are rows
pgt1ps.ary = reactive(FindPssFrstGtPrbsAry(pgtps.ary(), L(), n.sims()))

# Possible genopair probabilities for each kinship for each study, 3 x 3 x
# n_loci x n_sims
pgps.UPs.ary = reactive(FindPssGpPsUPsAry(pgt1ps.ary()))
pgps.SPs.ary = reactive(FindPssGpPsSPsAry(pgtps.ary(), L(), n.sims()))
pgps.POPs.ary = reactive({
  FindPssGpPsPOPsAry(pgt1ps.ary(), ale.frqs.ary(), L(), n.sims())
})
pgps.HSPs.ary = reactive(FindPssGpPsHSPsAry(pgps.UPs.ary(), pgps.POPs.ary()))

## Find full set of genopair log-probabilities for each kinship

# Function for normal serial computing, taking possible genopair probabilities
# given a certain kinship, over all studies, and sample-individual index pairs
# as inputs, and able to see reactive objects
FindGLPsSrl = function(hst.ind, pgps.ary, siips.lst) {
  # Update progress bar
  incProgress(1/n.sims())

  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  FindGLPs(
    pgps.ary[, , , hst.ind], 
    attributes(sim.lst()$hists.lst[[hst.ind]])$ind.gts,
    siips.lst[[hst.ind]], L(), sk = T
  )
}

# Function using "parallel" library to use multiple CPU cores, giving > 2x
# speedup for datasets >= 100 loci and 100 studies. Sample histories, genotypes,
# and possible genopair probabilities are exported to the new R sessions in
# advance
FindFGLPsPl = function(hst.ind) {
  # Genopair log-probabilities over all loci given each kinship, for each
  # pair to include in likelihood
  FindGLPs(
    pgps.ary[, , , hst.ind], 
    attributes(hsts.pl[[hst.ind]])$ind.gts, 
    siips.pl[[hst.ind]], L.pl, sk = T
  )
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
          lst = parLapply(cl(), 1:n.sims(), FindFGLPsPl)
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
        L(), "loci and", n.sims(), "studies in", Sys.time() - s, "\n\n"
      )

      # Return results
      lst
    }

    # Otherwise use normal code
    else {
      # Find actual genopair log-probabilities
      withProgress(
        {
          lst = lapply(1:n.sims(), FindGLPsSrl, pgps.ary, fsisyips()$siips)
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
fglps.UPs = MakeFGLPsRctv("unrelated")
fglps.SPs = MakeFGLPsRctv("self")
fglps.POPs = MakeFGLPsRctv("parent-offspring")
fglps.HSPs = MakeFGLPsRctv("half-sibling")
