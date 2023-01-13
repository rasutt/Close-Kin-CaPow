# Function to find possible genotype probabilities
FindPssGtPrbs = function(ale.frqs) {
  ales.1.inds = pss.gts[1, ]
  ales.2.inds = pss.gts[2, ]
  
  # Indexing the 2 x L allele frequencies matrix for each allele of each
  # possible genotype (globally defined for SNP genotypes), and multiplying by 2
  # possible cases for heterozygous genotypes
  matrix(
    ale.frqs[ales.1.inds, ] * ale.frqs[ales.2.inds, ] * 
      (1 + (ales.1.inds != ales.2.inds)), 
    nrow = n.pss.gts
  )
}

# Function to find possible first genotype probabilities for genopairs
FindPssFrstGtPrbs = function(pss.gt.prbs, L) {
  # 3 x 1 x L array which is indexed 3 times to fill the three columns.
  array(
    pss.gt.prbs[, rep(1:L, each = n.pss.gts)], c(n.pss.gts, n.pss.gts, L)
  )
}

# Function to find possible genopair probabilities for unrelated pairs
FindPssGpPsUPs = function(pss.gt.1.prbs) {
  # Conditional probabilities given that the pair are unrelated, the products of
  # the respective genotype probabilities, the second found by permuting the
  # array containing the first.
  pss.gt.1.prbs * aperm(pss.gt.1.prbs, c(2, 1, 3))
}

# Function to find possible genopair probabilities for parent-offspring pairs
FindPssGpPsPOPs = function(ale.frqs, L, pss.gt.1.prbs) {
  # Conditional probabilities given that the pair are parent and offspring
  # (unordered). Products of the first genotype probabilities and the
  # conditional probabilities of the second genotypes given that the pair are
  # parent and offspring.
  pss.gt.1.prbs * 
    
    # Conditional probabilities of second genotype given that the pair are parent
    # and offspring (unordered).  0.5 for each allele in first genotype being
    # inherited, multiplied by the probability of the second genotype in each
    # case, as in table 3, pg. 269, Bravingtion et al. (2016) Close-Kin
    # Mark-Recapture. Filled into an L x 3 x 3 array which is permuted to the
    # standard dimensions.
    aperm(
      array(
        # Note: order data enters array is down columns, not across rows
        c(
          ale.frqs[1, ], 0.5 * ale.frqs[1, ], rep(0, L), 
          ale.frqs[2, ], 0.5 * colSums(ale.frqs), ale.frqs[1, ],
          rep(0, L), 0.5 * ale.frqs[2, ], ale.frqs[2, ]
        ), 
        c(L, n.pss.gts, n.pss.gts)
      ), 
      c(2, 3, 1)
    )
}

# Function to find possible genopair probabilities for self-pairs
FindPssGpPsSPs = function(pss.gt.1.prbs, L) {
  # Conditional probabilities given that the pair are the same individual
  # sampled twice. Genotype probabilities when the genotypes are the same, and
  # zero otherwise.
  aperm(
    array(
      cbind(
        pss.gt.1.prbs[1, 1, ], 0, 0, 
        0, pss.gt.1.prbs[2, 2, ], 0, 
        0, 0, pss.gt.1.prbs[3, 3, ]
      ), 
      c(L, n.pss.gts, n.pss.gts)
    ), 
    c(2, 3, 1)
  )
}

# Function to find possible genopair probabilities for half-sibling pairs
FindPssGpPsHSPs = function(pss.gp.prbs.UP, pss.gp.prbs.POP) {
  # Conditional probabilities given that the pair are half-siblings.  Average of
  # probabilities for unrelated and parent-offspring pairs.
  (pss.gp.prbs.UP + pss.gp.prbs.POP) / 2
}

# Function to find possible genopair probabilities over multiple loci given
# multiple kinships
FindPssGPPsKPs = function(ale.frqs, L, knshp.st) {
  # Find possible genotype probabilities
  pss.gt.prbs = FindPssGtPrbs(ale.frqs)
  
  # Find possible first genotype probabilities for genopairs
  pss.gt.1.prbs = FindPssFrstGtPrbs(pss.gt.prbs, L)

  # Find possible genopair probabilities for unrelated pairs
  pss.gp.prbs.UP = FindPssGpPsUPs(pss.gt.1.prbs)
  
  # If parent-offspring pairs in kinship set selected
  if("Parent-offspring" %in% knshp.st) {
    # Find possible genopair probabilities for parent-offspring pairs
    pss.gp.prbs.POP = FindPssGpPsPOPs(ale.frqs, L, pss.gt.1.prbs)
  } else {
    pss.gp.prbs.POP = NULL
  }
  
  # If self-pairs in kinship set selected
  if("Self" %in% knshp.st) {
    # Function to find possible genopair probabilities for self-pairs
    pss.gp.prbs.SP = FindPssGpPsSPs(pss.gt.1.prbs, L)
  } else {
    pss.gp.prbs.SP = NULL
  }
  
  # If half-sibling pairs in kinship set selected
  if("Half-sibling" %in% knshp.st) {
    # Find possible genopair probabilities for half-sibling pairs
    pss.gp.prbs.HSP = FindPssGpPsHSPs(pss.gp.prbs.UP, pss.gp.prbs.POP)
  } else {
    pss.gp.prbs.HSP = NULL
  }
  
  # Possible genotypes for row and column-names
  pss.gts = c("00", "01", "11")
  
  # Combine genopair probabilities for selected kinships and return
  array(
    c(pss.gp.prbs.UP, pss.gp.prbs.HSP, pss.gp.prbs.POP, pss.gp.prbs.SP),
    dim = c(n.pss.gts, n.pss.gts, L, length(knshp.st) + 1),
    dimnames = list(
      gt.1 = pss.gts, gt.2 = pss.gts, Locus = paste0("L", 1:L),
      Kinship = c("UP", c("HSP", "POP", "SP")[rev(knshp.chcs) %in% knshp.st])
    )
  )
}