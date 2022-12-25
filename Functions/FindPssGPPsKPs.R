# Function to find possible genopair probabilities over multiple loci given
# multiple kinships
FindPssGPPsKPs = function(ale.frqs, L) {
  ales.1.inds = pss.gts[1, ]
  ales.2.inds = pss.gts[2, ]
  
  # First genotypes, found by indexing the 2 x L allele frequencies matrix for
  # each allele of each possible genotype (globally defined for SNP genotypes),
  # and multiplying by 2 possible cases for heterozygous genotypes.  Filled into
  # a 3 x 1 x L array which is indexed 3 times to fill the three columns.
  pss.gt.prbs = matrix(
    ale.frqs[ales.1.inds, ] * ale.frqs[ales.2.inds, ] * 
      (1 + (ales.1.inds != ales.2.inds)), 
    nrow = n.pss.gts
  )
  pss.gt.2.prbs = array(
    rep(pss.gt.prbs, each = n.pss.gts), 
    c(n.pss.gts, n.pss.gts, L)
  )
  pss.gt.1.prbs = aperm(pss.gt.2.prbs, c(2, 1, 3))
  
  # Conditional probabilities given that the pair are unrelated, the products of
  # the respective genotype probabilities, the second found by permuting the
  # array containing the first.
  pss.gp.prbs.UP = pss.gt.1.prbs * pss.gt.2.prbs
  
  # Conditional probabilities of second genotype given that the pair are parent
  # and offspring (unordered).  0.5 for each allele in first genotype being
  # inherited, multiplied by the probability of the second genotype in each
  # case, as in table 3, pg. 269, Bravingtion et al. (2016) Close-Kin
  # Mark-Recapture. Filled into an L x 3 x 3 array which is permuted to the
  # standard dimensions.
  pss.cnd.gt.2.prbs.POP = aperm(
    array(
      # Note: order data enters array is down columns, not across rows
      c(
        ale.frqs[1, ], 0.5 * ale.frqs[1, ], rep(0, L), 
        ale.frqs[2, ], 0.5 * colSums(ale.frqs), ale.frqs[1, ],
        rep(0, L), 0.5 * ale.frqs[2, ], ale.frqs[2, ]
      ), c(L, n.pss.gts, n.pss.gts)
    ), c(2, 3, 1)
  )
  
  # Conditional probabilities given that the pair are parent and offspring
  # (unordered). Products of the first genotype probabilities and the
  # conditional probabilities of the second genotypes given that the pair are
  # parent and offspring.
  pss.gp.prbs.POP = pss.gt.1.prbs * pss.cnd.gt.2.prbs.POP
  
  # Conditional probabilities given that the pair are the same individual
  # sampled twice. Genotype probabilities when the genotypes are the same, and
  # zero otherwise.
  pss.gp.prbs.SP = aperm(
    array(
      cbind(
        pss.gt.1.prbs[1, 1, ], 0, 0, 0, pss.gt.1.prbs[2, 2, ], 0, 0, 0, 
        pss.gt.1.prbs[3, 3, ]
      ), c(L, n.pss.gts, n.pss.gts)
    ), c(2, 3, 1)
  )
  
  # Conditional probabilities given that the pair are half-siblings.  Average of
  # probabilities for unrelated and parent-offspring pairs.
  pss.gp.prbs.HSP = (pss.gp.prbs.UP + pss.gp.prbs.POP) / 2
  
  pss.gts = c("00", "01", "11")
  array(
    c(pss.gp.prbs.UP, pss.gp.prbs.HSP, pss.gp.prbs.POP, pss.gp.prbs.SP),
    dims = c(n.pss.gts, n.pss.gts, L, 4),
    dimnames = list(
      gt.1 = pss.gts, gt.2 = pss.gts, Locus = paste0("L", 1:L), 
      Kinship = c("UP", "HSP", "POP", "SP")
    )
  )
}