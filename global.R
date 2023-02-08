# Set bounds on inputs and step size for input slider and NLL plots
min_rho = 0
max_rho = 0.16
step_rho = 0.005
min_phi = 0.9
max_phi = 1
step_phi = 0.005

# Number of points at which to plot likelihood surfaces
n.pts = 100

# Model choices
mdl.chcs = c(
  "Popan", "Full true kinship", "Offset true kinship", "Full genopair", 
  "Offset genopair"
)

# Kinship choices for models
knshp.chcs = c("Self", "Parent-offspring", "Half-sibling")
all.knshps.bln = rep(1, 3)

# Genopair probability kinships
gpkts = c("Unrelated", rev(knshp.chcs))

# Set whether numbers of births stochastic
stch.bths <- T 

# Probability of males permanently emigrating each year. Note, code is not yet
# consistent with non-zero value, see readme.
pmt.emgn <- 0 

# Number of years to try to check intermediate estimators in derivations.
# Actual number limited by length of simulations
n.yrs.try.chk.t = 10

# Types of kin-pairs to be analysed (not same as those in models)
kp.tps = c(
  "Population sizes", "All pairs", "Self-pairs (all)", 
  "Self-pairs (parents known)", "Parent-offspring pairs", 
  "Same-mother pairs", "Same-father pairs", "Full-sibling pairs", 
  "Half-sibling pairs"
)
kp.tps.pop.wtn = kp.tps[c(1:2, 5:9)]
kp.tps.pop.btn = kp.tps[2:9]
kp.tps.prb.wtn = kp.tps[5:9]
kp.tps.prb.btn = kp.tps[c(3, 5:6)]
kp.tps.cap.wtn = kp.tps[c(5:6, 9)]
kp.tps.cap.btn = kp.tps[c(3, 5)]

kp.tps.t = c(
  "SMP{t=fnl,b1,b2=fnl}", "SFP{t=fnl,b1,b2=fnl}", "SFP{t=fnl,b}"
  # "SMP{t,f.yr,tm2,tm1}", "SMP{t,f.yr,tm1,f.yrm1}",
  # "SMP{t,f.yr,tm1,btwn.t.f.yr}", "SMP{t,f.yr,fst.yr,btwn}"
)

# Numbers of types of kin-pairs
n.kp.tps.pop.wtn = length(kp.tps.pop.wtn)
n.kp.tps.pop.btn = length(kp.tps.pop.btn)
n.kp.tps.cap.wtn = length(kp.tps.cap.wtn) 
n.kp.tps.cap.btn = length(kp.tps.cap.btn)
n.kp.tps.prb.wtn = length(kp.tps.prb.wtn)
n.kp.tps.prb.btn = length(kp.tps.prb.btn)
n.kp.tps.t = length(kp.tps.t)

# Descriptions and headings for kin-pairs within and between survey-years
wtn_btn_descs = c(
  "Pairs in which both individuals are alive in the same survey-year.", 
  "Pairs in which one individual is alive in each of two different survey-years."
)
wtn_btn_headings = c("Within survey-years", "Between two survey-years")

# Kin-pair names and descriptions
kp.nms = c(
  "Population sizes", "All pairs", "Self-pairs", "Parent-offspring pairs", 
  "Same-mother pairs", "Same-father pairs", "Full-sibling pairs", 
  "Half-sibling pairs"
)
names(kp.nms) = c("N", "APs", "SPs", "POPs", "SMPs", "SFPs", "FSPs", "HSPs")
rglr.kp.dscs = c(
  APs = "Total numbers of pairs of individuals.",
  POPs = "Numbers of pairs of individuals that are parent and offspring.",
  SMPs = "Numbers of pairs of individuals with the same mother.",
  SFPs = "Numbers of pairs of individuals with the same father.",
  FSPs = "Numbers of pairs of individuals with the same parents.",
  HSPs = "Numbers of pairs of individuals that share exactly one parent."
)

# IDs for regular kin-pairs, those existing both within and between survey-years
rglr.kp.ids = c("APs", "POPs", "SMPs", "SFPs", "FSPs", "HSPs")

# Set possible genotypes at each locus, representing 0 and 1-coded SNP
# alleles as 1 and 2 respectively, to index corresponding allele frequencies
pss.gt.lbls = c("00", "01", "11")
pss.gts = cbind(c(1, 1), 1:2, c(2, 2))
ales.1.inds = pss.gts[1, ]
ales.2.inds = pss.gts[2, ]
n.pss.gts = 3

# Find possible genopairs by indexing possible genotypes to get all combinations
pss.gts.1 = pss.gts[, rep(1:n.pss.gts, n.pss.gts)]
pss.gts.2 = pss.gts[, rep(1:n.pss.gts, each = n.pss.gts)]


