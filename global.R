# Set bounds on inputs and step size for input slider and NLL plots
min_rho = 0
max_rho = 0.16
step_rho = 0.005
min_phi = 0.9
max_phi = 1
step_phi = 0.005

# Model choices
mod.choices = c("Popan", "Close-kin")

# Set simulation parameters for basic scenario
stch.bths <- T # Set whether numbers of births stochastic

# Probability of males permanently emigrating each year. Note, code is not yet
# consistent with non-zero value, see readme.
pmt.emgn <- 0 

# Number of years to try to check intermediate estimators in derivations.
# Actual number limited by length of simulations
n.yrs.try.chk.t = 20

# Types of kin-pairs to be analysed (not same as those in models)
kp.tps = c(
  "Population sizes", "All-pairs", "Self-pairs", "Parent-offspring pairs", 
  "Same-mother pairs", "Same-father pairs", "Full-sibling pairs", 
  "Half-sibling pairs"
)
kp.tps.pop.wtn = kp.tps[c(1:2, 4:8)]
kp.tps.pop.btn = kp.tps[c(2:3, 5)]
kp.tps.prb.wtn = kp.tps[4:8]
kp.tps.prb.btn = kp.tps[c(3, 5)]
kp.tps.cap.wtn = kp.tps[c(4:5, 8)]
kp.tps.cap.btn = kp.tps[3:4]

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
