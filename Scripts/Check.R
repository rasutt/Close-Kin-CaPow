# Check simulated studies

# Set number of studies to check
n.stds.chk <- 1000

# Load parameters and simulated studies
# load(file = paste0("Saved R data/", sim.name, "_sim.Rdata"))

# Load kin pair functions
source("Functions/FindNsKinPairs.R")
source("Functions/FindExpNsKPs.R")

# Create matrices for population trajectories (as columns), numbers captured,
# and expected and observed number of kin pairs
N.t.mat <- matrix(nrow = n.stds.chk, ncol = hist.len)
ns.caps.mat <- ns.clvng.caps.mat <- ns.clvng.mat <- ns.POPs.wtn.mat <- 
  ns.HSPs.wtn.mat <- exp.ns.HSPs.wtn.mat <- 
  exp.ns.POPs.wtn.mat <- matrix(nrow = n.stds.chk, ncol = k)
ns.POPs.btn.mat <- ns.SPs.btn.mat <- exp.ns.POPs.btn.mat <- 
  exp.ns.SPs.btn.mat <- matrix(nrow = n.stds.chk, ncol = n.srvy.prs)

# Create vectors for proportions with unknown parents, and superpopulation sizes
prpn.prnts.unkn.vec <- Ns.vec <- numeric(n.stds.chk)

# Display progress
cat("Checking study: ")

# Loop over histories
for (hist.ind in 1:n.stds.chk) {
  # Display progress
  if (hist.ind %% 100 == 1) cat(hist.ind, "")
  
  # Get simulated family and capture histories of population of animals over
  # time
  pop.cap.hist <- hists.lst[[hist.ind]]
  
  # Record population curve
  N.t.mat[hist.ind, ] <- attributes(pop.cap.hist)$N.t.vec
  
  # Record superpopulation size
  Ns.vec[hist.ind] <- attributes(pop.cap.hist)$Ns
  
  # Get numbers captured and calving in each survey
  ns.caps <- attributes(pop.cap.hist)$ns.caps
  ns.caps.mat[hist.ind, ] <- ns.caps
  ns.clvng.mat[hist.ind, ] <- attributes(pop.cap.hist)$ns.clvng
  ns.clvng.caps.mat[hist.ind, ] <- 
    colSums(pop.cap.hist[, 4:(3 + k)] * pop.cap.hist[, (4 + k):(3 + 2 * k)])
    
  # Find proportion captured with unknown parents
  prpn.prnts.unkn.vec[hist.ind] <- mean(is.na(pop.cap.hist$mum))
  
  # Find numbers of known kin pairs
  ns.kps.lst <- FindNsKinPairs()
  
  # Record in matrices
  ns.POPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.wtn
  ns.HSPs.wtn.mat[hist.ind, ] <- ns.kps.lst$ns.HSPs.wtn
  ns.POPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.POPs.btn
  ns.SPs.btn.mat[hist.ind, ] <- ns.kps.lst$ns.SPs.btn
  
  # Find expected numbers of kin pairs
  exp.ns.kps.lst <- FindExpNsKPs()

  # Record in matrices
  exp.ns.POPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.wtn
  exp.ns.HSPs.wtn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.HSPs.wtn
  exp.ns.POPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.POPs.btn
  exp.ns.SPs.btn.mat[hist.ind, ] <- exp.ns.kps.lst$exp.ns.SPs.btn
}

# Display progress
cat(hist.ind, "\n")

# Plot population trajectories and expected value
matplot(t(N.t.mat), type = 'l', col = rgb(0, 0, 0, alpha = 0.1), lty = 1, 
        xlab = 't', ylab = 'Nt', main = "Population size over time")
lines(exp.N.t, col = 'red', lwd = 2)

# Get vector of final population sizes
N.fin.vec <- N.t.mat[, hist.len]

# Plot final population sizes, mean, and expected value
boxplot(N.fin.vec, main = "Final population size")
abline(h = mean(N.fin.vec), col = 'blue')
abline(h = exp.N.t[hist.len], col = 'red')

# Plot superpopulation sizes
boxplot(Ns.vec, main = "Superpopulation size")
abline(h = mean(Ns.vec), col = 'blue')
abline(h = exp.Ns, col = "red")

# Check proportions captured match capture probabilities

# Non-calving animals
boxplot((ns.caps.mat - ns.clvng.caps.mat) / 
          (N.t.mat[, hist.len + srvy.yrs - f.year] - ns.clvng.mat),
        main = "Non-calving animals", xlab = "survey", 
        ylab = "proportion captured")
abline(h = p, col = 'red')

# Calving animals
boxplot(ns.clvng.caps.mat / ns.clvng.mat,
        main = "Calving animals", xlab = "survey", 
        ylab = "proportion captured")
abline(h = p + clvng.p, col = 'red')

# Proportions of animals calving in each survey year
boxplot(ns.clvng.mat / N.t.mat[, hist.len + srvy.yrs - f.year],
        main = "Proportions of animals calving", xlab = "survey", 
        ylab = "proportion")

# Display proportion of captures for which the parents are unknown
cat("Proportion of captures for which the parents are unknown:", 
    mean(prpn.prnts.unkn.vec), "\n")

# Plot differences between observed and expected numbers of kinpairs

# Parent-offspring pairs within samples
boxplot(ns.POPs.wtn.mat - exp.ns.POPs.wtn.mat,
        main = "POPs within samples", xlab = "survey", 
        ylab = "observed - expected")
abline(h = mean(ns.POPs.wtn.mat - exp.ns.POPs.wtn.mat), col = 'blue')
abline(h = 0, col = 'red')
cat("Mean difference between observed and expected numbers of 
    POPs within samples:",
    mean(ns.POPs.wtn.mat - exp.ns.POPs.wtn.mat), "\n")

# Half-sibling pairs within samples
boxplot(ns.HSPs.wtn.mat - exp.ns.HSPs.wtn.mat,
        main = "HSPs within samples", xlab = "survey", 
        ylab = "observed - expected")
abline(h = mean(ns.HSPs.wtn.mat - exp.ns.HSPs.wtn.mat), col = 'blue')
abline(h = 0, col = 'red')
cat("Mean difference between observed and expected numbers of 
    HSPs within samples:",
    mean(ns.HSPs.wtn.mat - exp.ns.HSPs.wtn.mat), "\n")

# Parent-offspring pairs between samples
boxplot(ns.POPs.btn.mat - exp.ns.POPs.btn.mat,
        main = "POPs between samples", xlab = "survey pair",
        ylab = "observed - expected")
abline(h = mean(ns.POPs.btn.mat - exp.ns.POPs.btn.mat), col = 'blue')
abline(h = 0, col = 'red')
cat("Mean difference between observed and expected numbers of 
    POPs between samples:",
    mean(ns.POPs.btn.mat - exp.ns.POPs.btn.mat), "\n")

# Self-pairs between samples
boxplot(ns.SPs.btn.mat - exp.ns.SPs.btn.mat,
        main = "SPs between samples", xlab = "survey pair",
        ylab = "observed - expected")
abline(h = mean(ns.SPs.btn.mat - exp.ns.SPs.btn.mat), col = 'blue')
abline(h = 0, col = 'red')
cat("Mean difference between observed and expected numbers of 
    SPs between samples:",
    mean(ns.SPs.btn.mat - exp.ns.SPs.btn.mat), "\n")

# The downward skew of observed numbers of kinpairs is associated with large
# population sizes
cor(N.fin.vec, (ns.POPs.wtn.mat - exp.ns.POPs.wtn.mat)[, 1])

boxplot((ns.POPs.wtn.mat - exp.ns.POPs.wtn.mat)[N.fin.vec < 3000, ],
        main = "POPs within samples", xlab = "survey",
        ylab = "observed - expected")

matplot(t(N.t.mat[(ns.POPs.wtn.mat - exp.ns.POPs.wtn.mat)[, 1] < -20, ]),
        type = 'l', col = rgb(1, 0, 0, alpha = 0.2), lty = 1,
        xlab = 't', ylab = 'Nt', 
        main = "Population size over time when many POPs within samples")
matlines(t(N.t.mat[(ns.POPs.wtn.mat - exp.ns.POPs.wtn.mat)[, 1] > -20, ]),
        col = rgb(0, 0, 1, alpha = 0.02))

# boxplot((ns.HSPs.wtn.mat - exp.ns.HSPs.wtn.mat)[N.fin.vec < 2200, ],
#         main = "HSPs within samples", xlab = "survey",
#         ylab = "observed - expected")
