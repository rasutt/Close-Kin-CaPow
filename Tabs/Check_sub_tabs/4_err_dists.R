# Outputs for error distributions sub-tab of checks tab

# Kin-pair probability estimates
exp.probs.KPs.lst = reactive({
  list(
    SPs.btn = est.ns.kps.pop.lst()$btn[2, ] / est.ns.kps.pop.lst()$btn[1, ],
    SMPs.wtn = est.ns.kps.pop.lst()$btn[3, ] / est.ns.kps.pop.lst()$btn[1, ]
  )
})

# Kin-pair probabilities observed
probs.KPs.lst = reactive({
  list(
    SPs.btn = checks.lst()$ns.kps.pop.btn.arr[, , 2] / 
      checks.lst()$ns.kps.pop.btn.arr[, , 1],
    SMPs.wtn = checks.lst()$ns.kps.pop.btn.arr[, , 3] / 
      checks.lst()$ns.kps.pop.btn.arr[, , 1]
  )
})

# Indices for types of kin-pairs
KP_pop_inds = c(
  "Survey", "Survey", "Survey-pair", "Survey-pair", "Survey", "Survey", 
  "Survey", "Survey", "Survey-pair"
)
KP_prob_inds = c("Survey-pair", "Survey")
KP_inds = c("Survey-pair", "Survey", "Survey-pair", "Survey", "Survey")

# Function to plot simulated versus expected numbers of kin-pairs for one type
# of kin-pair
nsKPsPlot = function(i, pop = F, prob = F) {
  if (pop) {
    diffs = ns.kps.pop.wtn.est.errs()[, , i]
    xlab = KP_pop_inds[i]
    main = kp.tps.pop.wtn[i]
  } else {
    diffs = checks.lst()$ns.KPs.lst[[i]] / checks.lst()$exp.ns.KPs.lst[[i]] - 1 
    xlab = KP_inds[i]
    main = KP_names[i]
  }
  if (prob) {
    diffs = 
      t(t(probs.KPs.lst()[[i]]) / exp.probs.KPs.lst()[[i]] - 1)
    xlab = KP_prob_inds[i]
    main = KP_prob_names[i]
  }
  if (xlab == "Survey") colnames(diffs) = srvy.yrs()
  else colnames(diffs) = apply(combn(srvy.yrs(), 2), 2, paste, collapse = "-")
  {
    boxplot(
      diffs, main = main, xlab = xlab, 
      ylab = "(Observed - expected numbers) / expected numbers"
    )
    abline(h = 0, col = 'red')
    abline(h = mean(diffs), col = 'blue')
    legend(
      "topleft", col = c(2, 4), lty = 1,
      legend = c(
        "Expected difference over all surveys (zero)", 
        "Average difference over all surveys"
      ),
    )
  }
}

# Apply plot function to each type of kin-pair

# All pairs within survey years for whole population
output$nsAPsWtnPop = renderPlot(nsKPsPlot(1, T))
# # All pairs between survey years for whole population
# output$nsAPsBtnPop = renderPlot(nsKPsPlot(2, T))
# # Self-pairs between survey years for whole population
# output$nsSPsBtnPop = renderPlot(nsKPsPlot(3, T))
# # Same-mother pairs within survey years for whole population
# output$nsSMPsWtnPop = renderPlot(nsKPsPlot(4, T))
# # Same-father pairs within survey years for whole population
# output$nsSFPsWtnPop = renderPlot(nsKPsPlot(5, T))
# # Same-mother pairs between survey years for whole population
# output$nsSMPsBtnPop = renderPlot(nsKPsPlot(9, T))

# # Self-pair probabilities between survey years for whole population
# output$probSPsBtnPop = renderPlot(nsKPsPlot(1, prob = T))
# # Same-mother pair probabilities within survey years for whole population
# output$probSMPsWtnPop = renderPlot(nsKPsPlot(2, prob = T))

# # Self-pairs between samples
# output$nsSPsBtn = renderPlot(nsKPsPlot(1))
# # Parent-offspring pairs within samples
# output$nsPOPsWtn = renderPlot(nsKPsPlot(2))
# # Parent-offspring pairs between samples
# output$nsPOPsBtn = renderPlot(nsKPsPlot(3))
# # Same-mother pairs within samples
# output$nsSMPsWtn = renderPlot(nsKPsPlot(4))
# # Half-sibling pairs within samples
# output$nsHSPsWtn = renderPlot(nsKPsPlot(5))
