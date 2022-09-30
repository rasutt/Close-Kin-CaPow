# Outputs for error distributions sub-tab of checks tab

# Function to plot simulated versus expected numbers of kin-pairs for one type
# of kin-pair
nsKPsPlot = function(errs, kp.type) {
  boxplot(
    errs, main = kp.type, xlab = names(dimnames(errs))[2],
    ylab = "Proportional errors", show.names = T
  )
  abline(h = 0, col = 'red')
  abline(h = mean(errs), col = 'blue')
  legend(
    "topleft", col = c(2, 4), lty = 1,
    legend = c("Estimated error (zero)", "Average error"),
  )
}

### Plot kin-pair estimate error distributions

## Numbers in whole population

# Temporal estimates

# Same-father pairs in the final year, with one born in the year indicated, and
# one born in the final year
output$nsSFPsFnlB1Fnl = renderPlot(
  nsKPsPlot(ns.kps.t.errs()[, , 2], kp.tps.t[2])
)

# Same-father pairs in the final year, with one born in the year indicated, and
# one born in the final year
output$nsSFPsFnlB = renderPlot(
  nsKPsPlot(ns.kps.t.errs()[, , 3], kp.tps.t[3])
)

## Probabilities

# # Self-pair probabilities between survey years 
# output$probSPsBtnPop = renderPlot(nsKPsPlot(1, prob = T))
# Same-mother pair probabilities within survey years 
output$probSMPsWtnPop = renderPlot(
  nsKPsPlot(ns.kps.prb.wtn.errs()[, , 1], kp.tps.prb.wtn[1])
)
