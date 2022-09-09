# Outputs for error distributions sub-tab of checks tab

# Function to plot simulated versus expected numbers of kin-pairs for one type
# of kin-pair
nsKPsPlot = function(errs, kp.type) {
  boxplot(
    errs, main = kp.type, xlab = names(dimnames(errs))[2],
    ylab = "Proportional errors"
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

# Same-mother pairs in the final year, with one born in the year indicated, and
# one born in the final year
output$nsSMPsFnlB1Fnl = renderPlot(
  nsKPsPlot(ns.kps.t.est.errs()[, , 1], kp.tps.t[1])
)

# Same-father pairs in the final year, with one born in the year indicated, and
# one born in the final year
output$nsSFPsFnlB1Fnl = renderPlot(
  nsKPsPlot(ns.kps.t.est.errs()[, , 2], kp.tps.t[2])
)

# Same-father pairs in the final year, with one born in the year indicated, and
# one born in the final year
output$nsSFPsFnlB = renderPlot(
  nsKPsPlot(ns.kps.t.est.errs()[, , 3], kp.tps.t[3])
)

# Within surveys
output$nsAPsWtnPop = renderPlot(
  nsKPsPlot(
    ns.kps.pop.wtn.est.errs()[, , 2],
    kp.tps.pop.wtn[2]
  )
)

# Between surveys
output$nsAPsBtnPop = renderPlot(
  nsKPsPlot(
    ns.kps.pop.btn.est.errs()[, , 1],
    kp.tps.pop.btn[1]
  )
)

# # Self-pairs between survey years 
# output$nsSPsBtnPop = renderPlot(nsKPsPlot(3, T))
# # Same-mother pairs within survey years 
# output$nsSMPsWtnPop = renderPlot(nsKPsPlot(4, T))
# # Same-father pairs within survey years 
# output$nsSFPsWtnPop = renderPlot(nsKPsPlot(5, T))
# Same-mother pairs between survey years
output$nsSMPsBtnPop = renderPlot(
  nsKPsPlot(
    ns.kps.pop.btn.est.errs()[, , 3],
    kp.tps.pop.btn[3]
  )
)

## Probabilities

# # Self-pair probabilities between survey years 
# output$probSPsBtnPop = renderPlot(nsKPsPlot(1, prob = T))
# Same-mother pair probabilities within survey years 
output$probSMPsWtnPop = renderPlot(
  nsKPsPlot(
    ns.kps.prb.wtn.est.errs()[, , 1],
    kp.tps.prb.wtn[1]
  )
)

## Numbers in samples

# # Self-pairs between surveys
# output$nsSPsBtn = renderPlot(nsKPsPlot(1))
# Parent-offspring pairs within surveys
# output$nsPOPsWtn = renderPlot(nsKPsPlot(2))
# # Parent-offspring pairs between surveys
# output$nsPOPsBtn = renderPlot(nsKPsPlot(3))
# # Same-mother pairs within surveys
# output$nsSMPsWtn = renderPlot(nsKPsPlot(4))
# # Half-sibling pairs within surveys
# output$nsHSPsWtn = renderPlot(nsKPsPlot(5))
