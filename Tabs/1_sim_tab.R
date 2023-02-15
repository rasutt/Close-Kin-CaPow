# Simulation and outputs bound to simulate button

# Plot expected population size over time
output$nextExpPop <- renderPlot({
  plot(
    sim.yrs.rct(), exp.N.t.rct(), 
    col = 'red', lwd = 2, type = 'l',
    xlab = 'Year', ylab = 'Population size', 
    main = "Expected population size over time"
  )
  # Base year
  abline(v = input$base.yr, h = input$exp.N.base, col = 2)
  # Surveys
  abline(v = srvy.yrs.rct(), lty = 2)
  # Add legend
  legend(
    "topleft", legend = c("Over time", "In base year", "Survey years"),
    col = c(2, 2, 1), lwd = c(2, 1, 1), lty = c(1, 1, 2)
  )
})

# if (beta < 0 | beta > 1)
#   cat("Implied birthrate for mature females:", round(beta, 3), "\n")
# if (beta < 0) stop("Negative birth rates impossible")
# if (beta > 1) stop("Maximum one calf at a time")

# Display implied parameter values
output$nextParsImpld <- renderTable({
  frmt.pars.impld(
    lambda.rct(), beta.rct(), exp.N.t.rct()[input$hist.len], exp.Ns.rct()
  )
}, digits = 3)

# Function to format table of integers
FrmtTbl = function(data, rw.nms, cl.nms) {
  mode(data) = "integer"
  df = data.frame(matrix(data, ncol = length(cl.nms)), row.names = rw.nms)
  names(df) = cl.nms
  df
}

# Predicted numbers of kin-pairs for whole population
pred.ns.kps.pop.rct = reactive({
  FindPredNsKPsPop(
    exp.N.t.rct(), s.yr.inds.rct(), input$phi, input$rho, lambda.rct(), 
    input$alpha, srvy.yrs.rct(), k.rct()
  )
})

# Predicted numbers of kin-pairs in population - repeats predictions for
# self-pairs as haven't implemented for parents unknown
frst.kp.preds = reactive({
  preds_wtn = t(pred.ns.kps.pop.rct()$wtn)
  cbind(
    rbind(
      preds_wtn[1:2, ], matrix(NA, nrow = 2, ncol = k.rct()), 
      preds_wtn[-(1:2), ]
    ),
    rbind(
      rep(NA, n.srvy.prs.rct()), 
      t(pred.ns.kps.pop.rct()$btn)[c(1:2, 2, 3:4, 6:8), ]
    )
  )
})

# Predicted numbers of kin-pairs among sampled individuals
output$firstEstNsKPsSmp = renderTable({
  # Predictions for probabilities for population multiplied by expected total
  # numbers of pairs sampled
  exp.n.smps = pred.ns.kps.pop.rct()$wtn[, 1] * input$p
  exp.n.APs.smps = c(
    choose(exp.n.smps, 2), combn(exp.n.smps, 2, function(ns) ns[1] * ns[2])
  )
  preds = rbind(
    c(exp.n.smps, rep(NA, n.srvy.prs.rct())), exp.n.APs.smps, 
    frst.kp.preds()[c(3, 5, 9), ] / frst.kp.preds()[2, ] * exp.n.APs.smps
  )
  
  FrmtTbl(
    preds, 
    c("Number of samples", "Total number of pairs", "Self-pairs (all)", 
      "Parent-offspring pairs", "Half-sibling pairs"), 
    c(srvy.yrs.rct(), srvy.prs.rct())
  )
}, rownames = T)

# Predicted numbers of kin-pairs among offset pairs
output$firstEstNsKPsOff = renderTable({
  # Predictions for probabilities for population multiplied by expected total
  # numbers of pairs sampled
  exp.n.smps = pred.ns.kps.pop.rct()$wtn[, 1] * input$p
  exp.n.APs.offst = c(
    exp.n.smps, combn(exp.n.smps, 2, function(ns) max(ns[1], ns[2]))
  )
  preds = rbind(
    c(exp.n.smps, rep(NA, n.srvy.prs.rct())), exp.n.APs.offst, 
    frst.kp.preds()[c(3, 5, 9), ] / frst.kp.preds()[2, ] * exp.n.APs.offst
  )
  
  FrmtTbl(
    preds, 
    c("Number of samples", "Total number of pairs", "Self-pairs (all)", 
      "Parent-offspring pairs", "Half-sibling pairs"), 
    c(srvy.yrs.rct(), srvy.prs.rct())
  )
}, rownames = T)

# Output predicted numbers in population
output$firstEstNsKPsPop = renderTable({
  FrmtTbl(frst.kp.preds(), kp.tps, c(srvy.yrs.rct(), srvy.prs.rct()))
}, rownames = T)

# Simulate population and capture histories
bindEvent(observe({
  # Initial population size
  N.init = round(exp.N.t()[1])
  
  # Create list for population and capture histories
  hists.lst <- vector("list", n.sims())
  # Create vectors for final and super-population sizes
  N.fin.vec <- Ns.vec <- numeric(n.sims())
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n.sims()) {
      # Simulate family and capture histories of population of animals over
      # time
      hists.lst[[hist.ind]] <- SimPopStud(
        phi(), lambda(), N.init, hist.len(), srvy.yrs(), k(), fnl.year(), p(),
        L(), clvng.p(), tmp.emgn(), alpha(), clvng.ints()
      )
      # Collect final and super-population sizes
      N.fin.vec[hist.ind] <- tail(attributes(hists.lst[[hist.ind]])$N.t.vec, 1)
      Ns.vec[hist.ind] <- attributes(hists.lst[[hist.ind]])$Ns
      
      # Update progress. Unexplained "Error in as.vector: object 'x' not
      # found" seen 19/12/2021 coming from incProgress...
      incProgress(1/n.sims())
    }
  }, value = 0, message = "Simulating populations")
  
  sim.lst(list(hists.lst = hists.lst, N.fin.vec = N.fin.vec, Ns.vec = Ns.vec))
}), input$simulate)

