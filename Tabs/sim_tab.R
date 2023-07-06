# Simulation and outputs bound to simulate button

# Plot expected population size over time
output$nextExpPop = renderPlot({
  # Population size
  plot(
    sim.yrs.rct(), exp.N.t.rct(), ylim = c(0, max(exp.N.t.rct())),
    col = 'red', lwd = 2, type = 'l',
    xlab = 'Year', ylab = 'Population size', 
    main = "Expected population size over time"
  )
  
  # Base year and population
  abline(v = input$base.yr, h = input$exp.N.base, col = 2)
  
  # Surveys
  abline(v = srvy.yrs.rct(), lty = 2)
  
  # Add legend
  legend(
    "topleft", legend = c("Over time", "In base year", "Survey years"),
    col = c(2, 2, 1), lwd = c(2, 1, 1), lty = c(1, 1, 2)
  )
})

# Display implied parameter values
output$nextParsImpld = renderTable({
  frmt.pars.impld(
    lambda.rct(), beta.rct(), N.init.rct(), exp.N.t.rct()[input$hist.len], 
    exp.Ns.rct()
  )
}, digits = 3)

# Function to format table with one data type
FrmtTbl = function(data, rw.nms, cl.nms, md = "integer") {
  mode(data) = md
  df = data.frame(matrix(data, ncol = length(cl.nms)), row.names = rw.nms)
  names(df) = cl.nms
  df
}

# Predicted numbers of kin-pairs in population, reactive version. Transformed as
# matrix combining numbers within surveys and between survey-pairs
pred.ns.kps.pop.rct = reactive({
  # Get numbers
  preds = FindPredNsKPsPop(
    exp.N.t.rct(), s.yr.inds.rct(), input$phi, input$rho, lambda.rct(), 
    input$alpha, srvy.yrs.rct(), k.rct()
  )
  
  # Combine numbers within and between surveys in one table
  preds_wtn = t(preds$wtn)
  cbind(
    # Within surveys
    rbind(
      # Population sizes and total numbers of pairs
      preds_wtn[1:2, ], 
      
      # Self-pairs don't apply
      rep(NA, k.rct()), 
      
      # Other close-kin pairs
      preds_wtn[-(1:2), ]
    ),
    
    # Between survey-pairs
    rbind(
      # Population sizes don't apply
      rep(NA, n.srvy.prs.rct()), 
      
      # Other close-kin pairs
      t(preds$btn)[-5, ]
    )
  )
})

# Expected numbers of animals sampled
exp.ns.smp = reactive(pred.ns.kps.pop.rct()[1, 1:k.rct()] * input$p)

# Function to find and combine predicted numbers of offset or regular close-kin
# pairs among sampled animals, given the predicted total numbers of pairs
preds = reactive({
  function(pred.ns.APs) {
    rbind(
      # Numbers sampled
      c(exp.ns.smp(), rep(NA, n.srvy.prs.rct())), 
      
      # Total numbers of pairs
      pred.ns.APs, 
      
      # Close-kin pairs. Probabilities from predicted numbers in population,
      # multiplied by predicted total numbers of pairs.
      t(t(pred.ns.kps.pop.rct()[-(1:2), ]) / pred.ns.kps.pop.rct()[2, ] * 
          pred.ns.APs)
    )
  }
})

# Survey-years and survey-pairs
srvy.yrs.prs = reactive(c(srvy.yrs.rct(), srvy.prs.rct()))

# Predicted numbers of kin-pairs among sampled animals
output$predNsKPsSmpRct = renderTable({
  # Expected total numbers of pairs among sampled animals, within surveys, and
  # between survey-pairs
  pred.ns.APs.smp = c(
    choose(exp.ns.smp(), 2), combn(exp.ns.smp(), 2, function(ns) ns[1] * ns[2])
  )

  # Format as table of integers and output
  FrmtTbl(preds()(pred.ns.APs.smp), kp.tps[-4], srvy.yrs.prs())
}, rownames = T)

# Predicted numbers of kin-pairs among offset pairs of sampled animals
output$predNsKPsOffRct = renderTable({
  # Expected total numbers of offset pairs among sampled animals, within
  # surveys and between survey-pairs
  pred.ns.APs.offst = c(
    exp.ns.smp(), combn(exp.ns.smp(), 2, function(ns) max(ns[1], ns[2]))
  )

  # Format as table of numeric and output
  FrmtTbl(
    preds()(pred.ns.APs.offst), kp.tps[-4], srvy.yrs.prs(), md = "numeric"
  )
}, rownames = T)

# Output predicted numbers in population
output$predNsKPsPopRct = renderTable({
  FrmtTbl(pred.ns.kps.pop.rct(), kp.tps[-4], srvy.yrs.prs())
}, rownames = T)

# Simulate population and capture histories
bindEvent(observe({
  # Create list for population and capture histories
  hists.lst = vector("list", input$n.sims.rqd)
  
  # Create vectors for final and super-population sizes, and numbers of sample
  # histories
  N.fin.vec = Ns.vec = n.smp.hsts = n.smp.hsts = numeric(input$n.sims.rqd)
  any.empty = no.pairs = logical(input$n.sims.rqd)
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:input$n.sims.rqd) {
      # Simulate family and capture histories of population of animals over
      # time
      hists.lst[[hist.ind]] = stdy = SimPopStud(
        phi(), lambda(), N.init(), hist.len(), srvy.yrs(), k(), fnl.year(), p(),
        L(), imaf(), clvng.p(), tmp.emgn(), alpha(), clvng.ints()
      )
      
      # Collect final and super-population sizes, and numbers of sample
      # histories
      N.fin.vec[hist.ind] = tail(attributes(stdy)$N.t.vec, 1)
      Ns.vec[hist.ind] = attributes(stdy)$Ns
      no.pairs[hist.ind] = sum(attributes(stdy)$ns.caps) < 2
      
      # Update progress. Unexplained "Error in as.vector: object 'x' not
      # found" seen 19/12/2021 coming from incProgress...
      incProgress(1 / input$n.sims.rqd)
    }
  }, value = 0, message = "Simulating populations")
  
  # Check for populations that went extinct
  extinct = N.fin.vec == 0

  # Find studies that were OK
  stdy.ok = !extinct & !no.pairs

  # Update number of successful simulations
  n.sims(sum(stdy.ok))
  
  # Combine and return for studies that were OK
  sim.lst(list(
    hists.lst = hists.lst[stdy.ok], N.fin.vec = N.fin.vec[stdy.ok], 
    Ns.vec = Ns.vec[stdy.ok], n.extinct = sum(extinct), 
    n.no.pairs = sum(no.pairs)
  ))
}), input$simulate)

