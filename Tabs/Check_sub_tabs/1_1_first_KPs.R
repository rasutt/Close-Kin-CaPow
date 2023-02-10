# Outputs for first study kin-pairs sub-tab of checks tab

# Get attributes of first study
FS.atts = reactive(attributes(frst.std()))

# Population sizes
FS.Nts.wtn = reactive(FS.atts()$N.t.vec[s.yr.inds()])

# Function to format table of integers
FrmtTbl = function(data, rw.nms, cl.nms) {
  mode(data) = "integer"
  df = data.frame(matrix(data, ncol = length(cl.nms)), row.names = rw.nms)
  names(df) = cl.nms
  df
}

## Numbers of kin-pairs in whole population

# Simulated
output$firstNsKPsPop = renderTable({
  SPs.prnts.kwn = find.SPs.prnts.kwn(FS.atts(), k())
  SMPs = c(
    find.SMPs.wtn(FS.atts(), k()), 
    find.SMSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  )
  SFPs = c(
    find.SFPs.wtn(FS.atts(), k()),
    find.SFSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  )
  FSPs = c(
    find.FSPs.wtn(FS.atts(), k()),
    find.FSSPs.btn(FS.atts(), k()) - SPs.prnts.kwn
  )
  
  FrmtTbl(
    rbind(
      # Population sizes
      c(FS.Nts.wtn(), rep(NA, n.srvy.prs())),
      
      # Total numbers of pairs
      c(
        choose(FS.Nts.wtn(), 2),
        combn(FS.Nts.wtn(), 2, function(N.s.pr) N.s.pr[1] * N.s.pr[2])
      ), 
      
      # Self-pairs with parents unknown and known
      c(rep(NA, k()), find.SPs(FS.atts(), k())), 
      c(rep(NA, k()), SPs.prnts.kwn), 
      
      # Other kin-pairs
      c(find.POPs.wtn(FS.atts(), k()), find.POPs.btn(FS.atts(), k())),
      SMPs, SFPs, FSPs, SMPs + SFPs - 2 * FSPs
    ), 
    kp.tps, c(srvy.yrs(), srvy.prs())
  )
}, rownames = T)

# Predicted - repeats predictions for self-pairs as haven't implemented for
# parents unknown
frst.kp.preds = reactive({
  preds_wtn = t(pred.ns.kps.pop()$wtn)
  cbind(
    rbind(
      preds_wtn[1:2, ], matrix(NA, nrow = 2, ncol = k()), preds_wtn[-(1:2), ]
    ),
    rbind(
      rep(NA, n.srvy.prs()), t(pred.ns.kps.pop()$btn)[c(1:2, 2, 3:4, 6:8), ]
    )
  )
})

output$firstEstNsKPsPop = renderTable({
  FrmtTbl(frst.kp.preds(), kp.tps, c(srvy.yrs(), srvy.prs()))
}, rownames = T)

## Numbers of kin-pairs among sampled individuals

# Find numbers of kin pairs
frst.ns.kps.lst = reactive(FindNsKinPairs(k(), n.srvy.prs(), frst.std()))

# Simulated
frst.kps = reactive({
  wtn = rbind(
    FS.atts()$ns.caps, choose(FS.atts()$ns.caps, 2), NA,
    frst.ns.kps.lst()$wtn
  )
  
  APs.btn = combn(
    FS.atts()$ns.caps, 2, function(ns.caps.pr) ns.caps.pr[1] * ns.caps.pr[2]
  )
  btn = rbind(NA, APs.btn, frst.ns.kps.lst()$btn)
  
  cbind(wtn, btn)
})
output$firstNsKPsSmp = renderTable({
  FrmtTbl(
    frst.kps(), 
    c(
      "Total number sampled", "All pairs", "Self-pairs (all)", 
      "Parent-offspring pairs", "Half-sibling pairs"
    ), 
    c(srvy.yrs(), srvy.prs())
  )
}, rownames = T)

# Predicted (based on numbers sampled)
output$firstEstNsKPsSmp = renderTable({
  # Predictions for probabilities for population multiplied by total numbers of
  # pairs sampled
  preds = frst.kp.preds()[c(3, 5, 9), ] / frst.kp.preds()[2, ] * frst.kps()[2, ]

  FrmtTbl(
    preds, 
    c("Self-pairs (all)", "Parent-offspring pairs", "Half-sibling pairs"), 
    c(srvy.yrs(), srvy.prs())
  )
}, rownames = T)

