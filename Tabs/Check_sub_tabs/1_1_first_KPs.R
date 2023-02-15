# Outputs for first study kin-pairs sub-tab of checks tab

# Get attributes of first study
FS.atts = reactive(attributes(frst.std()))

# Population sizes
FS.Nts.wtn = reactive(FS.atts()$N.t.vec[s.yr.inds()])

# Numbers of kin-pairs in population
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

# Numbers of kin-pairs among sampled individuals
frst.ns.kps.lst = reactive(FindNsKinPairs(k(), n.srvy.prs(), frst.std()))
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

