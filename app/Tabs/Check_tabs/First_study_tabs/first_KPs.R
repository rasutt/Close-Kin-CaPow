# Outputs for first study kin-pairs sub-tab of checks tab

# Get attributes of first study
FS.atts = reactive(attributes(frst.std()))

# Population sizes
FS.Nts.wtn = reactive(FS.atts()$N.t.vec[s.yr.inds()])

# Numbers of kin-pairs among all sampled individuals, and offset pairs
frst.ns.skps.lst = reactive(FindNsKinPairs(k(), n.srvy.prs(), frst.std()))
frst.ns.okps.lst = reactive({
  FindNsOKPs(k(), n.srvy.prs(), frst.std(), frst.osiips(), frst.osyips())
})

# Function to combine kinpair numbers within and between surveys
combnKPs = function(ns.kps.lst, offset = F) {
  wtn = rbind(
    FS.atts()$ns.caps, 
    if (offset) {
      FS.atts()$ns.caps
    } else {
      choose(FS.atts()$ns.caps, 2)
    },
    NA,
    ns.kps.lst$wtn
  )
  
  btn = rbind(
    NA, 
    if (offset) {
      combn(
        FS.atts()$ns.caps, 2, 
        function(ns.caps.pr) max(ns.caps.pr[1], ns.caps.pr[2])
      )
    } else {
      combn(
        FS.atts()$ns.caps, 2, function(ns.caps.pr) ns.caps.pr[1] * ns.caps.pr[2]
      )
    },
    ns.kps.lst$btn
  )
  
  FrmtTbl(
    cbind(wtn, btn), 
    c("Total number sampled", "All pairs", "Self-pairs (all)", 
      "Parent-offspring pairs", "Half-sibling pairs"), 
    c(srvy.yrs(), srvy.prs())
  )
}

# Show kinpair numbers among all sampled individuals, and offset pairs
output$firstNsKPsSmp = renderTable({
  combnKPs(frst.ns.skps.lst())
}, rownames = T)
output$firstNsKPsOfst = renderTable({
  combnKPs(frst.ns.okps.lst(), offset = T)
}, rownames = T)

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

