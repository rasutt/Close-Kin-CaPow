# Save/load tab

# Save objects
output$downloadData <- downloadHandler(
  filename = "ckc_saved_objs.Rdata",
  content = function(file) {
    saved.objs = list(
      phi = phi(),
      rho = rho(),
      lambda = lambda(),
      base.yr = base.yr(),
      exp.N.base = exp.N.base(),
      srvy.yrs = srvy.yrs(),
      p = p(),
      L = L(),
      clvng.ints = clvng.ints(),
      clvng.p = clvng.p(),
      tmp.emgn = tmp.emgn(),
      alpha = alpha(),
      hist.len = hist.len(),
      n.sims = n.sims(),
      srvy.gaps = srvy.gaps(),
      k = k(),
      n.srvy.prs = n.srvy.prs(),
      srvy.prs = srvy.prs(),
      fnl.year = fnl.year(),
      fst.year = fst.year(),
      sim.yrs = sim.yrs(),
      n.yrs.chk.t = n.yrs.chk.t(),
      yrs.chk.t = yrs.chk.t(),
      s.yr.inds = s.yr.inds(),
      exp.N.t = exp.N.t(),
      exp.N.fin = exp.N.fin(),
      exp.Ns = exp.Ns(),
      par.names = par.names(),
      par.vals = par.vals(),
      est.par.names = est.par.names(),
      sim.opts.bio.scen = sim.opts.bio.scen(),
      sim.lst = sim.lst(),
      checks.lst = checks.lst(),
      fit.lst = fit.lst(),
      N.t.mat = N.t.mat(),
      avg.phi.obs = avg.phi.obs(),
      ns.SPs = ns.SPs(),
      ns.POPs = ns.POPs(),
      ns.SMPs = ns.SMPs(),
      ns.SFPs = ns.SFPs(),
      ns.SMPs.t = ns.SMPs.t(),
      ns.SFPs.t = ns.SFPs.t(),
      ns.SibPs = ns.SibPs(),
      pns.UPs = pns.UPs(),
      lg.gp.prbs.KP = lg.gp.prbs.KP()
    )
    save(saved.objs, file = file)
  }
)

# Increase maximum file upload size to 50MB
options(shiny.maxRequestSize=100*1024^2)

# Load saved objects
observeEvent(input$file, {
  load(input$file$datapath)
  
  phi(saved.objs$phi)
  rho(saved.objs$rho)
  lambda(saved.objs$lambda)
  base.yr(saved.objs$base.yr)
  exp.N.base(saved.objs$exp.N.base)
  srvy.yrs(saved.objs$srvy.yrs)
  p(saved.objs$p)
  L(saved.objs$L)
  clvng.ints(saved.objs$clvng.ints)
  clvng.p(saved.objs$clvng.p)
  tmp.emgn(saved.objs$tmp.emgn)
  alpha(saved.objs$alpha)
  hist.len(saved.objs$hist.len)
  n.sims(saved.objs$n.sims)
  srvy.gaps(saved.objs$srvy.gaps)
  k(saved.objs$k)
  n.srvy.prs(saved.objs$n.srvy.prs)
  srvy.prs(saved.objs$srvy.prs)
  fnl.year(saved.objs$fnl.year)
  fst.year(saved.objs$fst.year)
  sim.yrs(saved.objs$sim.yrs)
  n.yrs.chk.t(saved.objs$n.yrs.chk.t)
  yrs.chk.t(saved.objs$yrs.chk.t)
  s.yr.inds(saved.objs$s.yr.inds)
  exp.N.t(saved.objs$exp.N.t)
  exp.N.fin(saved.objs$exp.N.fin)
  exp.Ns(saved.objs$exp.Ns)
  par.names(saved.objs$par.names)
  par.vals(saved.objs$par.vals)
  est.par.names(saved.objs$est.par.names)
  sim.opts.bio.scen(saved.objs$sim.opts.bio.scen)
  sim.lst(saved.objs$sim.lst)
  checks.lst(saved.objs$checks.lst)
  fit.lst(saved.objs$fit.lst)
  N.t.mat(saved.objs$N.t.mat)
  avg.phi.obs(saved.objs$avg.phi.obs)
  ns.SPs(saved.objs$ns.SPs)
  ns.POPs(saved.objs$ns.POPs)
  ns.SMPs(saved.objs$ns.SMPs)
  ns.SFPs(saved.objs$ns.SFPs)
  ns.SMPs.t(saved.objs$ns.SMPs.t)
  ns.SFPs.t(saved.objs$ns.SFPs.t)
  ns.SibPs(saved.objs$ns.SibPs)
  pns.UPs(saved.objs$pns.UPs)
  lg.gp.prbs.KP(saved.objs$lg.gp.prbs.KP)
})

# Show number of files uploaded
up.msg = reactiveVal("No files uploaded")
observeEvent(input$file, {
  up.msg(paste0("File uploaded!"))
})
output$upMsg = renderText(up.msg())
