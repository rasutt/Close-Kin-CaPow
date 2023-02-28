# Functions for app outputs

# Function to make data frame of parameter values for display
par.vals.df = function(par.vals, par.names, lng.nmbr.inds = 3:4) {
  par.df = data.frame(matrix(par.vals, nrow = 1))
  names(par.df) = par.names
  par.df[, lng.nmbr.inds] = as.integer(par.df[, lng.nmbr.inds])
  par.df
}

# Function to format parameters implied by inputs for output
frmt.pars.impld = function(lambda, beta, N.init, exp.N.lst, exp.Ns) {
  par.vals.df(
    c(lambda, beta, N.init, exp.N.lst, exp.Ns),
    c("Population growth rate", 
      "Birthrate among mature females", "Initial population size", 
      "Expected final population size", "Expected superpopulation size"),
    3:5
  )
}

# Function to prepare proportion to print as percentage with 1 D.P.
perc = function(prpn) paste0(round(prpn * 100, 1), "%")

# Function to find biases over all surveys from array of proportional errors for
# multiple estimators
find.bias = function(errs) {
  df = data.frame(matrix(perc(colMeans(errs, dims = 2)), nrow = 1))
  names(df) = dimnames(errs)[["kp.type"]]
  df
}

# Function to find biases in each survey from matrix of proportional errors for
# a single estimator
find.bias.srvy = function(errs) {
  df = data.frame(matrix(perc(colMeans(errs)), nrow = 1))
  names(df) = colnames(errs)
  df
}

# Function to plot simulated and predicted values ----
plot.vals = function(vals, preds, par.name) {
  boxplot(
    vals, main = par.name, xlab = names(dimnames(vals))[2], show.names = T,
    ylab = "Observed values"
  )
  points(preds, col = 'red', lwd = 2)
  points(colMeans(vals), col = 'blue', lwd = 2)
  legend(
    "topleft", col = c(2, 4), lty = 1, lwd = 2,
    legend = c("Predicted value", "Average value"),
  )
}

# Function to plot proportional differences between simulated and predicted
# values ----
plot.errs = function(errs, var.name) {
  boxplot(
    errs, main = var.name, xlab = names(dimnames(errs))[2], show.names = T,
    ylab = "Proportional errors"
  )
  abline(h = 0, col = 'red')
  abline(h = mean(errs), col = 'blue')
  legend(
    "topleft", col = c(2, 4), lty = 1,
    legend = c("Estimated error (zero)", "Average error"),
  )
}

# Function to make outputs for values, biases, and errors
VBE.ui <- function(id, types, descs) {
  ns <- NS(id)
  tagList(lapply(1:length(types), function(i) {
    type = types[i]
    list(
      h3(type),
      p(descs[i]),
      h4("Values"),
      plotOutput(ns(paste0("vals", type))),
      h4("Biases"),
      tableOutput(ns(paste0("bss", type))),
      h4("Errors"),
      plotOutput(ns(paste0("errs", type)))
    )
  }))
}

# Function to make outputs for titles, descriptions, values, biases, and errors
TDVBE.ui <- function(
    id, ttl_dsc, types = wtn_btn_headings, typ_dscs = wtn_btn_descs
) {
  tagList(h2(kp.nms[id]), p(ttl_dsc), VBE.ui(id, types, typ_dscs))
}

# Function to make kin-pair tab-panel
KP.tab.ui = function(id) {
  tabPanel(
    title = kp.nms[id],
    value = paste0(id, ".tab"),
    TDVBE.ui(id, rglr.kp.dscs[id])
  )
}

# Functions to create servers for values, biases, and errors modules, in general
# case, and for regular estimates within and between survey-years
VPE.srvr <- function(id, vals, preds, errs, types = wtn_btn_headings) {
  moduleServer(id, function(input, output, session) {
    lapply(1:length(types), function(i) {
      output[[paste0("vals", types[i])]] = 
        renderPlot(plot.vals(vals()[[i]], preds()[[i]], kp.nms[id]))
      output[[paste0("bss", types[i])]] = 
        renderTable(find.bias.srvy(errs()[[i]]))
      output[[paste0("errs", types[i])]] = 
        renderPlot(plot.errs(errs()[[i]], kp.nms[id]))
    })
  })
}
VPE.srvr.rglr <- function(id, vals, preds.lst, errs, types = wtn_btn_headings) {
  VPE.srvr(id, vals, reactive(preds.lst()[[id]]), errs, types)
}

