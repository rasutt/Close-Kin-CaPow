# Functions for app outputs

# Function to make data frame of parameter values for display
par.vals.df = function(par.vals, par.names, lg.nm.inds = 3:4) {
  par.df = data.frame(matrix(par.vals, nrow = 1))
  names(par.df) = par.names
  par.df[, lg.nm.inds] = as.integer(par.df[, lg.nm.inds])
  par.df
}

# Function to prepare proportion to print as percentage
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

# Function to create servers for values, biases, and errors modules
VPE.srvr <- function(id, vals, preds, errs, types, var.name) {
  moduleServer(id, function(input, output, session) {
    lapply(1:length(types), function(i) {
      output[[paste0("vals", types[i])]] = 
        renderPlot(plot.vals(vals()[[i]], preds()[[i]], var.name))
      output[[paste0("bss", types[i])]] = 
        renderTable(find.bias.srvy(errs()[[i]]))
      output[[paste0("errs", types[i])]] = 
        renderPlot(plot.errs(errs()[[i]], var.name))
    })
  })
}

