CombNLL <- function(pars) {
  # print(pars)
  
  # Unpack parameters and transform for survey gaps
  rho <- pars[1]
  phi <- pars[2]
  N.fin <- pars[3]
  lambda <- rho + phi
  lambda.gaps <- lambda^srvy.gaps
  phi.gaps <- phi^srvy.gaps
  
  # Find Ns as first N_t divided by first entry proportion
  Ns <- (N.fin / prod(lambda.gaps)) / pent_func(lambda.gaps, phi.gaps)[1]
  
  # Find and return sum of POPAN and close kin negative log likelihoods
  popan_nll(c(pars[1:2], Ns, pars[4:(3 + k)])) + CloseKinNLL(pars[1:3])
}