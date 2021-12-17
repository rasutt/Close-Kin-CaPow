// POPAN negative log-likelihood function for TMB

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Declare data inputs
  DATA_INTEGER(k);
  DATA_IVECTOR(srvygaps);
  DATA_SCALAR(ncaphists);
  DATA_VECTOR(firsttab);
  DATA_VECTOR(lasttab);
  DATA_VECTOR(caps);
  DATA_VECTOR(noncaps);
  DATA_VECTOR(survives);
  
  // Declare parameter input
  PARAMETER_VECTOR(pars);
  
  // Unpack parameters and request standard error for lambda
  Type rho = pars(0);
  Type phi = pars(1);
  Type Ns = pars(2);
  vector<Type> pvec(k);
  for(int i = 0; i < k; i++) {
    pvec(i) = pars(3 + i);
  }
  Type lambda = rho + phi;
  ADREPORT(lambda);
  
  // Find lambda and phi with respect to gaps between surveys
  vector<Type> lambdagaps(k - 1);
  vector<Type> phigaps(k - 1);
  for(int i = 0; i < k - 1; i++) {
    lambdagaps(i) = pow(lambda, srvygaps(i));
    phigaps(i) = pow(phi, srvygaps(i));
  }
  
  // Find entry proportions
  vector<Type> pentvec(k);
  pentvec(0) = Type(1.0);
  Type cumlambda = Type(1.0);
  for(int t = 1; t < k; t++) {
    pentvec(t) = (lambdagaps(t - 1) - phigaps(t - 1)) * cumlambda;
    cumlambda *= lambdagaps(t - 1);
  }
  Type sumpentvec = pentvec.sum();
  for(int i = 0; i < k; i++) {
    pentvec(i) = pentvec(i) / sumpentvec;
  }

  // Find E(Nt) and specify it as a derived parameter to calculate and report
  // the standard error
  vector<Type> exp_N_t(k);
  exp_N_t(0) = Ns * pentvec(0);
  for(int i = 1; i < k; i++) {
    exp_N_t(i) = exp_N_t(i - 1) * phigaps(i - 1) + Ns * pentvec(i);
  }
  ADREPORT(exp_N_t(k - 1));

  // Find probability that animal not consequenty observed given alive in
  // population at each survey
  vector<Type> chivec(k);
  chivec(k - 1) = Type(1.0);
  for(int i = k - 2; i >= 0; i--) {
    chivec(i) = Type(1.0) - phigaps(i) +
      phigaps(i) * (Type(1.0) - pvec(i + 1)) * chivec(i + 1);
  }

  // Find probability that animal alive in population and not previously
  // observed at each survey
  vector<Type> psivec(k);
  psivec(0) = pentvec(0);
  for(int i = 1; i < k; i++) {
    psivec(i) = psivec(i - 1) * (1 - pvec(i - 1)) * phigaps(i - 1) + pentvec(i);
  }

  // Find 1 - pvec
  vector<Type> pveccomp(k);
  for(int i = 0; i < k; i++) {
    pveccomp(i) = Type(1.0) - pvec(i);
  }

  // Find log(P(unseen))
  Type logpunseen = log((pentvec * pveccomp * chivec).sum());

  // Find and return negative log likelihood
  Type nllike = -lgamma(Ns + Type(1.0)) + lgamma(Ns - ncaphists + Type(1.0)) -
    (Ns - ncaphists) * logpunseen;
  nllike = nllike -
    (firsttab * log(psivec)).sum() -
    (lasttab * log(chivec)).sum() -
    (caps * log(pvec)).sum() -
    (noncaps * log(pveccomp)).sum() -
    (survives * log(phigaps)).sum();
  return nllike;
}
