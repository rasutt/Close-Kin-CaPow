// Genopair negative log-likelihood function for TMB

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Declare data inputs
  DATA_MATRIX(gpprobs);
  DATA_IMATRIX(sampyrinds);
  DATA_INTEGER(npairs);
  DATA_INTEGER(k);
  DATA_IVECTOR(srvygaps);
  DATA_INTEGER(fyear);
  DATA_IVECTOR(srvyyrs);
  DATA_INTEGER(alpha);
  
  // Declare parameter input
  PARAMETER_VECTOR(pars);
  
  // Unpack parameters and request standard error for lambda
  Type rho = pars(0);
  Type phi = pars(1);
  Type Nfinal = pars(2);
  Type lambda = rho + phi;
  ADREPORT(lambda);
  
  // Find superpopulation size to compare with POPAN estimates
  
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
  
  // Find Ns as first N_t divided by first entry proportion and specify it as a
  // derived parameter to calculate and report the standard error
  Type Ns = (Nfinal / cumlambda) / pentvec(0);
  ADREPORT(Ns);
  
  // Set negative log likelihood to zero
  Type nll = Type(0.0);
  
  // Parent-offspring pair probabilities within samples
  
  // Declare variables for loop
  Type expNsurvyr;
  matrix<Type> prbPOPsmat(k, k);
  prbPOPsmat.setZero();

  // Loop over surveys
  for(int srvyind = 0; srvyind < k; srvyind++) {
    // Find the expected number alive at the sample year
    expNsurvyr = Nfinal / pow(lambda, fyear - srvyyrs(srvyind));
    
    // Probability of POPs within sample
    prbPOPsmat(srvyind, srvyind) = Type(2.0) / (expNsurvyr - Type(1.0)) * rho * 
      (Type(1.0) + phi) / (lambda - pow(phi, 2));
  }
  
  // Self and parent-offspring pairs between samples

  // Declare variables for loop
  Type expNsurvyr1, expNsurvyr2, prbPOPsbrnbtn;
  matrix<Type> prbSPsmat(k, k);
  prbSPsmat.setZero();
  int srvyyr1, srvyyr2, srvygap;

  // Loop over all but last survey
  for(int srvyind1 = 0; srvyind1 < k - 1; srvyind1++) {
    // Find first survey year
    srvyyr1 = srvyyrs(srvyind1);

    // Loop over surveys with greater indices than first
    for(int srvyind2 = srvyind1 + 1; srvyind2 < k; srvyind2++) {
      // Find second survey year
      srvyyr2 = srvyyrs(srvyind2);

      // Find gap between surveys.  Remember not necessarily consecutive
      srvygap = srvyyr2 - srvyyr1;

      // Find the expected numbers alive in survey years
      expNsurvyr1 = Nfinal / pow(lambda, fyear - srvyyr1);
      expNsurvyr2 = Nfinal / pow(lambda, fyear - srvyyr2);

      // Probability of SPs between samples
      prbSPsmat(srvyind1, srvyind2) = pow(phi, srvygap) / expNsurvyr2;
      
      // New probability of POPs between samples
      if(srvyyr1 + alpha < srvyyr2) {
        prbPOPsbrnbtn = (srvyyr2 - (srvyyr1 + alpha)) * 
          pow(lambda / phi, srvyyr1 + alpha);
      } else {
        prbPOPsbrnbtn = Type(0.0);
      }
      prbPOPsmat(srvyind1, srvyind2) = prbPOPsmat(srvyind1, srvyind1) * 
        (expNsurvyr1 - 1) / expNsurvyr2 * pow(phi, srvygap) + 
        Type(2.0) / expNsurvyr1 * (Type(1.0) - phi / lambda) * 
        pow(phi / lambda, srvyyr2) *
        ((pow(lambda / phi, srvyyr1 + Type(1.0)) - 
        pow(lambda / phi, std::min(srvyyr1 + alpha, srvyyr2) + Type(1.0))) / 
        (Type(1.0) - lambda / phi) + prbPOPsbrnbtn);
    }
  }
  ADREPORT(prbPOPsmat);
  ADREPORT(prbSPsmat);
  
  int inds1, inds2;
  Type prbPOP, prbSP;
  
  // Loop over genopairs
  for(int gpind = 0; gpind < npairs; gpind++) {
    // Get sample-year indices and kinship probabilities
    inds1 = sampyrinds(gpind, 0);
    inds2 = sampyrinds(gpind, 1);
    prbPOP = prbPOPsmat(inds1, inds2);
    prbSP = prbSPsmat(inds1, inds2);
    
    // Add negative log likelihood from genopair probabilities given kinships
    // and kinship probabilities in terms of parameters.
    nll = nll - log(
      (Type(1.0) - prbPOP - prbSP) * gpprobs(gpind, 0) +
        prbPOP * gpprobs(gpind, 1) + prbSP * gpprobs(gpind, 2)
    );
  }
  
  // Return negative log likelihood
  return nll;
}  
