// Unified negative log-likelihood function for TMB

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Declare data inputs
  // Study features
  DATA_INTEGER(k);
  DATA_IVECTOR(srvygaps);
  DATA_INTEGER(fyear);
  DATA_IVECTOR(srvyyrs);
  DATA_INTEGER(alpha);
  
  // // Popan model inputs
  // DATA_SCALAR(ncaphists);
  // DATA_VECTOR(firsttab);
  // DATA_VECTOR(lasttab);
  // DATA_VECTOR(caps);
  // DATA_VECTOR(noncaps);
  // DATA_VECTOR(survives);
  
  // True kinship model inputs
  DATA_IVECTOR(nsSPsbtn);
  DATA_IVECTOR(nsPOPsbtn);
  DATA_IVECTOR(nsPOPswtn);
  // DATA_IVECTOR(nsHSPswtn);
  DATA_IVECTOR(nscaps);

  // Genopair model inputs
  DATA_MATRIX(gpprobs);
  DATA_IMATRIX(sampyrinds);
  DATA_INTEGER(npairs);
  
  // Model type
  DATA_STRING(mdltp)
  
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
  
  // Temporary variables for kinpair probabilities
  Type prbPOP, prbSP;
  
  // Set negative log likelihood to zero
  Type nll = Type(0.0);
  
  // If fitting genopair model
  if(mdltp == "true kinship") {
    // Numbers of unrelated pairs
    // int nonPOPsHSPswtn;
    int nonPOPswtn, nonPOPsSPsbtn;
    
    // Parent-offspring pairs within samples
    
    // Loop over number of surveys
    for(int srvyind = 0; srvyind < k; srvyind++) {
      // Find number of non-POP-non-HSPs within sample
      // nonPOPsHSPswtn = nscaps(srvyind) * (nscaps(srvyind) - 1) / 2 -
      //   nsPOPswtn(srvyind) - nsHSPswtn(srvyind);
      nonPOPswtn = nscaps(srvyind) * (nscaps(srvyind) - 1) / 2 -
        nsPOPswtn(srvyind);
      
      // Get parent-offspring pair probability
      prbPOP = prbPOPsmat(srvyind, srvyind);
      
      // Add negative log likelihood from numbers of POPs and HSPs within sample
      // nll = nll - nsPOPswtn(srvyind) * log(prbPOPswtn) -
      //   nsHSPswtn(srvyind) * log(prbHSPswtn) -
      //   nonPOPsHSPswtn * log(Type(1.0) - prbPOPswtn - prbHSPswtn);
      nll = nll - nsPOPswtn(srvyind) * log(prbPOP) -
        nonPOPswtn * log(Type(1.0) - prbPOP);
    }
    
    // Self and parent-offspring pairs between samples
    
    // Set pair counter to zero
    int prcnt = 0;
    
    // Loop over all but last survey
    for(int srvyind1 = 0; srvyind1 < k - 1; srvyind1++) {
      // Loop over surveys with greater indices than first
      for(int srvyind2 = srvyind1 + 1; srvyind2 < k; srvyind2++) {
        // Probability of SPs between samples
        prbSP = prbSPsmat(srvyind1, srvyind2);
        
        // Get parent-offspring pair probability
        prbPOP = prbPOPsmat(srvyind1, srvyind2);
        
        // Find number of non-POP-non-SPs within sample
        nonPOPsSPsbtn = nscaps(srvyind1) * nscaps(srvyind2) -
          nsSPsbtn(prcnt) - nsPOPsbtn(prcnt);
        
        // Add negative log likelihood from numbers of SPs and POPs observed.
        // Omitting multinomial coefficient as only adds a constant w.r.t. the
        // parameters.
        nll = nll -
          nsSPsbtn(prcnt) * log(prbSP) -
          nsPOPsbtn(prcnt) * log(prbPOP) -
          nonPOPsSPsbtn * log(Type(1.0) - prbPOP - prbSP);
        
        // Increment pair counter.  Remember first index is zero in cpp
        prcnt++;
      }
    }
  }
  
  // If fitting genopair model
  if(mdltp == "genopair") {
    // Temporary variables for genopairs
    int inds1, inds2;

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
  }
  
  // Return negative log likelihood
  return nll;
}  
