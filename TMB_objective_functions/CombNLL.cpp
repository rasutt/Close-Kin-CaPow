// Combined POPAN and close kin negative log-likelihood function for TMB

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Declare data inputs
  // POPAN
  DATA_INTEGER(k);
  DATA_IVECTOR(srvygaps);
  DATA_SCALAR(ncaphists);
  DATA_VECTOR(firsttab);
  DATA_VECTOR(lasttab);
  DATA_VECTOR(caps);
  DATA_VECTOR(noncaps);
  DATA_VECTOR(survives);
  
  // Close kin
  DATA_IVECTOR(nsPOPswtn);
  DATA_IVECTOR(nsHSPswtn);
  DATA_IVECTOR(nsPOPsbtn);
  DATA_IVECTOR(nsSPsbtn);
  DATA_INTEGER(fyear);
  DATA_IVECTOR(srvyyrs);
  DATA_IVECTOR(nscaps);
  DATA_INTEGER(alpha);
  
  // Declare parameter input
  PARAMETER_VECTOR(pars);
  
  // Unpack parameters
  Type rho = pars(0);
  Type phi = pars(1);
  Type Nfinal = pars(2);
  vector<Type> pvec(k);
  for(int i = 0; i < k; i++) {
    pvec(i) = pars(3 + i);
  } 
  Type lambda = rho + phi;

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
  pentvec.setZero();
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
  
  // Find E(Nt) and specify it as a derived parameter to calculate and report
  // the standard error
  vector<Type> exp_N_t(k);
  exp_N_t(0) = Ns * pentvec(0);
  for(int i = 1; i < k; i++) {
    exp_N_t(i) = exp_N_t(i - 1) * phigaps(i - 1) + Ns * pentvec(i);
  }
  ADREPORT(exp_N_t);

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

  // Find popan model negative log likelihood
  Type nll = -lgamma(Ns + Type(1.0)) + lgamma(Ns - ncaphists + Type(1.0)) -
    (Ns - ncaphists) * logpunseen;
  nll = nll -
    (firsttab * log(psivec)).sum() -
    (lasttab * log(chivec)).sum() -
    (caps * log(pvec)).sum() -
    (noncaps * log(pveccomp)).sum() -
    (survives * log(phigaps)).sum();

  // Parent-offspring pairs within samples
  
  // Declare variables for loop
  Type expNsurvyr, prbPOPswtn, prbHSPswtn;
  vector<Type> prbPOPswtnvec(k);
  int nonPOPsHSPswtn;
  
  // Loop over number of surveys
  for(int srvyind = 0; srvyind < k; srvyind++) {
    // Find the expected number alive at the sample year
    expNsurvyr = Nfinal / pow(lambda, fyear - srvyyrs(srvyind));
    
    // Probability of POPs within sample
    prbPOPswtn = Type(4.0) / (expNsurvyr - Type(1.0)) * phi * rho /
      (lambda - pow(phi, 2));
    
    // Store probability of POPs between samples
    prbPOPswtnvec(srvyind) = prbPOPswtn;
    
    // Probability of HSPs within sample
    prbHSPswtn = Type(4.0) * rho / (expNsurvyr - Type(1.0)) * 
      (Type(1.0) - phi / lambda) * pow(lambda / phi, alpha) * 
      (lambda * (phi + lambda) / pow(lambda - pow(phi, 2), 2) -
      Type(4.0) * rho * pow(lambda / phi, alpha) * 
      pow(phi, 2) * lambda / pow(lambda - pow(phi, 3), 2));
    
    // Find number of non-POP-non-HSPs within sample
    nonPOPsHSPswtn = nscaps(srvyind) * (nscaps(srvyind) - 1) / 2 -
      nsPOPswtn(srvyind) - nsHSPswtn(srvyind);
    
    // Add negative log likelihood from numbers of POPs and HSPs within sample
    nll = nll - nsPOPswtn(srvyind) * log(prbPOPswtn) -
      nsHSPswtn(srvyind) * log(prbHSPswtn) -
      nonPOPsHSPswtn * log(Type(1.0) - prbPOPswtn - prbHSPswtn);
  }
  
  // Self and parent-offspring pairs between samples
  
  // Declare variables for loop
  Type expNsurvyr1, expNsurvyr2, prbPOPsbrnbtn, prbPOPsbrnbtnfctr;
  Type prbPOPsbtn, prbSPsbtn;
  int srvyyr1, srvyyr2, srvygap, nonPOPsSPsbtn;
  
  // Set pair counter to zero
  int prcnt = 0;
  
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
      prbSPsbtn = pow(phi, srvygap) / expNsurvyr2;
      
      // Set probability factor term to zero
      prbPOPsbrnbtnfctr = Type(0.0);
      
      // Loop over years between surveys
      for(int gapyrind = 0; gapyrind < srvygap; gapyrind++) {
        // Find which factor term to add, not sure of the biological
        // significance lol
        if(gapyrind > srvygap - alpha - 1) {
          // Add to probability factor term
          prbPOPsbrnbtnfctr = prbPOPsbrnbtnfctr +
            pow(phi / lambda, gapyrind);
        } else {
          // Add to probability factor term
          prbPOPsbrnbtnfctr = prbPOPsbrnbtnfctr +
            pow(phi / lambda, srvygap - alpha - 1);
        }
      }
      
      // Probability of POPs where one born between samples
      prbPOPsbrnbtn = Type(2.0) * rho / lambda / expNsurvyr1 *
        prbPOPsbrnbtnfctr;
      
      // Add modified probability from POPS that existed in the first survey.
      // Strangely simplified with self-pair probability lol
      prbPOPsbtn = prbPOPsbrnbtn + prbPOPswtnvec(srvyind1) * prbSPsbtn *
        (expNsurvyr1 - 1);
      
      // Find number of non-POP-non-SPs within sample
      nonPOPsSPsbtn = nscaps(srvyind1) * nscaps(srvyind2) -
        nsSPsbtn(prcnt) - nsPOPsbtn(prcnt);
      
      // Add negative log likelihood from numbers of SPs and POPs observed
      nll = nll -
        nsSPsbtn(prcnt) * log(prbSPsbtn) -
        nsPOPsbtn(prcnt) * log(prbPOPsbtn) -
        nonPOPsSPsbtn * log(Type(1.0) - prbPOPsbtn - prbSPsbtn);
      
      // Increment pair counter.  Remember first index is zero in cpp
      prcnt++;
    }
  }
  
  // Return negative log likelihood
  return nll;
}  
