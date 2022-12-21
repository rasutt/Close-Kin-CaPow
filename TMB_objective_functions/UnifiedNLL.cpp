// Unified negative log-likelihood function for TMB

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Declare data inputs
  
  // Model type
  DATA_STRING(mdl_tp)
  
  // Study features
  DATA_INTEGER(k);
  DATA_IVECTOR(srvy_gaps);
  DATA_INTEGER(fnl_year);
  DATA_IVECTOR(srvy_yrs);
  
  // Popan model inputs
  DATA_SCALAR(n_cap_hists);
  DATA_VECTOR(first_tab);
  DATA_VECTOR(last_tab);
  DATA_VECTOR(caps);
  DATA_VECTOR(non_caps);
  DATA_VECTOR(survives);
  
  // Close-kin model inputs
  DATA_INTEGER(alpha);
  
  // True kinship model inputs
  DATA_IVECTOR(ns_SPs_btn);
  DATA_IVECTOR(ns_POPs_btn);
  DATA_IVECTOR(ns_POPs_wtn);
  // DATA_IVECTOR(ns_HSPs_wtn);
  DATA_IVECTOR(ns_caps);

  // Genopair model inputs
  DATA_MATRIX(gp_probs);
  DATA_IMATRIX(smp_yr_ind_prs);
  DATA_INTEGER(n_pairs);
  
  // Declare parameter input
  PARAMETER_VECTOR(pars);
  
  // Unpack parameters and request standard error for lambda
  Type rho = pars(0);
  Type phi = pars(1);
  Type lambda = rho + phi;
  ADREPORT(lambda);
  
  // Find superpopulation size to compare with POPAN estimates
  
  // Find lambda and phi with respect to gaps between surveys
  vector<Type> lambdagaps(k - 1);
  vector<Type> phigaps(k - 1);
  for(int i = 0; i < k - 1; i++) {
    lambdagaps(i) = pow(lambda, srvy_gaps(i));
    phigaps(i) = pow(phi, srvy_gaps(i));
  }
  
  // Find entry proportions
  vector<Type> pentvec(k);
  pentvec(0) = Type(1.0);
  Type cumlambda = Type(1.0);
  for(int i = 1; i < k; i++) {
    pentvec(i) = (lambdagaps(i - 1) - phigaps(i - 1)) * cumlambda;
    cumlambda *= lambdagaps(i - 1);
  }
  Type sumpentvec = pentvec.sum();
  for(int i = 0; i < k; i++) {
    pentvec(i) = pentvec(i) / sumpentvec;
  }
  
  // Set negative log likelihood to zero
  Type nll = Type(0.0);
  
  if(mdl_tp == "true kinship" || mdl_tp == "genopair") {
    // Find Ns as first N_t divided by first entry proportion and specify it as a
    // derived parameter to calculate and report the standard error
    Type Nfinal = pars(2);
    Type Ns = (Nfinal / cumlambda) / pentvec(0);
    ADREPORT(Ns);
    
    // Parent-offspring pair probabilities within samples
    
    // Declare variables for loop
    Type expNsurvyr, expnsSMPs;
    matrix<Type> prbPOPsmat(k, k), prbHSPsmat(k, k);
    prbPOPsmat.setZero();
    prbHSPsmat.setZero();
    
    // Lambda minus phi-squared
    Type lmb_m_ph_sq = lambda - pow(phi, 2);
    // Probability not new-born (phi over lambda)
    Type p_o_l = phi / lambda;
    // Reciprocal of probability not new-born (phi over lambda)
    Type l_o_p = lambda / phi;
    // Reciprocal of probability that an animal is mature
    Type rcl_prb_mtr = pow(l_o_p, alpha);
    // Birth rate among mature females
    Type beta = 2 * (1 - p_o_l) * rcl_prb_mtr;
    
    // Loop over surveys
    for(int srvyind = 0; srvyind < k; srvyind++) {
      // Find the expected number alive at the sample year
      expNsurvyr = Nfinal / pow(lambda, fnl_year - srvy_yrs(srvyind));
      
      // Probability of POPs within sample
      prbPOPsmat(srvyind, srvyind) = Type(2.0) / (expNsurvyr - Type(1.0)) * 
        rho * (Type(1.0) + phi) / (lambda - pow(phi, 2));
      
      // Same-mother pairs within survey years
      // exp.ns.SMPs.wtn = exp.N.s.yrs * beta * rho * phi^2 / lmb_m_ph_sq^2
      expnsSMPs = expNsurvyr * beta * rho * pow(phi, 2) / pow(lmb_m_ph_sq, 2);
      
      // // Same-father pairs within survey years, split into same and different
      // // birth years as used separately later
      // exp.ns.SFPs.diff.b.yrs.wtn = phi * exp.ns.SMPs.wtn
      // exp.ns.SFPs.same.b.yr.wtn = beta^2 * phi^(alpha + 1) / 4 * 
      //   (exp.N.s.yrs / (lambda^(alpha - 1) * lmb_m_ph_sq) - 1 / (1 - phi^2))
      // exp.ns.SFPs.wtn = exp.ns.SFPs.diff.b.yrs.wtn + exp.ns.SFPs.same.b.yr.wtn
      // 
      // // Full-sibling pairs within survey years, constant over time but gets
      // // repeated by cbind when returned
      // exp.ns.FSPs.wtn = 2 * beta * rclprbmtr * rho * phi^4 / 
      //   (lambda * (lambda - phi^3) * (1 - phi^2))
      //   
      // // Half-sibling pairs within survey years
      // exp.ns.HSPs.wtn = exp.ns.SMPs.wtn + exp.ns.SFPs.wtn - 2 * exp.ns.FSPs.wtn
        
        
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
      srvyyr1 = srvy_yrs(srvyind1);
      
      // Loop over surveys with greater indices than first
      for(int srvyind2 = srvyind1 + 1; srvyind2 < k; srvyind2++) {
        // Find second survey year
        srvyyr2 = srvy_yrs(srvyind2);
        
        // Find gap between surveys.  Remember not necessarily consecutive
        srvygap = srvyyr2 - srvyyr1;
        
        // Find the expected numbers alive in survey years
        expNsurvyr1 = Nfinal / pow(lambda, fnl_year - srvyyr1);
        expNsurvyr2 = Nfinal / pow(lambda, fnl_year - srvyyr2);
        
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
    
    // If fitting genopair model
    if(mdl_tp == "true kinship") {
      // Numbers of unrelated pairs
      // int nonPOPsHSPswtn;
      int nonPOPswtn, nonPOPsSPsbtn;
      
      // Parent-offspring pairs within samples
      
      // Loop over number of surveys
      for(int srvyind = 0; srvyind < k; srvyind++) {
        // Find number of non-POP-non-HSPs within sample
        // nonPOPsHSPswtn = ns_caps(srvyind) * (ns_caps(srvyind) - 1) / 2 -
        //   ns_POPs_wtn(srvyind) - ns_HSPs_wtn(srvyind);
        nonPOPswtn = ns_caps(srvyind) * (ns_caps(srvyind) - 1) / 2 -
          ns_POPs_wtn(srvyind);
        
        // Get parent-offspring pair probability
        prbPOP = prbPOPsmat(srvyind, srvyind);
        
        // Add negative log likelihood from numbers of POPs and HSPs within sample
        // nll = nll - ns_POPs_wtn(srvyind) * log(prbPOPswtn) -
        //   ns_HSPs_wtn(srvyind) * log(prbHSPswtn) -
        //   nonPOPsHSPswtn * log(Type(1.0) - prbPOPswtn - prbHSPswtn);
        nll = nll - ns_POPs_wtn(srvyind) * log(prbPOP) -
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
          nonPOPsSPsbtn = ns_caps(srvyind1) * ns_caps(srvyind2) -
            ns_SPs_btn(prcnt) - ns_POPs_btn(prcnt);
          
          // Add negative log likelihood from numbers of SPs and POPs observed.
          // Omitting multinomial coefficient as only adds a constant w.r.t. the
          // parameters.
          nll = nll -
            ns_SPs_btn(prcnt) * log(prbSP) -
            ns_POPs_btn(prcnt) * log(prbPOP) -
            nonPOPsSPsbtn * log(Type(1.0) - prbPOP - prbSP);
          
          // Increment pair counter.  Remember first index is zero in cpp
          prcnt++;
        }
      }
    }
    
    // If fitting genopair model
    if(mdl_tp == "genopair") {
      // Temporary variables for genopairs
      int inds1, inds2;
      
      // Loop over genopairs
      for(int gpind = 0; gpind < n_pairs; gpind++) {
        // Get sample-year indices and kinship probabilities
        inds1 = smp_yr_ind_prs(gpind, 0);
        inds2 = smp_yr_ind_prs(gpind, 1);
        prbPOP = prbPOPsmat(inds1, inds2);
        prbSP = prbSPsmat(inds1, inds2);
        
        // Add negative log likelihood from genopair probabilities given kinships
        // and kinship probabilities in terms of parameters.
        nll = nll - log(
          (Type(1.0) - prbPOP - prbSP) * gp_probs(gpind, 0) +
            prbPOP * gp_probs(gpind, 1) + prbSP * gp_probs(gpind, 2)
        );
      }
    }
  } 
  
  if(mdl_tp == "popan") {
    Type Ns = pars(2);
    vector<Type> pvec(k);
    for(int i = 0; i < k; i++) {
      pvec(i) = pars(3 + i);
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
      psivec(i) = 
        psivec(i - 1) * (1 - pvec(i - 1)) * phigaps(i - 1) + pentvec(i);
    }
    
    // Find 1 - pvec
    vector<Type> pveccomp(k);
    for(int i = 0; i < k; i++) {
      pveccomp(i) = Type(1.0) - pvec(i);
    }
    
    // Find log(P(unseen))
    Type logpunseen = log((pentvec * pveccomp * chivec).sum());
    
    // Find and return negative log likelihood
    nll = -lgamma(Ns + Type(1.0)) + lgamma(Ns - n_cap_hists + Type(1.0)) -
      (Ns - n_cap_hists) * logpunseen;
    nll = nll -
      (first_tab * log(psivec)).sum() -
      (last_tab * log(chivec)).sum() -
      (caps * log(pvec)).sum() -
      (non_caps * log(pveccomp)).sum() -
      (survives * log(phigaps)).sum();
  }
  
  // Return negative log likelihood
  return nll;
}  
