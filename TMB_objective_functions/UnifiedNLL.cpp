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
  vector<Type> lmbd_gaps(k - 1);
  vector<Type> phi_gaps(k - 1);
  for(int i = 0; i < k - 1; i++) {
    lmbd_gaps(i) = pow(lambda, srvy_gaps(i));
    phi_gaps(i) = pow(phi, srvy_gaps(i));
  }
  
  // Find entry proportions
  vector<Type> p_ent(k);
  p_ent(0) = Type(1.0);
  Type cum_lmbd = Type(1.0);
  for(int i = 1; i < k; i++) {
    p_ent(i) = (lmbd_gaps(i - 1) - phi_gaps(i - 1)) * cum_lmbd;
    cum_lmbd *= lmbd_gaps(i - 1);
  }
  Type sump_ent = p_ent.sum();
  for(int i = 0; i < k; i++) {
    p_ent(i) = p_ent(i) / sump_ent;
  }
  
  // Set negative log likelihood to zero
  Type nll = Type(0.0);
  
  if(mdl_tp == "true kinship" || mdl_tp == "genopair") {
    // Find Ns as first N_t divided by first entry proportion and specify it as a
    // derived parameter to calculate and report the standard error
    Type N_fnl = pars(2);
    Type Ns = (N_fnl / cum_lmbd) / p_ent(0);
    ADREPORT(Ns);
    
    // Parent-offspring pair probabilities within samples
    
    // Declare variables for loop
    Type exp_N_s_yr, exp_ns_SMPs;
    matrix<Type> prbs_POPs(k, k), prbs_HSPs(k, k);
    prbs_POPs.setZero();
    prbs_HSPs.setZero();
    
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
    for(int s_ind = 0; s_ind < k; s_ind++) {
      // Expected population sizes in survey years
      exp_N_s_yr = N_fnl / pow(lambda, fnl_year - srvy_yrs(s_ind));

      // Probability of POPs within sample
      prbs_POPs(s_ind, s_ind) = Type(2.0) / (exp_N_s_yr - Type(1.0)) * 
        rho * (Type(1.0) + phi) / (lambda - pow(phi, 2));
      
      // Same-mother pairs within survey years
      // exp.ns.SMPs.wtn = exp.N.s.yrs * beta * rho * phi^2 / lmb_m_ph_sq^2
      exp_ns_SMPs = exp_N_s_yr * beta * rho * 
        pow(phi, 2) / pow(lmb_m_ph_sq, 2);
      
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
    Type exp_N_s_yr_1, exp_N_s_yr_2, prb_POPs_brn_btn;
    matrix<Type> prbs_SPs(k, k);
    prbs_SPs.setZero();
    int s_yr_1, s_yr_2, s_gap;
    
    // Loop over all but last survey
    for(int s_ind_1 = 0; s_ind_1 < k - 1; s_ind_1++) {
      // Find first survey year
      s_yr_1 = srvy_yrs(s_ind_1);
      
      // Loop over surveys with greater indices than first
      for(int s_ind_2 = s_ind_1 + 1; s_ind_2 < k; s_ind_2++) {
        // Find second survey year
        s_yr_2 = srvy_yrs(s_ind_2);
        
        // Find gap between surveys.  Not necessarily consecutive surveys so
        // different from survey gaps passed in as data.
        s_gap = s_yr_2 - s_yr_1;
        
        // Find the expected numbers alive in survey years
        exp_N_s_yr_1 = N_fnl / pow(lambda, fnl_year - s_yr_1);
        exp_N_s_yr_2 = N_fnl / pow(lambda, fnl_year - s_yr_2);
        
        // Probability of SPs between samples
        prbs_SPs(s_ind_1, s_ind_2) = pow(phi, s_gap) / exp_N_s_yr_2;
        
        // New probability of POPs between samples
        if(s_yr_1 + alpha < s_yr_2) {
          prb_POPs_brn_btn = (s_yr_2 - (s_yr_1 + alpha)) * 
            pow(lambda / phi, s_yr_1 + alpha);
        } else {
          prb_POPs_brn_btn = Type(0.0);
        }
        prbs_POPs(s_ind_1, s_ind_2) = prbs_POPs(s_ind_1, s_ind_1) * 
          (exp_N_s_yr_1 - 1) / exp_N_s_yr_2 * pow(phi, s_gap) + 
          Type(2.0) / exp_N_s_yr_1 * (Type(1.0) - phi / lambda) * 
          pow(phi / lambda, s_yr_2) *
          ((pow(lambda / phi, s_yr_1 + Type(1.0)) - 
          pow(lambda / phi, std::min(s_yr_1 + alpha, s_yr_2) + Type(1.0))) / 
          (Type(1.0) - lambda / phi) + prb_POPs_brn_btn);
      }
    }
    ADREPORT(prbs_POPs);
    ADREPORT(prbs_SPs);
    
    // Temporary variables for kinpair probabilities
    Type prb_POP, prb_SP;
    
    // If fitting true kinship model
    if(mdl_tp == "true kinship") {
      // Numbers of "unrelated" pairs
      int ns_UPs_wtn, ns_UPs_btn;
      
      // Parent-offspring pairs within samples
      
      // Loop over number of surveys
      for(int s_ind = 0; s_ind < k; s_ind++) {
        // Find number of non-POP-non-HSPs within sample
        // ns_UPs_wtn = ns_caps(s_ind) * (ns_caps(s_ind) - 1) / 2 -
        //   ns_POPs_wtn(s_ind) - ns_HSPs_wtn(s_ind);
        ns_UPs_wtn = ns_caps(s_ind) * (ns_caps(s_ind) - 1) / 2 -
          ns_POPs_wtn(s_ind);
        
        // Get parent-offspring pair probability
        prb_POP = prbs_POPs(s_ind, s_ind);
        
        // Add negative log likelihood from numbers of POPs and HSPs within sample
        // nll = nll - ns_POPs_wtn(s_ind) * log(prb_POPswtn) -
        //   ns_HSPs_wtn(s_ind) * log(prbHSPswtn) -
        //   ns_UPs_wtn * log(Type(1.0) - prb_POPswtn - prbHSPswtn);
        nll = nll - ns_POPs_wtn(s_ind) * log(prb_POP) -
          ns_UPs_wtn * log(Type(1.0) - prb_POP);
      }
      
      // Self and parent-offspring pairs between samples
      
      // Set pair counter to zero
      int pr_cnt = 0;
      
      // Loop over all but last survey
      for(int s_ind_1 = 0; s_ind_1 < k - 1; s_ind_1++) {
        // Loop over surveys with greater indices than first
        for(int s_ind_2 = s_ind_1 + 1; s_ind_2 < k; s_ind_2++) {
          // Probability of SPs between samples
          prb_SP = prbs_SPs(s_ind_1, s_ind_2);
          
          // Get parent-offspring pair probability
          prb_POP = prbs_POPs(s_ind_1, s_ind_2);
          
          // Find number of non-POP-non-SPs within sample
          ns_UPs_btn = ns_caps(s_ind_1) * ns_caps(s_ind_2) -
            ns_SPs_btn(pr_cnt) - ns_POPs_btn(pr_cnt);
          
          // Add negative log likelihood from numbers of SPs and POPs observed.
          // Omitting multinomial coefficient as only adds a constant w.r.t. the
          // parameters.
          nll = nll -
            ns_SPs_btn(pr_cnt) * log(prb_SP) -
            ns_POPs_btn(pr_cnt) * log(prb_POP) -
            ns_UPs_btn * log(Type(1.0) - prb_POP - prb_SP);
          
          // Increment pair counter.  Remember first index is zero in cpp
          pr_cnt++;
        }
      }
    }
    
    // If fitting genopair model
    if(mdl_tp == "genopair") {
      // Temporary variables for genopairs
      int smp_yr_ind_1, smp_yr_ind_2;
      
      // Loop over genopairs
      for(int gpind = 0; gpind < n_pairs; gpind++) {
        // Get sample-year indices and kinship probabilities
        smp_yr_ind_1 = smp_yr_ind_prs(gpind, 0);
        smp_yr_ind_2 = smp_yr_ind_prs(gpind, 1);
        prb_POP = prbs_POPs(smp_yr_ind_1, smp_yr_ind_2);
        prb_SP = prbs_SPs(smp_yr_ind_1, smp_yr_ind_2);
        
        // Add negative log likelihood from genopair probabilities given kinships
        // and kinship probabilities in terms of parameters.
        nll = nll - log(
          (Type(1.0) - prb_POP - prb_SP) * gp_probs(gpind, 0) +
            prb_POP * gp_probs(gpind, 1) + prb_SP * gp_probs(gpind, 2)
        );
      }
    }
  } 
  
  if(mdl_tp == "popan") {
    // Superpopulation size
    Type Ns = pars(2);
    
    // Capture probabilities and their complements
    vector<Type> p(k), p_comp(k);
    for(int i = 0; i < k; i++) {
      p(i) = pars(3 + i);
      p_comp(i) = Type(1.0) - p(i);
    }
    
    // Find probability that animal not consequently observed given alive in
    // population at each survey
    vector<Type> chi(k);
    chi(k - 1) = Type(1.0);
    for(int i = k - 2; i >= 0; i--) {
      chi(i) = Type(1.0) - phi_gaps(i) +
        phi_gaps(i) * (Type(1.0) - p(i + 1)) * chi(i + 1);
    }
    
    // Find probability that animal alive in population and not previously
    // observed at each survey, and expected population sizes in survey years
    // (probably request SEs for all years later, so leave in)
    vector<Type> psi(k), exp_N_s_yrs(k);
    psi(0) = p_ent(0);
    exp_N_s_yrs(0) = Ns * p_ent(0);
    for(int i = 1; i < k; i++) {
      psi(i) = psi(i - 1) * (1 - p(i - 1)) * phi_gaps(i - 1) + p_ent(i);
      exp_N_s_yrs(i) = exp_N_s_yrs(i - 1) * phi_gaps(i - 1) + Ns * p_ent(i);
    }
    
    // Request standard error for expected population size in final year
    ADREPORT(exp_N_s_yrs(k - 1));
    
    // Find log(P(unseen))
    Type lg_p_unseen = log((p_ent * p_comp * chi).sum());
    
    // Find negative log likelihood
    nll = -lgamma(Ns + Type(1.0)) + lgamma(Ns - n_cap_hists + Type(1.0)) -
      (Ns - n_cap_hists) * lg_p_unseen -
      (first_tab * log(psi)).sum() -
      (last_tab * log(chi)).sum() -
      (caps * log(p)).sum() -
      (non_caps * log(p_comp)).sum() -
      (survives * log(phi_gaps)).sum();
  }
  
  // Return negative log likelihood
  return nll;
}  
