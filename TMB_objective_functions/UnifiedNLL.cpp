// Unified negative log-likelihood function for TMB

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Declare data inputs
  
  // Model type
  DATA_STRING(mdl_tp);
  
  // Study features
  DATA_INTEGER(k);
  DATA_VECTOR(srvy_gaps);
  DATA_SCALAR(fnl_year);
  DATA_VECTOR(srvy_yrs);
  
  // Popan model inputs
  DATA_SCALAR(n_cap_hists);
  DATA_VECTOR(first_tab);
  DATA_VECTOR(last_tab);
  DATA_VECTOR(caps);
  DATA_VECTOR(non_caps);
  DATA_VECTOR(survives);
  
  // Close-kin model inputs
  DATA_SCALAR(alpha);
  DATA_IVECTOR(knshp_st_bool);
  
  // True kinship model inputs
  DATA_IVECTOR(ns_SPs_btn);
  DATA_IVECTOR(ns_POPs_wtn);
  DATA_IVECTOR(ns_POPs_btn);
  DATA_IVECTOR(ns_HSPs_wtn);
  DATA_IVECTOR(ns_HSPs_btn);
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
  Type sum_p_ent = p_ent.sum();
  for(int i = 0; i < k; i++) {
    p_ent(i) = p_ent(i) / sum_p_ent;
  }
  
  // Set negative log likelihood to zero
  Type nll = Type(0.0);
  
  // If fitting true kinship or genopair model
  if(mdl_tp == "true kinship" || mdl_tp == "genopair") {
    // Find super-population size to compare with popan models, first N_t
    // divided by first entry proportion, and request standard error
    Type N_fnl = pars(2);
    Type Ns = (N_fnl / cum_lmbd) / p_ent(0);
    ADREPORT(Ns);
    
    // Kinship probabilities. Using vectorised code as simpler to avoid
    // recomputing terms for each element in expressions, though slightly less
    // tidy.
    
    // Unpack kinship set variables
    bool incld_SPs, incld_POPs, incld_HSPs;
    incld_SPs = knshp_st_bool(0) == 1;
    incld_POPs = knshp_st_bool(1) == 1;
    incld_HSPs = knshp_st_bool(2) == 1;
    
    // Declare variables
    Type lmb_m_ph_sq, p_o_l, l_o_p, rcl_prb_mtr, beta, s_yr_1, s_yr_2, s_gap, 
      p_t_s_gap, exp_N_s_yr_1, exp_N_s_yr_2, prb_POPs_brn_btn;
    vector<Type> exp_N_s_yrs(k);
    matrix<Type> exp_ns_APs(k, k), prbs_SPs(k, k), prbs_POPs(k, k), 
      exp_ns_SMPs(k, k), exp_ns_SFPs_diff_b_yrs(k, k), 
      exp_ns_SFPs_same_b_yr(k, k), exp_ns_SFPs(k, k), exp_ns_FSPs(k, k), 
      exp_ns_HSPs(k, k), prbs_HSPs(k, k);

    // Lambda minus phi-squared
    lmb_m_ph_sq = lambda - pow(phi, 2);
    // Probability not new-born (phi over lambda)
    p_o_l = phi / lambda;
    // Reciprocal of probability not new-born (phi over lambda)
    l_o_p = lambda / phi;
    // Reciprocal of probability that an animal is mature
    rcl_prb_mtr = pow(l_o_p, alpha);
    // Birth rate among mature females
    beta = 2 * (1 - p_o_l) * rcl_prb_mtr;
    
    // Expected population sizes in survey years
    exp_N_s_yrs = 
      N_fnl / pow(lambda, vector<Type>::Constant(k, fnl_year) - srvy_yrs);
    
    // Expected total number of pairs in survey years
    exp_ns_APs = (exp_N_s_yrs * 
      (exp_N_s_yrs - Type(1.0)) / Type(2.0)).matrix().asDiagonal();
    
    if(incld_SPs) {
      // Initialise probabilities of SPs
      prbs_SPs.setZero();
    }
    
    if(incld_POPs) {
      // Probability of POPs within samples
      prbs_POPs = (Type(2.0) / (exp_N_s_yrs - Type(1.0)) * 
        rho * (Type(1.0) + phi) / lmb_m_ph_sq).matrix().asDiagonal();
    }
    
    if(incld_HSPs) {
      // Same-mother pairs within survey years
      exp_ns_SMPs = (exp_N_s_yrs * beta * rho * 
        pow(phi, 2) / pow(lmb_m_ph_sq, 2)).matrix().asDiagonal();
      
      // Same-father pairs within survey years, split into same and different
      // birth years as used separately later
      exp_ns_SFPs_same_b_yr = 
        (pow(beta, 2) * pow(phi, alpha + Type(1.0)) / Type(4.0) *
        (exp_N_s_yrs / (pow(lambda, alpha - Type(1.0)) * lmb_m_ph_sq) - 
        Type(1.0) / (Type(1.0) - pow(phi, 2)))).matrix().asDiagonal();
      exp_ns_SFPs_diff_b_yrs = phi * exp_ns_SMPs;
      exp_ns_SFPs = exp_ns_SFPs_diff_b_yrs + exp_ns_SFPs_same_b_yr;
      
      // Full-sibling pairs within survey years, constant over time
      exp_ns_FSPs = matrix<Type>::Constant(
        k, 1, 
        Type(2.0) * beta * rcl_prb_mtr * rho * pow(phi, 4) /
          (lambda * (lambda - pow(phi, 3)) * (Type(1.0) - pow(phi, 2)))
      ).asDiagonal();
      
      // Half-sibling pairs within survey years
      exp_ns_HSPs.diagonal() = exp_ns_SMPs.diagonal() + exp_ns_SFPs.diagonal() - 
        Type(2.0) * exp_ns_FSPs.diagonal();
      prbs_HSPs.setZero();
      prbs_HSPs.diagonal() = (exp_ns_HSPs.diagonal().array() / 
        exp_ns_APs.diagonal().array()).matrix();
    }
    
    // Self and parent-offspring pairs between samples
    
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
        exp_N_s_yr_1 = exp_N_s_yrs(s_ind_1);
        exp_N_s_yr_2 = exp_N_s_yrs(s_ind_2);
        
        // Phi to the power of the survey-gap
        p_t_s_gap = pow(phi, s_gap);

        // All pairs between pairs of surveys
        exp_ns_APs(s_ind_1, s_ind_2) = exp_N_s_yr_1 * exp_N_s_yr_2;
        
        // Probability of SPs between samples
        if(incld_SPs) {
          prbs_SPs(s_ind_1, s_ind_2) = p_t_s_gap / exp_N_s_yr_2;
        }
        
        // Probability of POPs between samples
        if(incld_POPs) {
          // POPs born between samples
          if(s_yr_1 + alpha < s_yr_2) {
            prb_POPs_brn_btn = (s_yr_2 - (s_yr_1 + alpha)) * 
              pow(lambda / phi, s_yr_1 + alpha);
          } else {
            prb_POPs_brn_btn = Type(0.0);
          }
          
          // All POPs
          prbs_POPs(s_ind_1, s_ind_2) = prbs_POPs(s_ind_1, s_ind_1) * 
            (exp_N_s_yr_1 - 1) / exp_N_s_yr_2 * p_t_s_gap + 
            Type(2.0) / exp_N_s_yr_1 * (Type(1.0) - phi / lambda) * 
            pow(phi / lambda, s_yr_2) *
            ((pow(lambda / phi, s_yr_1 + Type(1.0)) - 
            pow(lambda / phi, std::min(s_yr_1 + alpha, s_yr_2) + Type(1.0))) / 
            (Type(1.0) - lambda / phi) + prb_POPs_brn_btn);
        }
        
        // Half-sibling pairs
        if(incld_HSPs) {
          // Same-mother pairs between survey years
          exp_ns_SMPs(s_ind_1, s_ind_2) = Type(2.0) * 
            exp_ns_SMPs(s_ind_1, s_ind_1) * p_t_s_gap +
            s_gap * exp_N_s_yr_2 * beta * (1 - p_o_l) *
            lambda / lmb_m_ph_sq * pow(p_o_l, s_gap);
          
          // Same-father pairs between survey years
          exp_ns_SFPs(s_ind_1, s_ind_2) = phi * exp_ns_SMPs(s_ind_1, s_ind_2) +
            Type(2.0) * p_t_s_gap * exp_ns_SFPs_same_b_yr(s_ind_1, s_ind_1);
          
          // Full-sibling pairs between survey years, note the predicted number
          // within surveys is constant
          exp_ns_FSPs(s_ind_1, s_ind_2) = 2 * exp_ns_FSPs(s_ind_1, s_ind_1) * 
            p_t_s_gap + 2 * beta * (1 - p_o_l) * rcl_prb_mtr * 
            pow(phi, s_yr_1 + s_yr_2 + 1) *
            (pow(p_o_l, s_yr_1 + 1) - pow(p_o_l, s_yr_2 + 1)) /
              (pow(pow(phi, 3) / lambda, s_yr_1) * (1 - pow(phi, 3) / lambda) * 
                (1 - p_o_l));
          
          // Half-sibs
          prbs_HSPs(s_ind_1, s_ind_2) = (exp_ns_SMPs(s_ind_1, s_ind_2) + 
            exp_ns_SFPs(s_ind_1, s_ind_2) - 
            Type(2.0) * exp_ns_FSPs(s_ind_1, s_ind_2)) / 
            exp_ns_APs(s_ind_1, s_ind_2);
        }
      }
    }
    
    // Request values and derivatives for kinship probabilities, need values
    // separately to check for true parameter values
    REPORT(prbs_SPs);
    REPORT(prbs_POPs);
    REPORT(prbs_HSPs);
    ADREPORT(prbs_SPs);
    ADREPORT(prbs_POPs);
    ADREPORT(prbs_HSPs);
    
    // Temporary variables for kinpair probabilities
    Type prb_SP, prb_POP, prb_HSP, prb_UP;
    
    // If fitting true kinship model
    if(mdl_tp == "true kinship") {
      // Numbers of "unrelated" pairs, not related in another way
      int ns_UPs_wtn, ns_UPs_btn;
      
      if(incld_POPs | incld_HSPs) {
        // Pairs within samples
      
        // Loop over number of surveys
        for(int s_ind = 0; s_ind < k; s_ind++) {
          // Find number of "unrelated pairs" within sample, updated below
          ns_UPs_wtn = ns_caps(s_ind) * (ns_caps(s_ind) - 1) / 2;
          
          prb_UP = Type(1.0);
          
          if(incld_POPs) {
            ns_UPs_wtn = ns_UPs_wtn - ns_POPs_wtn(s_ind);
            prb_POP = prbs_POPs(s_ind, s_ind);
            nll = nll - ns_POPs_wtn(s_ind) * log(prb_POP);
            prb_UP = prb_UP - prb_POP;
          }
          
          if(incld_HSPs) {
            ns_UPs_wtn = ns_UPs_wtn - ns_HSPs_wtn(s_ind);
            prb_HSP = prbs_HSPs(s_ind, s_ind);
            nll = nll - ns_HSPs_wtn(s_ind) * log(prb_HSP);
            prb_UP = prb_UP - prb_HSP;
          }
          
          nll = nll - ns_UPs_wtn * log(prb_UP);
        }
      }
      
      // Pairs between samples
      
      // Set pair counter to zero
      int pr_cnt = 0;
      
      // Loop over all but last survey
      for(int s_ind_1 = 0; s_ind_1 < k - 1; s_ind_1++) {
        // Loop over surveys with greater indices than first
        for(int s_ind_2 = s_ind_1 + 1; s_ind_2 < k; s_ind_2++) {
          // Find number of non-POP-non-SPs within sample
          ns_UPs_btn = ns_caps(s_ind_1) * ns_caps(s_ind_2);

          prb_UP = Type(1.0);
          
          if(incld_SPs) {
            // Probability of SPs between samples
            prb_SP = prbs_SPs(s_ind_1, s_ind_2);
            ns_UPs_btn = ns_UPs_btn - ns_SPs_btn(pr_cnt);
            nll = nll - ns_SPs_btn(pr_cnt) * log(prb_SP);
            prb_UP = prb_UP - prb_SP;
          }
          if(incld_POPs) {
            // Get parent-offspring pair probability
            prb_POP = prbs_POPs(s_ind_1, s_ind_2);
            ns_UPs_btn = ns_UPs_btn - ns_POPs_btn(pr_cnt);
            nll = nll - ns_POPs_btn(pr_cnt) * log(prb_POP);
            prb_UP = prb_UP - prb_POP;
          }
          if(incld_HSPs) {
            // Get parent-offspring pair probability
            prb_HSP = prbs_HSPs(s_ind_1, s_ind_2);
            ns_UPs_btn = ns_UPs_btn - ns_HSPs_btn(pr_cnt);
            nll = nll - ns_HSPs_btn(pr_cnt) * log(prb_HSP);
            prb_UP = prb_UP - prb_HSP;
          }
          
          // Add negative log likelihood from numbers of UPs observed. Omitting
          // multinomial coefficients as only adds a constant w.r.t. the
          // parameters.
          nll = nll - ns_UPs_btn * log(prb_UP);
          
          // Increment pair counter.  Remember first index is zero in cpp
          pr_cnt++;
        }
      }
    }
    
    // If fitting genopair model
    if(mdl_tp == "genopair") {
      // Temporary variables for genopairs
      int smp_yr_ind_1, smp_yr_ind_2;
      
      Type gp_prb;
      
      // Find indices of kinships in genopair probabilities matrix
      int knshp_ind = 1, HSP_ind = 0, POP_ind = 0, SP_ind = 0;
      if(incld_HSPs) {
        HSP_ind = knshp_ind;
        knshp_ind++;
      }
      if(incld_POPs) {
        POP_ind = knshp_ind;
        knshp_ind++;
      }
      if(incld_SPs) {
        SP_ind = knshp_ind;
      }
      
      // Loop over genopairs
      for(int gpind = 0; gpind < n_pairs; gpind++) {
        // Get sample-year indices
        smp_yr_ind_1 = smp_yr_ind_prs(gpind, 0);
        smp_yr_ind_2 = smp_yr_ind_prs(gpind, 1);
        
        // Reset genopair and unrelated pair probabilities
        gp_prb = Type(0.0);
        prb_UP = Type(1.0);
        
        // Add genopair probabilities over kinships and update unrelated pair
        // probability
        if(incld_SPs) {
          prb_SP = prbs_SPs(smp_yr_ind_1, smp_yr_ind_2);
          gp_prb = gp_prb + prb_SP * gp_probs(gpind, SP_ind);
          prb_UP = prb_UP - prb_SP;
        }
        if(incld_POPs) {
          prb_POP = prbs_POPs(smp_yr_ind_1, smp_yr_ind_2);
          gp_prb = gp_prb + prb_POP * gp_probs(gpind, POP_ind);
          prb_UP = prb_UP - prb_POP;
        }
        if(incld_HSPs) {
          prb_HSP = prbs_HSPs(smp_yr_ind_1, smp_yr_ind_2);
          gp_prb = gp_prb + prb_HSP * gp_probs(gpind, HSP_ind);
          prb_UP = prb_UP - prb_HSP;
        }
        
        // Add negative log likelihood from genopair probabilities given kinships
        // and kinship probabilities in terms of parameters.
        nll = nll - log(gp_prb + prb_UP * gp_probs(gpind, 0));
      }
    }
    
    // Return negative log likelihood
    return nll;
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
    
    // Find and return negative log likelihood
    return -lgamma(Ns + Type(1.0)) + lgamma(Ns - n_cap_hists + Type(1.0)) -
      (Ns - n_cap_hists) * lg_p_unseen -
      (first_tab * log(psi)).sum() -
      (last_tab * log(chi)).sum() -
      (caps * log(p)).sum() -
      (non_caps * log(p_comp)).sum() -
      (survives * log(phi_gaps)).sum();
  }
  
  return nll;
}  
