#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
 compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");
   
  //arrays to be filled by results from ODE Solver
  Vector deltaCDM_array(n_k * n_x);
  Vector deltaB_array(n_k * n_x);
  Vector vCDM_array(n_k * n_x);
  Vector vB_array(n_k * n_x);
  Vector Phi_array(n_k * n_x);
  Vector Theta0_array(n_k * n_x);
  Vector Theta1_array(n_k * n_x);
  Vector Theta2_array(n_k * n_x);
  Vector Theta3_array(n_k * n_x);
  Vector Theta4_array(n_k * n_x);
  Vector Theta5_array(n_k * n_x);
  Vector Theta6_array(n_k * n_x);
  Vector Theta7_array(n_k * n_x);

  //Psi-array to be filled later
  Vector Psi_array(n_k * n_x);
  
  //Set up x-array to integrate over
  Vector x_array =  Utils::linspace(x_start, x_end, n_x);
  
  //Set up k_array with logarithmic spacing
  Vector k_array;
  k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  for(int ik = 0; ik < n_k; ik++){
     k_array[ik] = exp(k_array[ik]);
 }

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){
     
    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to and show it
    double x_end_tight = get_tight_coupling_time(k);
    int index_end_tc;
    double eps = 1e-2;
    for(int iix = 0; iix < n_x; iix++){
       if(x_array[iix] > x_end_tight-eps and x_array[iix] < x_end_tight+eps ){
           index_end_tc = iix;
           std::cout << "x_end_tight: " << x_end_tight << " with index: " << index_end_tc  << std::endl;
        }
    }
    
    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODESolver y_tight_coupling_ode;
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };
    Vector x_array_tight_coupling = Utils::linspace(x_start, x_end_tight, index_end_tc+1); //make sure we have n_x points all together (x_end_tight occurs two times)
    y_tight_coupling_ode.solve(dydx_tight_coupling, x_array_tight_coupling, y_tight_coupling_ini);

    auto deltaCDM_tight_coupling = y_tight_coupling_ode.get_data_by_component(0);
    auto deltab_tight_coupling = y_tight_coupling_ode.get_data_by_component(1);
    auto vCDM_tight_coupling = y_tight_coupling_ode.get_data_by_component(2);
    auto vb_tight_coupling = y_tight_coupling_ode.get_data_by_component(3);
    auto Phi_tight_coupling = y_tight_coupling_ode.get_data_by_component(4);
    auto Theta0_tight_coupling = y_tight_coupling_ode.get_data_by_component(5);
    auto Theta1_tight_coupling = y_tight_coupling_ode.get_data_by_component(6);

    //Create vector of initial values for integration after tight coupling regime
    int last_ind = index_end_tc;
    Vector y_tight_coupling{deltaCDM_tight_coupling[last_ind], deltab_tight_coupling[last_ind], vCDM_tight_coupling[last_ind], vb_tight_coupling[last_ind], Phi_tight_coupling[last_ind], Theta0_tight_coupling[last_ind], Theta1_tight_coupling[last_ind]};
    
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODESolver y_full_ode;
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };
    Vector x_array_after = Utils::linspace(x_end_tight, x_end, n_x-last_ind); 
    y_full_ode.solve(dydx_full, x_array_after, y_full_ini);

    auto deltaCDM_after = y_full_ode.get_data_by_component(0);
    auto deltab_after = y_full_ode.get_data_by_component(1);
    auto vCDM_after = y_full_ode.get_data_by_component(2);
    auto vb_after = y_full_ode.get_data_by_component(3);
    auto Phi_after = y_full_ode.get_data_by_component(4);
    auto Theta0_after = y_full_ode.get_data_by_component(5);
    auto Theta1_after = y_full_ode.get_data_by_component(6);
    auto Theta2_after = y_full_ode.get_data_by_component(7);
    auto Theta3_after = y_full_ode.get_data_by_component(8);
    auto Theta4_after = y_full_ode.get_data_by_component(9);
    auto Theta5_after = y_full_ode.get_data_by_component(10);
    auto Theta6_after = y_full_ode.get_data_by_component(11);
    auto Theta7_after = y_full_ode.get_data_by_component(12);
    
    //Fill arrays corresponding to all quantities
    int index;
    for(int ix = 0; ix < n_x; ix++){
       index = ix + n_x * ik;
       if(ix < last_ind){ //we are in tight coupling regime
           deltaCDM_array[index] = deltaCDM_tight_coupling[ix];
           deltaB_array[index] =  deltab_tight_coupling[ix] ;
           vCDM_array[index] = vCDM_tight_coupling[ix];
           vB_array[index] = vb_tight_coupling[ix];
           Phi_array[index] = Phi_tight_coupling[ix];
           Theta0_array[index] = Theta0_tight_coupling[ix];
           Theta1_array[index] = Theta1_tight_coupling[ix];
           //all higher multipoles are expressed directly in terms of the lower ones
           double x = x_array[ix];
           double c = Constants.c;
           double Hp;
           double dtaudx;
           dtaudx = rec->dtaudx_of_x(x);
           Hp = cosmo->Hp_of_x(x);
           Theta2_array[index] = -(20.*c*k*Theta1_array[index])/(45.*Hp*dtaudx);
           Theta3_array[index] = -3.*c*k*Theta2_array[index]/(7.*Hp*dtaudx);
           Theta4_array[index] = -4.*c*k*Theta3_array[index]/(9.*Hp*dtaudx);
           Theta5_array[index] = -5.*c*k*Theta4_array[index]/(11.*Hp*dtaudx);
           Theta6_array[index] = -6.*c*k*Theta5_array[index]/(13.*Hp*dtaudx);
           Theta7_array[index] = -7.*c*k*Theta6_array[index]/(15.*Hp*dtaudx);
       } else { //after tight coupling regime
           deltaCDM_array[index] = deltaCDM_after[ix-last_ind];
           deltaB_array[index] =  deltab_after[ix-last_ind] ;
           vCDM_array[index] = vCDM_after[ix-last_ind];
           vB_array[index] = vb_after[ix-last_ind];
           Phi_array[index] = Phi_after[ix-last_ind];
           Theta0_array[index] = Theta0_after[ix-last_ind];
           Theta1_array[index] = Theta1_after[ix-last_ind];
           Theta2_array[index] = Theta2_after[ix-last_ind]; 
           Theta3_array[index] = Theta3_after[ix-last_ind]; 
           Theta4_array[index] = Theta4_after[ix-last_ind]; 
           Theta5_array[index] = Theta5_after[ix-last_ind]; 
           Theta6_array[index] = Theta6_after[ix-last_ind]; 
           Theta7_array[index] = Theta7_after[ix-last_ind]; 
           }
     }

  
  }// end loop over wavenumbers
  Utils::EndTiming("integrateperturbation");
  
  //Create splines of variables calculated by ODE Solver
  delta_cdm_spline.create(x_array,k_array,deltaCDM_array, "delta_cdm_spline");
  delta_b_spline.create(x_array,k_array,deltaB_array, "delta_b_spline");
  v_cdm_spline.create(x_array,k_array,vCDM_array, "v_cdm_spline");
  v_b_spline.create(x_array,k_array,vB_array, "v_b_spline");
  Phi_spline.create(x_array,k_array,Phi_array, "v_cdm_spline");
  Theta0_spline.create(x_array,k_array, Theta0_array, "Theta0_spline");
  Theta1_spline.create(x_array,k_array, Theta1_array, "Theta1_spline");
  Theta2_spline.create(x_array,k_array, Theta2_array, "Theta2_spline");
  Theta3_spline.create(x_array,k_array, Theta3_array, "Theta3_spline");
  Theta4_spline.create(x_array,k_array, Theta4_array, "Theta4_spline");
  Theta5_spline.create(x_array,k_array, Theta5_array, "Theta5_spline");
  Theta6_spline.create(x_array,k_array, Theta6_array, "Theta6_spline");
  Theta7_spline.create(x_array,k_array, Theta7_array, "Theta7_spline");
  
  //Calculate Psi
  int index;
  for(int ix = 0; ix < n_x; ix++){
       double x = x_array[ix];
   for(int ik = 0; ik < n_k; ik++){
      double k = k_array[ik];
      index = ix + n_x * ik;
      if(x == x_start){
        Psi_array[index] = -2./3.;
      } else {
      double c = Constants.c;
      double H0; //Hubble parameter today (1/s)
      double Hp;
      double dtaudx;
      dtaudx = rec->dtaudx_of_x(x);
      double OmegaR_today;
      H0 = cosmo->get_H0();
      Hp = cosmo->Hp_of_x(x);
      OmegaR_today = cosmo->get_OmegaR(0.);
      double a = exp(x);
      double Phi= Phi_array[index];
      double Theta1 = Theta1_array[index]; 
      double Theta2;
      Theta2 = Theta2_array[index]; 
      double Psi = -Phi-((12.*H0*H0*OmegaR_today*Theta2)/(c*c*k*k*a*a));
      Psi_array[index] = Psi;
      }//end else
   }//end for ik

  }//end for ix
   
  Psi_spline.create(x_array,k_array, Psi_array, "Psi_spline");
} 
//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  double c = Constants.c;
  double Hp;
  double dtaudx;
  Hp = cosmo->Hp_of_x(x);
  dtaudx = rec->dtaudx_of_x(x);

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double &Theta0        = y_tc[Constants.ind_start_theta_tc];
  double &Theta1        = y_tc[Constants.ind_start_theta_tc+1];

  
  //Scalar quantities (Gravitational potential, baryons and CDM)
  double Psi = -2./3.;
  Phi = -Psi; 
  delta_cdm = (-3./2.)*Psi;
  delta_b = (-3./2.)*Psi;
  v_cdm = (-c*k/(2.*Hp))*Psi;
  v_b = (-c*k/(2.*Hp))*Psi;
  
  
  //Photon temperature perturbations (Theta_ell)
  Theta0 = -0.5*Psi;
  Theta1 = (c*k/(6.*Hp))*Psi;
  
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  
  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double &Theta0_tc        = y_tc[Constants.ind_start_theta_tc];
  const double &Theta1_tc        = y_tc[Constants.ind_start_theta_tc+1];
  

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double &Theta0           = y[Constants.ind_start_theta_tc];
  double &Theta1           = y[Constants.ind_start_theta_tc+1];
  double &Theta2           = y[Constants.ind_start_theta_tc+2];
  //Theta3, Theta4, Theta5, Theta6, Theta7 will be set directly, without references
  
  double c = Constants.c;
  double Hp;
  double dtaudx;
  Hp = cosmo->Hp_of_x(x);
  dtaudx = rec->dtaudx_of_x(x);

  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm  = v_cdm_tc;
  v_b = v_b_tc;
  Phi = Phi_tc;
  Theta0 = Theta0_tc;
  Theta1 = Theta1_tc;
  Theta2 = -20.*c*k*Theta1/(45.*Hp*dtaudx);
  double Theta_l, Theta_lminus1;
  for(int l = 3; l < Constants.n_ell_theta ; l++){
    Theta_lminus1 = y[Constants.ind_start_theta_tc+l-1];
    Theta_l = -l*c*k*Theta_lminus1/((2.*l + 1)*Hp*dtaudx);
    y[Constants.ind_start_theta_tc+l] = Theta_l;
  }
  
  
  return y;
}

//====================================================
// The time when tight coupling ends
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end;
  int counter = 0;
  double z_recombination = 1700.;
  double a_recombination = 1./(z_recombination + 1.);
  double x_recombination = log(a_recombination);
  Vector x_array =  Utils::linspace(x_start, x_end, n_x);
  double Hp;
  double dtaudx;
  double c = Constants.c;
  double minimum_of_two;
  for(int i = 0; i < n_x; i++){
     Hp = cosmo->Hp_of_x(x_array[i]);
     dtaudx = rec->dtaudx_of_x(x_array[i]);
     if(c*k/Hp < 1.){
        minimum_of_two = c*k/Hp;
     } else {
        minimum_of_two = 1.;
     }
     if(abs(dtaudx) < 10.*minimum_of_two){
        if(counter == 0){
           x_tight_coupling_end = x_array[i];
        }
        counter += 1;
     }
  
  } //end for loop
  //std::cout << "x_rec: " << x_recombination << std::endl;
  if(x_tight_coupling_end > x_recombination){ //if the condition results in time after recombination
     x_tight_coupling_end = x_recombination;
  }
  
 
  return x_tight_coupling_end;
}

//====================================================
// After integrating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //Set up x-array
  Vector x_array =  Utils::linspace(x_start, x_end, n_x);

  //Set up k_array with logarithmic spacing
  Vector k_array;
  k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  for(int ik = 0; ik < n_k; ik++){
  k_array[ik] = exp(k_array[ik]);
 }
 
  // Make storage for the source function (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  double c = Constants.c;  
  // Compute source function
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
   for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];
      const int index = ix + n_x * ik;
      const double Hp       = cosmo->Hp_of_x(x);
      const double Hp_prime = cosmo->dHpdx_of_x(x);
      const double Hp_prime_prime = cosmo->ddHpddx_of_x(x);
      const double tau      = rec->tau_of_x(x);
      const double g_tilde = rec->g_tilde_of_x(x);
      const double g_tilde_prime = rec->dgdx_tilde_of_x(x);  
      const double g_tilde_prime_prime = rec->ddgddx_tilde_of_x(x); 
      const double Theta0 = get_Theta(x,k,0);
      const double Pi = get_Theta(x,k,2); 
      double Pi_prime = Theta2_spline.deriv_x(x,k);
      double Pi_prime_prime = Theta2_spline.deriv_xx(x,k);  
      const double v_b = get_v_b(x,k);
      double v_b_prime = v_b_spline.deriv_x(x,k);
      double Hpgv_prime = Hp_prime*g_tilde*v_b + Hp*g_tilde_prime*v_b + Hp*g_tilde*v_b_prime;
      double last_term1 = (Hp_prime*Hp_prime + Hp*Hp_prime_prime)*g_tilde*Pi;
      double last_term2 = 3.*Hp*Hp_prime*(g_tilde_prime*Pi + g_tilde*Pi_prime); 
      double last_term3 =  Hp*Hp*(g_tilde_prime_prime*Pi + 2.*g_tilde_prime*Pi_prime + g_tilde*Pi_prime_prime);
      double last_term = last_term1 + last_term2 + last_term3;
      const double Psi = get_Psi(x,k);
      double Psi_prime = Psi_spline.deriv_x(x,k);
      double Phi_prime = Phi_spline.deriv_x(x,k);
      double source_f1 = g_tilde*(Theta0 + Psi + Pi/4.);
      double source_f2 = exp(-tau)*(Psi_prime - Phi_prime);
      double source_f3 = -Hpgv_prime/(c*k);
      double source_f4 = 3.*last_term/(4.*k*k*c*c);
      double source_f = source_f1 + source_f2 + source_f3 + source_f4;
      
      ST_array[index] = source_f;
      
    }
  }

  // Spline the source function
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  
  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double &Theta_0 = y[Constants.ind_start_theta_tc];
  const double &Theta_1 = y[Constants.ind_start_theta_tc+1];
 

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double &dThetadx_0 = dydx[Constants.ind_start_theta_tc];
  double &dThetadx_1 = dydx[Constants.ind_start_theta_tc+1];
  
  double c = Constants.c;
  double H0; //Hubble parameter today (1/s)
  double Hp;
  double dHpdx;
  double dtaudx;
  double ddtauddx;
  double OmegaR_today, OmegaCDM_today, OmegaB_today;
  H0 = cosmo->get_H0();
  Hp = cosmo->Hp_of_x(x);
  dHpdx = cosmo->dHpdx_of_x(x);
  dtaudx = rec->dtaudx_of_x(x);
  ddtauddx = rec->ddtauddx_of_x(x);
  OmegaR_today = cosmo->get_OmegaR(0.);
  OmegaCDM_today = cosmo->get_OmegaCDM(0.);
  OmegaB_today = cosmo->get_OmegaB(0.);
  double a = exp(x);
  
  double Theta_2 = -(20.*c*k*Theta_1)/(45.*Hp*dtaudx);
  double Psi = -Phi-((12.*H0*H0*OmegaR_today*Theta_2)/(c*c*k*k*a*a));
  
  double sum1 = OmegaCDM_today*(1./a)*delta_cdm + OmegaB_today*(1./a)*delta_b + 4.*OmegaR_today*(1./(a*a))*Theta_0;
  dPhidx = Psi - ((c*c*k*k*Phi)/(3.*Hp*Hp)) + ((H0*H0*sum1)/(2.*Hp*Hp));
  ddelta_cdmdx  = c*k*v_cdm/Hp -3.*dPhidx;
  ddelta_bdx = c*k*v_b/Hp -3.*dPhidx;
  dv_cdmdx = -v_cdm - c*k*Psi/Hp;
  dThetadx_0 = -c*k*Theta_1/Hp - dPhidx;
  double R = 4.*OmegaR_today/(3.*OmegaB_today*a);
  double q, q_numerator, q_denominator;
  q_numerator = -(dtaudx*(1.-R) + ddtauddx*(1.+R))*(3.*Theta_1 + v_b) - c*k*Psi/Hp + (1. - dHpdx/Hp)*(c*k/Hp)*(-Theta_0 + 2.*Theta_2) - c*k*dThetadx_0/Hp;
  q_denominator = (1.+R)*dtaudx + dHpdx/Hp - 1.;
  q = q_numerator/q_denominator;
  dv_bdx = (1./(1. + R))*(-v_b - c*k*Psi/Hp + R*(q+ (c*k/Hp)*(-Theta_0 + 2.*Theta_2) - c*k*Psi/Hp));
  dThetadx_1 = (q-dv_bdx)/3.;
 
  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double &Theta0           = y[Constants.ind_start_theta];
  const double &Theta1           = y[Constants.ind_start_theta+1];
  const double &Theta2           = y[Constants.ind_start_theta+2];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double &dThetadx0        = dydx[Constants.ind_start_theta];
  double &dThetadx1        = dydx[Constants.ind_start_theta+1];
 

  double c = Constants.c;
  double H0; //Hubble parameter today (1/s)
  double Hp;
  double dHpdx;
  double dtaudx;
  double ddtauddx;
  double OmegaR_today, OmegaCDM_today, OmegaB_today;
  double eta;
  H0 = cosmo->get_H0();
  Hp = cosmo->Hp_of_x(x);
  dHpdx = cosmo->dHpdx_of_x(x);
  dtaudx = rec->dtaudx_of_x(x);
  ddtauddx = rec->ddtauddx_of_x(x);
  OmegaR_today = cosmo->get_OmegaR(0.);
  OmegaCDM_today = cosmo->get_OmegaCDM(0.);
  OmegaB_today = cosmo->get_OmegaB(0.);
  eta = cosmo->eta_of_x(x);
  double a = exp(x);

  double Psi = -Phi-((12.*H0*H0*OmegaR_today*Theta2)/(c*c*k*k*a*a));
  double sum1 = OmegaCDM_today*(1./a)*delta_cdm + OmegaB_today*(1./a)*delta_b + 4.*OmegaR_today*(1./(a*a))*Theta0;
  dPhidx = Psi - ((c*c*k*k*Phi)/(3.*Hp*Hp)) + ((H0*H0*sum1)/(2.*Hp*Hp));
  ddelta_cdmdx  = c*k*v_cdm/Hp -3.*dPhidx;
  ddelta_bdx = c*k*v_b/Hp -3.*dPhidx;
  dv_cdmdx = -v_cdm - c*k*Psi/Hp;
  double R = 4.*OmegaR_today/(3.*OmegaB_today*a);
  dv_bdx = -v_b - c*k*Psi/Hp + dtaudx*R*(3.*Theta1 + v_b);

  dThetadx0 = -c*k*Theta1/Hp - dPhidx;
  dThetadx1 = c*k*Theta0/(3.*Hp) - 2.*c*k*Theta2/(3.*Hp) + c*k*Psi/(3.*Hp) + dtaudx*(Theta1 + v_b/3.);
  int l_max = Constants.n_ell_theta - 1;
  double dThetadxl, Thetalminus1, Thetalplus1, Thetal;
  for(int l = 2; l < l_max ; l++){
     Thetalminus1 = y[Constants.ind_start_theta+l-1];
     Thetalplus1 = y[Constants.ind_start_theta+l+1];
     Thetal = y[Constants.ind_start_theta+l];
     if ( l == 2){
        dThetadxl = l*c*k*Thetalminus1/(Hp*(2.*l + 1.)) - (l + 1.)*c*k*Thetalplus1/(Hp*(2.*l + 1.)) + dtaudx*(Thetal - Theta2/10.);
     } else {
        dThetadxl = l*c*k*Thetalminus1/(Hp*(2.*l + 1.)) - (l + 1.)*c*k*Thetalplus1/(Hp*(2.*l + 1.)) + dtaudx*Thetal;
     }
     dydx[Constants.ind_start_theta+l] = dThetadxl;
  }
  Thetal = y[Constants.ind_start_theta+l_max];
  dydx[Constants.ind_start_theta+l_max] = c*k*y[Constants.ind_start_theta+l_max-1]/Hp -c*(l_max+1)*Thetal/(Hp*eta) + dtaudx*Thetal;

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
  
}

double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}

double Perturbations::get_Theta(const double x, const double k, const int ell) const{
   if(ell == 0){
      return Theta0_spline(x,k);
   } else if (ell == 1) {
     return Theta1_spline(x,k);
   } else if (ell == 2) {
     return Theta2_spline(x,k);
   } else if (ell == 3) {
     return Theta3_spline(x,k);
   } else if (ell == 4) {
     return Theta4_spline(x,k);
   } else if (ell == 5) {
     return Theta5_spline(x,k);
   } else if (ell == 6) {
     return Theta6_spline(x,k);
   } else if (ell == 7) {
     return Theta7_spline(x,k);
   }
}



//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 1500;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    //double arg = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k) << " ";
    fp << get_v_cdm(x,k) << " ";
    fp << get_v_b(x,k) << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Psi(x,k)       << " ";
    //fp << get_Source_T(x,k)  << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

