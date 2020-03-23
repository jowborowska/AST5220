#include"RecombinationHistory.h"
# define M_PI           3.14159265358979323846  /* pi */

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    double h,
    double OmegaB,
    double TCMB,
    BackgroundCosmology *cosmo, 
    double Yp) :
  h(h),
  OmegaB(OmegaB),
  TCMB(TCMB),
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  Vector x_array;
  Vector Xe_array(npts_rec_arrays);
  Vector ne_array(npts_rec_arrays);
  x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts_rec_arrays);
  
  // initial value for Peebles integration
  double Xe_ini;
  double counter;
  counter = 0;
  bool saha_regime = true;
 
  for(int i = 0; i < npts_rec_arrays; i++){
    
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    //Find Saha prediction for halfway recombination
    if(Xe_current > 0.49 and Xe_current < 0.51){
       std:: cout << "Saha half-way" << " " << Xe_current << "  " << x_array[i] << '\n' ; //0.502406  -7.23077

     };

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit){
      saha_regime = false;
      counter = counter += 1.;
      if(counter == 1.0){
      std:: cout << "Peebles regime starts" << " " << Xe_current << " " << x_array[i] << '\n'; //first one at 0.988998 -7.36781
      };
    };

    if(saha_regime){
      Xe_array[i] = Xe_current;
      ne_array[i] = ne_current;
      
    } else {
       
      Xe_ini = Xe_array[i-1];
      
      // The Peebles ODE equation
      ODESolver Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      Vector Xe_ic{Xe_ini};
      //array of previous and present x  used for integration
      Vector x_array_Peebles;
      x_array_Peebles = {x_array[i-1],x_array[i]};
 
      Xe_ode.solve(dXedx, x_array_Peebles, Xe_ic);
      auto result1 = Xe_ode.get_data_by_component(0);
      Xe_array[i] = result1[1];
      ne_array[i] = Xe_array[i]*compute_nb(x_array[i]);

    } //end else - Peebles
    }//end for loop
  Vector log_Xe_array(npts_rec_arrays);
  Vector log_ne_array(npts_rec_arrays);
  for(int i = 0; i < npts_rec_arrays; i++){
      log_Xe_array[i] = log(Xe_array[i]);
      log_ne_array[i] =  log(ne_array[i]);

  };
  log_Xe_of_x_spline.create(x_array,log_Xe_array, "Xe");
  log_ne_of_x_spline.create(x_array,log_ne_array, "ne");


  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;
  
  double Tb; //baryon temperature
  Tb = TCMB/a;
  double nb; //baryon density (assume all baryons are protons)
  nb = compute_nb(x);
  
  double A; double B; double C; //quaadratic eq. coefficients
  double D; 
  A = -1.0;
  C = (1.0/nb)*pow(m_e*Tb*k_b/(2.0*M_PI*hbar*hbar), 1.5)*exp(-epsilon_0/(k_b*Tb));
  B = -C;
  D = -4.0*A*C/(B*B);
  double Xe1; double Xe2; //solutions of quadratic equation
 
 if(abs(D) < 1e-8){ 
     //Xe1 = (-B + B*(1.0+0.5*D))/(2.0*A);
     //Xe2 = (-B - B*(1.0+0.5*D))/(2.0*A);
     Xe1 = 1.;
     Xe2 = 1.;
 } else {
    //Xe1 = (-B + B*sqrt(1.0+D))/(2.0*A);
    //Xe2 = (-B - B*sqrt(1.0+D))/(2.0*A);  
    Xe1 = (-C + C*sqrt(1.0+4.0/C))/(2.0);
  };
  
  // Electron fraction and number density
  double Xe = Xe1; //Xe2 gives negative number
  double ne = nb*Xe;
  
  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;  

 
  double Tb; //baryon temperature
  Tb = TCMB/a;
  double nb; //baryon density (assume all baryons are protons)
  nb = compute_nb(x);

  double alpha = 1.0/137.0; //fine-structure constant
  double H = cosmo->H_of_x(x);
  double phi_2 = 0.448*log(epsilon_0/(Tb*k_b));
  double alpha_2 = phi_2*sqrt(epsilon_0/(k_b*Tb))*(64.0*M_PI/sqrt(27.0*M_PI))*alpha*alpha*hbar*hbar/(m_e*m_e*c);
  double beta = alpha_2*pow(m_e*k_b*Tb/(2.0*M_PI*hbar*hbar),1.5)*exp(-epsilon_0/(k_b*Tb));
  double beta_2;
  if(3.0*epsilon_0/(4.0*k_b*Tb) > 200){
     beta_2 = 0.0;
  } else {
  beta_2 = beta*exp(3.0*epsilon_0/(4.0*k_b*Tb));
  };
  double n1_s = (1.0 - X_e)*nb;
  double Lambda_alpha = H*27.0*epsilon_0*epsilon_0*epsilon_0/(64.0*M_PI*M_PI*n1_s*hbar*hbar*hbar*c*c*c);
  
  double Cr = (lambda_2s1s + Lambda_alpha)/(lambda_2s1s + Lambda_alpha + beta_2);

  dXedx[0] = (Cr/H)*(beta*(1.0 - X_e) -nb*alpha_2*X_e*X_e);
  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector tau_array(npts);
  Vector dtau_array(npts);
  
  ODESolver tau_ode;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    const double sigma_T     = Constants.sigma_T;
    const double c           = Constants.c;
    double ne = ne_of_x(x);
    double H = cosmo->H_of_x(x);

    dtaudx[0] = -ne*sigma_T*c/H;

    return GSL_SUCCESS;
  };
   
   //random initial value
   double tau_ini = 100.;
   Vector tau_ic{tau_ini};
   tau_ode.solve(dtaudx, x_array, tau_ic);
   auto result2 = tau_ode.get_data_by_component(0);
   double tau_today = result2[npts-1];

   const double sigma_T     = Constants.sigma_T;
   const double c           = Constants.c;
    for(int i = 0; i < npts; i++){
        tau_array[i] = result2[i]  - tau_today;
        double ne = ne_of_x(x_array[i]);
        double H = cosmo->H_of_x(x_array[i]);
        dtau_array[i] = -ne*sigma_T*c/H;
  }
  
  tau_of_x_spline.create(x_array,tau_array, "tau");
  dtau_of_x_spline.create(x_array,dtau_array, "dtau");
  

  Vector g_tilde_array(npts);
  for(int i = 0; i < npts; i++){
  g_tilde_array[i] = -dtaudx_of_x(x_array[i])*exp(-tau_of_x(x_array[i]));
  }
  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g");

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return dtau_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return dtau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
   return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
   return g_tilde_of_x_spline.deriv_xx(x); //wanted to use the same stratedy for ddtau, but gave some weird stuff
}

double RecombinationHistory::Xe_of_x(double x) const{
   return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

double RecombinationHistory::compute_nb(double x) const{
  const double a           = exp(x);

  // Physical constants in SI units
  const double H0_over_h   = Constants.H0_over_h;  
  const double G           = Constants.G;
  const double m_H         = Constants.m_H;
  
  double H0; //Hubble parameter today (1/s)
  H0 = H0_over_h*h; 
  double rho_crit; // Critical density today
  rho_crit = 3.0*H0*H0/(8.0*M_PI*G);
  double nb; //baryon density (assume all baryons are protons)
  nb = OmegaB*rho_crit/(m_H*a*a*a);
   return nb;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 3000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

