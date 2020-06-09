#include"PowerSpectrum.h"
# define M_PI           3.14159265358979323846  /* pi */

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){
  Vector k_array;
  k_array = Utils::linspace(Constants.k_min, Constants.k_max, n_k);
  Vector log_k_array = log(k_array);

  // Make splines for j_ell 
  generate_bessel_function_splines();

  //Line of sight integration to get Theta_ell(k)
  line_of_sight_integration(k_array);

  //Integration to get Cell 
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
   
  double eta_0 = cosmo->eta_of_x(0.0);
  Vector z_array;
  z_array = Utils::linspace(k_min*eta_0, k_max*eta_0, n_k);
  Vector j_ell_array(z_array.size());
  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    for(int j = 0; j < z_array.size(); j++){
           j_ell_array[j] = Utils::j_ell(ell, z_array[j]);
    }

    // Make the j_ell_splines[i] spline
    Spline ell_spline{"ell"};
    ell_spline.create(z_array, j_ell_array, "ell");
    j_ell_splines[i] = ell_spline;
    
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  double eta_0 = cosmo->eta_of_x(0.0);
  Vector x_array;
  const int npts = 500;
  x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts);

  for(size_t ik = 0; ik < k_array.size(); ik++){
     for(size_t ih = 0; ih < ells.size(); ih++){
  
        // The ODE for dTheta/dx
        ODEFunction dThetadx = [&](double x, const double *Theta, double *dThetadx){
        double z = k_array[ik]*(eta_0-cosmo->eta_of_x(x));
        dThetadx[0] = source_function(x, k_array[ik])*j_ell_splines[ih](z);
        return GSL_SUCCESS;
        };

       double Theta0_ini = 0.0; 
       Vector Theta_ic{Theta0_ini};
       ODESolver ode;
       ode.solve(dThetadx, x_array, Theta_ic);
       auto result_single = ode.get_data_by_component(0);
       // Store the result for a given ell and k
       result[ih][ik] = result_single[npts-1]; //for x=0
     }
  }
  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  
  //Make and store the splines for each ell
  for(size_t i = 0; i < nells; i++){
     Spline ThetaT_spline{"ThetaT"};
     ThetaT_spline.create(k_array,  thetaT_ell_of_k[i] , "ThetaT");
     thetaT_ell_of_k_spline[i] = ThetaT_spline;
  }  

}

//======================================
// Compute Cell (will be used for TT) 
//======================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
  Vector result(ells.size());
  for(size_t ih = 0; ih < ells.size(); ih++){
     // The ODE for dCell/dlogk
     ODEFunction dCelldlogk = [&](double log_k, const double *Cell, double *dCelldlogk){
     double k = exp(log_k);
     dCelldlogk[0] = 4.*M_PI*primordial_power_spectrum(k)*f_ell_spline[ih](k)*g_ell_spline[ih](k);
     return GSL_SUCCESS;
     };  
     double Cell_ini = 0.0; 
     Vector Cell_ic{Cell_ini};
     ODESolver ode2;
     ode2.solve(dCelldlogk, log_k_array, Cell_ic);
     auto result_single2 = ode2.get_data_by_component(0);
     result[ih] = result_single2[log_k_array.size()-1]; 
     } 
return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}


//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double k = k_mpc/Constants.Mpc;
  double P_primordial = 2.*M_PI*M_PI*primordial_power_spectrum(k)/(k*k*k);
  double c = Constants.c;
  double Phi = pert->get_Phi(x,k);
  double Omega_M_0 = cosmo->get_OmegaCDM(0.0) + cosmo->get_OmegaB(0.0);
  double a = exp(x);
  double H_0 = cosmo->get_H0();
  double Delta_M = (c*c*k*k*Phi*(2./3.)*a)/(Omega_M_0*H_0*H_0);
  double pofk = Delta_M*Delta_M*P_primordial;
  return pofk;
}
//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

//====================================================
// Output functions of k to file
//====================================================
void PowerSpectrum::output2(std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int k_number = 1000;
  auto kvalues = Utils::linspace(Constants.k_min*Constants.Mpc, Constants.k_max*Constants.Mpc, k_number);
  auto print_data = [&] (const double k_mpc) {
    
    fp << k_mpc/0.7                                 << " ";
    fp << get_matter_power_spectrum(0.0, k_mpc)*0.7*0.7*0.7/(Constants.Mpc*Constants.Mpc*Constants.Mpc) << " ";
    fp << thetaT_ell_of_k_spline[10](k_mpc/Constants.Mpc) << " ";   //ell=20
    fp << thetaT_ell_of_k_spline[32](k_mpc/Constants.Mpc) << " ";   //ell=500
    fp << thetaT_ell_of_k_spline[46](k_mpc/Constants.Mpc) << " ";   //ell=1200
    fp << thetaT_ell_of_k_spline[10](k_mpc/Constants.Mpc)*thetaT_ell_of_k_spline[10](k_mpc/Constants.Mpc)/(k_mpc/Constants.Mpc) << " ";   
    fp << thetaT_ell_of_k_spline[32](k_mpc/Constants.Mpc)*thetaT_ell_of_k_spline[32](k_mpc/Constants.Mpc)/(k_mpc/Constants.Mpc) << " ";  
    fp << thetaT_ell_of_k_spline[46](k_mpc/Constants.Mpc)*thetaT_ell_of_k_spline[46](k_mpc/Constants.Mpc)/(k_mpc/Constants.Mpc) << " ";    
    fp << "\n";
  };
  std::for_each(kvalues.begin(), kvalues.end(), print_data);
}






