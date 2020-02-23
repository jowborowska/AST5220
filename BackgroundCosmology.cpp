#include "BackgroundCosmology.h"

# define M_PI           3.14159265358979323846  /* pi */

    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaLambda,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaLambda(OmegaLambda),
  Neff(Neff), 
  TCMB(TCMB)
{
 H0 = Constants.H0_over_h*h; //Hubble parameter today (1/s)
 //std:: cout << H0 << "H0 parameter";
 
 double rho_crit; // Critical density today
 rho_crit = 3.0*H0*H0/(8.0*M_PI*Constants.G);
 OmegaR = M_PI*M_PI*pow(Constants.k_b*TCMB,4)/(rho_crit*15.0*pow(Constants.hbar,3)*pow(Constants.c,5)); 
 OmegaNu = 0;
 OmegaK = 0;
 //OmegaK = 1.0 - OmegaCDM - OmegaR - OmegaNu - OmegaLambda;
}

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
 
  Vector x_array;
  const int npts = 500;
  x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
  detadx[0] = Constants.c/Hp_of_x(x);
  return GSL_SUCCESS;
  };

  double eta0_ini = 0.0; //approximate eta for initial x=ln(10e-7)
  Vector eta_ic{eta0_ini};
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);
  auto result = ode.get_data();

  Vector y_array(npts); //gather result into y_array
  for(int i = 0; i < npts; i++){
      y_array[i] = result[i][0];
  }
  
  eta_of_x_spline.create(x_array, y_array, "eta"); //create spline

  Utils::EndTiming("Eta");
}


double BackgroundCosmology::H_of_x(double x) const{
  double a = exp(x);
  double H = get_H0()*sqrt((OmegaB+OmegaCDM)*pow(a,-3) + (OmegaR+OmegaNu)*pow(a,-4) + OmegaLambda);
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  double a = exp(x);
  double Hx = H_of_x(x);
  double Hp = a*Hx;
  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  double a = exp(x);
  double under_root = (OmegaB+OmegaCDM)*pow(a,-1) + (OmegaR+OmegaNu)*pow(a,-2) + OmegaLambda*pow(a, 2);
  double dunder_root_dx = -(OmegaB+OmegaCDM)*pow(a,-1) -2*(OmegaR+OmegaNu)*pow(a,-2) + 2*OmegaLambda*pow(a, 2);
  double dHpdx = H0*0.5*pow(under_root, -0.5)*dunder_root_dx;
  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double a = exp(x);
  double under_root = (OmegaB+OmegaCDM)*pow(a,-1) + (OmegaR+OmegaNu)*pow(a,-2) + OmegaLambda*pow(a, 2);
  double dunder_root_dx = -(OmegaB+OmegaCDM)*pow(a,-1) -2*(OmegaR+OmegaNu)*pow(a,-2) + 2*OmegaLambda*pow(a, 2);
  double ddunder_root_dx2 = (OmegaB+OmegaCDM)*pow(a,-1) +4*(OmegaR+OmegaNu)*pow(a,-2) + 4*OmegaLambda*pow(a, 2);
  double ddHpdx = H0*0.5*(-0.5)*pow(under_root, -1.5)*pow(dunder_root_dx,2) + H0*0.5*pow(under_root, -0.5)*ddunder_root_dx2;
  return ddHpdx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  double a = exp(x);
  double OmegaB_now = H0*H0*OmegaB*pow(a,-3)*pow(H_of_x(x),-2);
  return OmegaB_now;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  double a = exp(x);
  double OmegaR_now = H0*H0*OmegaR*pow(a,-4)*pow(H_of_x(x),-2);
  return OmegaR_now;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;
  return 0.0; //always 0, since we ignore neutrinos
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  double a = exp(x);
  double OmegaCDM_now = H0*H0*OmegaCDM*pow(a,-3)*pow(H_of_x(x),-2);
  return OmegaCDM_now;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  double a = exp(x);
  double OmegaLambda_now = H0*H0*OmegaLambda*pow(H_of_x(x),-2);
  return OmegaLambda_now; 
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;
   double a = exp(x);
   double OmegaK_now = H0*H0*OmegaK*pow(a,-2)*pow(H_of_x(x),-2);
   return OmegaK_now; //this will always be zero, since we set OmegaK = 0
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB() const{ 
  return TCMB; 
}


//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = Constants.x_start;
  const double x_max = Constants.x_end;
  const int    n_pts =  10000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

