#ifndef _PERTURBATIONS_HEADER
#define _PERTURBATIONS_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <vector>
#include <fstream>
#include <algorithm>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class Perturbations{
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
   
    // The scales we integrate over
    const int n_k        = 100;
    const double k_min   = Constants.k_min;
    const double k_max   = Constants.k_max;
    
    // Start and end of the time-integration
    const int n_x        = 1000;
    const double x_start = Constants.x_start;
    const double x_end   = Constants.x_end;

   
   // Splines of scalar perturbations quantities
    Spline2D delta_cdm_spline{"delta_cdm_spline"};
    Spline2D delta_b_spline{"delta_b_spline"};
    Spline2D v_cdm_spline{"v_cdm_spline"};
    Spline2D v_b_spline{"v_b_spline"};
    Spline2D Phi_spline{"Phi_spline"};
    Spline2D Theta0_spline{"Theta0_spline"};
    Spline2D Theta1_spline{"Theta1_spline"};
    Spline2D Theta2_spline{"Theta2_spline"};
    Spline2D Theta3_spline{"Theta3_spline"};
    Spline2D Theta4_spline{"Theta4_spline"};
    Spline2D Theta5_spline{"Theta5_spline"};
    Spline2D Theta6_spline{"Theta6_spline"};
    Spline2D Theta7_spline{"Theta7_spline"};
    Spline2D Psi_spline{"Psi_spline"};
   
    // Splines of source function (ST for temperature)
    Spline2D ST_spline{"ST_spline"};
    
    
    //==========================================================
    // [1] Tight coupling ODE system
    //==========================================================

    // Set the initial conditions at the start (which is in tight coupling)
    Vector set_ic(
        const double x, 
        const double k) const;
    
    // Right hand side of the ODE in the tight coupling regime
    int rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx);
    
    // Compute the time when tight coupling ends
    double get_tight_coupling_time(const double k) const;
    
    //==========================================================
    // [2] The full ODE system 
    //==========================================================
    
    // Set initial condition after tight coupling
    Vector set_ic_after_tight_coupling(
        const Vector &y_tight_coupling, 
        const double x, 
        const double k) const;

    // Right hand side of the ODE in the full regime
    int rhs_full_ode(double x, double k, const double *y, double *dydx);
    
    //==========================================================
    // [3] Integrate the full system
    //==========================================================
    
    // Integrate perturbations and spline the result
    void integrate_perturbations();
    
    //==========================================================
    // [4] Compute source functions from the result
    //==========================================================
    
    // Compute source functions and spline the result
    void compute_source_functions();

  public:

    // Constructors
    Perturbations() = default;
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec); 

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output info to file
    void output(const double k, const std::string filename) const;

    // Get the quantities we have integrated
    double get_delta_cdm(const double x,  const double k) const;
    double get_delta_b(const double x, const double k) const;
    double get_v_cdm(const double x, const double k) const;
    double get_v_b(const double x, const double k) const;
    double get_Phi(const double x, const double k) const;
    double get_Psi(const double x, const double k) const;
    double get_Theta(const double x, const double k, const int ell) const;
    double get_Source_T(const double x, const double k) const;
};

#endif
