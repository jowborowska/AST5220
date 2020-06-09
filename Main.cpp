#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");
  
  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.7;
  double OmegaB      = 0.046;
  double OmegaCDM    = 0.224;
  double OmegaLambda = 0.72995;
  double Neff        = 0; //3.046 or 0 if ignoring neutrinos
  double TCMB        = 2.725;

  // Recombination parameters
  double Yp          = 0;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaLambda, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("cosmology.txt");

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(h, OmegaB, TCMB, &cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.15/Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.15.txt");
  
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert);
  power.solve();
  
  // Output CMB power spectrum and k-dependent functions
  power.output("cells.txt");
  power.output2("ks.txt");
  

  Utils::EndTiming("Everything");
}
