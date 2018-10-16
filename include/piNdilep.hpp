#ifndef PINDILEP
#define PINDILEP

#include <map>
#include <string>
#include "Vectors.hpp"
using namespace Vectors;
#include "Spinors.hpp"
using namespace Spinors;
#include "units.hpp"
using namespace units_GeV;
#include "Array.hpp"
#include "Vrancx.hpp"

extern const double Gamma_ome0;

int geti(std::string);
double getd(std::string);

class DiffSigma_piN_Ndilep {
public:
  DiffSigma_piN_Ndilep(double srt, double M);
  Array<dcomplex,2> rho_prod(double costh) const;
  double B_coeff(double costh) const;
  double B_coeff_pipi(double costh) const;
  double dsig_dMdcosth_dcosthe_dphie(double costh, double costh_e, double phi_e) const;
  double dsig_dMdcosth(double costh) const;
  double dsig_dM() const;
  double dsig_pipi_dMdcosth(double costh) const;
  double dsig_pipi_dM() const;
private:
  const double srt; ///< srt(s)
  const double M; ///< dilepton mass
  double s;
  double M2;
  double mpi_pm;
  double mpi_pm2;
  double mpi_0;
  double mpi;
  double mpi2;
  double mrho;
  double mn;
  double mn2;
  double EN_in;
  double Epi_in;
  double pin_abs;
  FourVector p1;
  FourVector p2;
  FourVector q;
  FourVector P;
  double EN_out;
  double k0;
  double pout_abs;
  dcomplex F_rho;
  dcomplex F_ome;
};

#endif // PINDILEP
