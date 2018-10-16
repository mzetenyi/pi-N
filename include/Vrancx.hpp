#ifndef VRANCX_HPP
#define VRANCX_HPP

#include "Vectors.hpp"
using namespace Vectors;
#include "Spinors.hpp"
using namespace Spinors;
#include "units.hpp"
using namespace units_GeV;

double O32(FourVector p, uint mu, uint nu, uint la);
DiracMatrix O32(FourVector p, uint mu, uint la);

double O52(FourVector p, uint mu, uint nu, uint la, uint ro, uint si, uint ta);
double O52(FourVector p, uint mu, uint nu, uint la);

dcomplex vertexrhopipi(double g, uint nu, FourVector q1, FourVector q2);

dcomplex vertexsipipi(double g, FourVector q1, FourVector q2);

DiracMatrix vertex1hNpi(double g, halfint spinParity, FourVector q);

DiracMatrix vertex1hNsi(double g, halfint spinParity);

DiracMatrix vertex1hNrho(double g1, halfint spinParity, FourVector pR, uint nu, FourVector k);

DiracMatrix vertex3hNpi(double g, halfint spinParity, uint muR, FourVector pR, FourVector q);

DiracMatrix vertex3hNsi(double g, halfint spinParity, uint muR, FourVector pR, FourVector q);

DiracMatrix vertex3hNrho(double g1, double g2, double g3, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);

DiracMatrix vertex5hNpi(double g, halfint spinParity, uint mu, uint nu, FourVector pR, FourVector q);

DiracMatrix vertex5hNsi(double g, halfint spinParity, uint mu, uint nu, FourVector pR, FourVector q);

DiracMatrix vertex5hNrho(double g1, double g2, double g3, halfint spinParity, uint si, uint ta, FourVector pR, uint la, FourVector k);

DiracMatrix vertex3h3hpi(double g, 
			 halfint spinParityR, uint muR, FourVector pR, 
			 halfint spinParityD, uint muD, FourVector pD, 
			 FourVector q);

/** spin-3/2 propagator, without the terms that give 0 contrib. in the case of Vrancx Lagrangians. */
DiracMatrix P3h(FourVector p, double m, uint mu, uint nu);

/** spin-5/2 propagator, without the terms that give 0 contrib. in the case of Vrancx Lagrangians. */
DiracMatrix P5h(FourVector p, double m, uint mu, uint nu, uint la, uint ro);

DiracMatrix pro3half(FourVector p, double m, uint mu, uint nu);

double _G(FourVector p, double m, uint mu, uint nu);

DiracMatrix _T(FourVector p, double m, uint mu, uint nu);

DiracMatrix pro5half(FourVector p, double m, uint mu1, uint mu2, uint nu1, uint nu2);


#endif // VRANCX_HPP
