#include "Vrancx.hpp"
#include "Config.hpp"

#include <cassert>

double O32(FourVector p, uint mu, uint nu, uint la) {
  return p(mu)*g_(nu,la) - p(nu)*g_(mu,la);
}


DiracMatrix O32(FourVector p, uint mu, uint la) {
  return p(mu)*gamma_(la) - g_(mu,la)*gamma_(p);
}

double O52(FourVector p, uint mu, uint nu, uint la, uint ro, uint si, uint ta) {
  return 1./4. * (
                  + O32(p,mu,la,si)*O32(p,nu,ro,ta) + O32(p,mu,ro,si)*O32(p,nu,la,ta)
                  + O32(p,mu,la,ta)*O32(p,nu,ro,si) + O32(p,mu,ro,ta)*O32(p,nu,la,si)
                  );
}


double O52(FourVector p, uint mu, uint nu, uint si, uint ta) {
  return 1./2. * (
                  + 2*p(mu)*p(nu)*g_(si,ta) + (p*p)*(g_(mu,si)*g_(nu,ta) + g_(mu,ta)*g_(nu,si))
                  - p(mu)*p(si)*g_(nu,ta) - p(nu)*p(ta)*g_(mu,si) - p(mu)*p(ta)*g_(nu,si) - p(nu)*p(si)*g_(mu,ta)
                  );
}

dcomplex vertexrhopipi(double g, uint mu, FourVector q1, FourVector q2) {
  return i_ * g * (q1(mu) - q2(mu));
}

dcomplex vertexsipipi(double g, FourVector q1, FourVector q2) {
  const double mpi = Config::get<double>("pi_pm.mass");
  return -i_ * g/(2*mpi) * (q1*q2);
}


DiracMatrix vertex1hNpi(double g, halfint spinParity, FourVector q) {
  const double mpi_pm = Config::get<double>("pi_pm.mass");
  const double mpi_0 = Config::get<double>("pi_0.mass");
  const double mpi = mpi_pm;   // we consider only the charged pion case
  DiracMatrix ret = -g/mpi*gamma_(q);
  return ((spinParity > 0) ? gamma5_*ret : ret);
}

DiracMatrix vertex1hNsi(double g, halfint spinParity) {
  const double mpi = Config::get<double>("pi_pm.mass");
  DiracMatrix ret = -i_*g*gamma_unit;
  return ((spinParity > 0) ? ret : ret*gamma5_);
}

DiracMatrix vertex1hNrho(double g1, halfint spinParity, FourVector pR, uint nu, FourVector k) {
  const double mrho = Config::get<double>("rho.mass");
  DiracMatrix ret = g1/mrho * sigma_(k,nu);
  return ((spinParity > 0) ? ret : ((pR(0) > 0) ? -gamma5_*ret : gamma5_*ret));
}

DiracMatrix vertex3hNpi(double g, halfint spinParity, uint muR, FourVector pR, FourVector q) {
  const double mpi_pm = Config::get<double>("pi_pm.mass");
  const double mpi_0 = Config::get<double>("pi_0.mass");
  const double mpi = mpi_pm;   // we consider only the charged pion case
  DiracMatrix ret = -i_*g/(mpi*mpi) * ((pR*q)*gamma_(muR) - q(muR)*gamma_(pR)); // ORIGINAL !!!!!!!!
  //DiracMatrix ret = g/(mpi*mpi) * ((pR*q)*gamma_(muR) - q(muR)*gamma_(pR));
  return ((spinParity > 0) ? ret : ret*gamma5_);
}

DiracMatrix vertex3hNsi(double g, halfint spinParity, uint muR, FourVector pR, FourVector q) {
  const double msi = Config::get<double>("sigma.mass");
  DiracMatrix ret = -i_*g/(msi*msi) * ((pR*q)*gamma_(muR) - q(muR)*gamma_(pR)); // ORIGINAL !!!!!!!!
  //DiracMatrix ret = g/(msi*msi) * ((pR*q)*gamma_(muR) - q(muR)*gamma_(pR));
  return ((spinParity > 0) ? ret*gamma5_ : ret);
}

DiracMatrix vertex3h3hpi(double g, 
			 halfint spinParityR, uint muR, FourVector pR, 
			 halfint spinParityD, uint muD, FourVector pD, 
			 FourVector q) {
  const double mpi = Config::get<double>("pi_pm.mass"); // we consider only the charged pion case
  DiracMatrix Gamma = ((spinParityR == spinParityD) ? gamma5_ : gamma_unit);
  DiracMatrix ret = gamma_null;
  for (uint al(0); al<4; al++) {
    if (pR(0)>0) { // incoming resonance
      ret += g/POW<3>(mpi) * O32(pD,al,muD) * Gamma * gamma_(q) * O32(pR,al,muR) * sign_(al);
    } else {
      ret += g/POW<3>(mpi) * O32(pR,al,muR) * Gamma * gamma_(q) * O32(pD,al,muD) * sign_(al);
    }
  }
  return ret;
}

DiracMatrix vertex3hNrho(double g1, double g2, double g3, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k) {
  const double mrho = Config::get<double>("rho.mass");
  FourVector pN(-pR-k);
  DiracMatrix ret1(gamma_null);
  DiracMatrix ret2(gamma_null);
  DiracMatrix ret3(gamma_null);
  for (uint ro : {0,1,2,3}) {
    if (g1 != 0) {
      if (pR(0) > 0) { // incoming resonance - only valid if antibaryons are not present!
        ret1 += O32(k,ro,nu)*O32(pR,ro,muR) * sign_(ro);
      } else { // outgiong resonance
        ret1 += O32(pR,ro,muR)*O32(k,ro,nu) * sign_(ro);
      }
    }
    if (g2 != 0) {
      ret2 += O32(pR,ro,muR)*((pN*k)*g_(ro,nu) - k(ro)*pN(nu)) * sign_(ro);
    }
    if (g3 != 0) {
      ret3 += O32(pR,ro,muR)*((k*k)*g_(ro,nu) - k(ro)*k(nu)) * sign_(ro);
    }
  }
  if (spinParity > 0) {
    ret1 = ret1 * ( (pR(0) > 0) ? gamma5_ : -gamma5_ );
    ret2 = ret2 * gamma5_;
    ret3 = ret3 * gamma5_;
  }
  return i_*g1/(4.*mrho*mrho)*ret1 - 1./(8.*mrho*mrho*mrho) * (g2*ret2 + g3*ret3);
}

DiracMatrix vertex5hNpi(double g, halfint spinParity, uint mu, uint nu, FourVector pR, FourVector q) {
  const double mpi_pm = Config::get<double>("pi_pm.mass");
  const double mpi_0 = Config::get<double>("pi_0.mass");
  const double mpi = mpi_pm;   // we consider only the charged pion case
  DiracMatrix ret = g/POW<4>(mpi) * (
                                     + g_(mu,nu)*POW<2>(pR*q) + q(mu)*q(nu)*(pR*pR)
                                     - (pR(mu)*q(nu) + q(mu)*pR(nu))*(pR*q)
                                     ) * gamma_unit;
  if (spinParity > 0) { ret = ret * ( (pR(0) > 0) ? -gamma5_ : gamma5_ ); }
  return ret;
}

DiracMatrix vertex5hNsi(double g, halfint spinParity, uint mu, uint nu, FourVector pR, FourVector q) {
  const double msi = Config::get<double>("sigma.mass");
  DiracMatrix ret = g/POW<4>(msi) * (
                                     + g_(mu,nu)*POW<2>(pR*q) + q(mu)*q(nu)*(pR*pR)
                                     - (pR(mu)*q(nu) + q(mu)*pR(nu))*(pR*q)
                                     ) * gamma_unit;
  if (spinParity < 0) { ret = ret * ( (pR(0) > 0) ? -gamma5_ : gamma5_ ); }
  return ret;
}

DiracMatrix vertex5hNrho(double g1, double g2, double g3, halfint spinParity, uint si, uint ta, FourVector pR, uint la, FourVector k) {
  const double mrho = Config::get<double>("rho.mass");
  const double mn = Config::get<double>("Nucleon.mass");
  FourVector pN(-pR-k);
  DiracMatrix ret1(gamma_null);
  double scal2(0);
  double scal3(0);
  for (uint mu : {0,1,2,3}) {
    for (uint nu : {0,1,2,3}) {
      if (g1!=0) { ret1 += O32(k,nu,la) * pN(mu) * O52(pR,mu,nu,si,ta) * sign_(mu)*sign_(nu); }
      if (g2!=0) { scal2 += (pN(la)*k(nu) - g_(la,nu)*(pN*k)) * pN(mu) * O52(pR,mu,nu,si,ta) * sign_(mu)*sign_(nu); }
      if (g3!=0) { scal3 += (k(la)*k(nu) - g_(la,nu)*(k*k)) * pN(mu) * O52(pR,mu,nu,si,ta) * sign_(mu)*sign_(nu); }
    }
  }
  DiracMatrix ret2;
  DiracMatrix ret3;
  if (spinParity < 0) {
    ret1 = ret1 * gamma5_;
    ret2 = gamma5_ * scal2;
    ret3 = gamma5_ * scal3;
  } else {
    ret2 = gamma_unit * scal2;
    ret3 = gamma_unit * scal3;
    if (pR(0) > 0) {
      ret2 *= -1.;
      ret3 *= -1.;
    }
  }
  return -i_/POW<4>(2.*mrho) * (
                              + g1 * ret1
                              + g2/(2.*mrho) * ret2
                              + g3/(2.*mrho) * ret3
                              );
} 

/**
   Projectors used in propagators. Terms proportional to the momentum are dropped,
   because the vertices give zero when contracted with the momentum.
*/
DiracMatrix P3h(FourVector p, double m, uint mu, uint nu) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  return - g_(mu,nu)*gamma_unit + gamma_(mu)*gamma_(nu)/3.;
}


DiracMatrix P5h(FourVector p, double m, uint mu, uint nu, uint la, uint ro) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  assert(isIndex4(la));
  assert(isIndex4(ro));
  return 
    ( - (g_(mu,la)*g_(nu,ro) + g_(mu,ro)*g_(nu,la))/2. + g_(mu,nu)*g_(la,ro)/5. ) * gamma_unit
    + ( g_(mu,la)*gamma_(nu)*gamma_(ro) + g_(mu,ro)*gamma_(nu)*gamma_(la)
        + g_(nu,la)*gamma_(mu)*gamma_(ro) + g_(nu,ro)*gamma_(mu)*gamma_(la) )/10.;
}


/**
   Spin 3/2 projector.
*/
DiracMatrix pro3half(FourVector p, double m, uint mu, uint nu) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  double m2 = m*m;
  return (- g_(mu,nu) + 2./(3.*m2)*p(mu)*p(nu))*gamma_unit + gamma_(mu)*gamma_(nu)/3.
    + (gamma_(mu)*p(nu) - gamma_(nu)*p(mu))/(3.*m);
}

double _G(FourVector p, double m, uint mu, uint nu) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  return - g_(mu,nu) + p(mu)*p(nu)/(m*m);
}

DiracMatrix _T(FourVector p, double m, uint mu, uint nu) {
  const double m2(m*m); 
  return -1./2. * (gamma_(mu)*gamma_(nu) - gamma_(nu)*gamma_(mu))
    + 1./(2.*m2) * p(mu)*(gamma_(p)*gamma_(nu) - gamma_(nu)*gamma_(p))
    - 1./(2.*m2) * p(nu)*(gamma_(p)*gamma_(mu) - gamma_(mu)*gamma_(p));
}

/**
   Spin 5/2 projector.
*/
DiracMatrix pro5half(FourVector p, double m, uint mu1, uint mu2, uint nu1, uint nu2) {
  assert(isIndex4(mu1));
  assert(isIndex4(mu2));
  assert(isIndex4(nu1));
  assert(isIndex4(nu2));
  double m2 = p*p;
  return - 3./10. * (_G(p,m,mu1,nu1)*_G(p,m,mu2,nu2) + _G(p,m,mu1,nu2)*_G(p,m,mu2,nu1)) * gamma_unit
    + 1./5. * _G(p,m,mu1,mu2)*_G(p,m,nu1,nu2) * gamma_unit
    + 1./10. * (_T(p,m,mu1,nu1)*_G(p,m,mu2,nu2) + _T(p,m,mu2,nu2)*_G(p,m,mu1,nu1) + _T(p,m,mu1,nu2)*_G(p,m,mu2,nu1) + _T(p,m,mu2,nu1)*_G(p,m,mu1,nu2));
}

