#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "piNdilep.hpp"
#include "Array.hpp"
#include "wavefunc.hpp"
#include "RandomNumberGenerator.hpp"
#include "Config.hpp"
#include "Histogram.hpp"

using namespace std;

/**
   Calculate diff. cross section of the process pi + N -> N + e+ + e-

*/

const double alpha(1/137.036);
const double e(sqrt(4.*pi_*alpha));

const double grho(4.96);

bool sch(true);
bool uch(true);
bool D1232(false);
bool N1440(false);
bool N1535(false);
bool N1520(false);
bool D1600(false);
bool D1620(false);
bool N1650(false);
bool N1675(false);
bool N1680(false);

//*
double Gamma_rho(double m) {
  //static const double delta(0.3*GeV);
  static const double mpi_pm = Config::get<double>("pi_pm.mass");
  static const double mpi_0 = Config::get<double>("pi_0.mass");
  //static const double mpi = (2.*mpi_pm + mpi_0)/3.;
  static const double mpi = mpi_pm;
  static const double mrho = Config::get<double>("rho.mass");
  const double Gamma0 = Config::get<double>("rho.width");
  double x = m*m - 4.*mpi*mpi;
  if (x<=0) { return 0; }
  double q = m*sqrt(x);  // pion momentum
  double q0 = mrho*sqrt(mrho*mrho - 4.*mpi*mpi); // pion momentum at the pole
  //return Gamma0 * mrho/m * POW<3>(q/q0) * (q0*q0 + delta*delta)/(q*q + delta*delta);
  return Gamma0 * pow(q/q0,3);
}
//*/

double Gamma_R(double s, double m0, double G0, int l) {
  static const double mpi_pm = Config::get<double>("pi_pm.mass");
  static const double mpi_0 = Config::get<double>("pi_0.mass");
  //static const double mpi = (2.*mpi_pm + mpi_0)/3.;
  static const double mpi = mpi_pm;
  static const double mn = Config::get<double>("Nucleon.mass");
  if ((s>0) and (sqrt(s)>mn+mpi)) {
    double q  = sqrt(lambda(s,    mn*mn,mpi*mpi))/(2.*sqrt(s));
    double q0 = sqrt(lambda(m0*m0,mn*mn,mpi*mpi))/(2.*m0);
    double d2 = POW<2>(m0-mn-mpi) + G0*G0/4.;
    return G0 * m0/sqrt(s) * pow(q/q0,2*l+1) * pow((q0*q0+d2)/(q*q+d2),l+1);
  }
  return 0;
}

dcomplex BW(double s, double M, double Gamma) {
  if (s>0) {
    return 1./(s - M*M + i_*sqrt(s)*Gamma);
    //return 1./(s - M*M + i_*M*Gamma);
  }
  return 1./(s - M*M);
}

DiffSigma_piN_Ndilep::DiffSigma_piN_Ndilep(double srt, double M) : 
  srt(srt),
  M(M) {
  s = srt*srt;
  M2 = M*M;
  mpi_pm = Config::get<double>("pi_pm.mass");
  mpi_pm2 = mpi_pm*mpi_pm;
  mpi_0 = Config::get<double>("pi_0.mass");
  mpi = mpi_pm;
  //mpi = (2.*mpi_pm + mpi_0)/3.;
  mpi2 = mpi*mpi;
  mrho = Config::get<double>("rho.mass");
  mn = Config::get<double>("Nucleon.mass");
  mn2 = mn*mn;
  if (srt<mn+mpi) { cerr << "sqrt(s) is below threshold of pi+N" << endl; exit(0); }
  if (srt<mn+M) { cerr << "sqrt(s) is below threshold of N + dilep of mass M" << endl; exit(0); }
  EN_in = (s-mpi2+mn2)/(2.*srt);
  Epi_in = (s+mpi2-mn2)/(2.*srt);
  pin_abs = sqrt(lambda(s,mpi2,mn2))/(2.*srt);
  p1 = FourVector(EN_in,0,0,-pin_abs);
  q = FourVector(Epi_in,0,0,pin_abs);
  P = p1+p2;
  EN_out = (s-M2+mn2)/(2.*srt);
  k0 = (s+M2-mn2)/(2.*srt);
  pout_abs = sqrt(lambda(s,M2,mn2))/(2.*srt);
  F_rho = -e/grho * (M2)/(M2 - mrho*mrho + i_*sqrt(M2)*Gamma_rho(M));
  //F_rho = -e/grho * (M2)/(M2 - mrho*mrho + i_*mrho*Gamma_rho(M));
  }


//*

/**
   Calculate the production density matrix
*/
Array<dcomplex,2> DiffSigma_piN_Ndilep::rho_prod(double costh) const {
  // isospin factors:
  const double isoNs = -sqrt(2);
  const double isoNu =  sqrt(2);
  const double isoDs = sqrt(2)/3;
  const double isoDu = sqrt(2)/3;

  double sinth = sqrt(1-costh*costh);
  ThreeVector kk = pout_abs*ThreeVector(sinth,0,costh);
  FourVector p2 = FourVector(EN_out,-kk);
  FourVector k = FourVector(k0,kk);

  //  PR(p1+q-p2-k);

  double ch = k0/M;
  double sh = -pout_abs/M;
  FourTensor Boost(ch, 0,  0,  sh,
                   0, -1,  0,  0,
                   0,  0, -1,  0,
                  -sh, 0,  0, -ch );

  FourTensor Rot( 1,     0,   0,    0,
                  0,  -costh, 0,  -sinth,
                  0,    0,   -1,    0,
                  0,   sinth, 0,  -costh );

  //  PR(costh); PR(sinth); PR(costh*costh + sinth*sinth);
  // FourVector ktilde(M,0,0,0);
  // FourVector kprime = Boost*ktilde;
  // FourVector knew = Rot*kprime;

  Array<FourVector,1> epsRe(idx_s1);
  Array<FourVector,1> epsIm(idx_s1);

  // photon wave vectors:
  FourVector epsRe0(0,0,0,1.);
  FourVector epsIm0(0,0,0,0);
  epsRe(_0) = Rot*Boost*epsRe0;
  epsIm(_0) = Rot*Boost*epsIm0;

  FourVector epsRep(0,-1./sqrt(2.),0,0);
  FourVector epsImp(0,0,-1./sqrt(2.),0);
  epsRe(_1) = Rot*Boost*epsRep;
  epsIm(_1) = Rot*Boost*epsImp;
  //  PR(epsRep); PR(epsImp); PR(epsRe(1)); PR(epsIm(1));

  FourVector epsRem(0,1./sqrt(2.),0,0);
  FourVector epsImm(0,0,-1./sqrt(2.),0);
  epsRe(-_1) = Rot*Boost*epsRem;
  epsIm(-_1) = Rot*Boost*epsImm;
  //  PR(epsRem); PR(epsImm); PR(epsRe(-1)); PR(epsIm(-1));

  double uu = (p1-k)*(p1-k);

  //  if (uu>0) PR(uu);

  //  cerr << "BWratio = " << (2.*real(BWs*conj(BWu))) / (BWs*conj(BWs)) << endl;
  //PR(BW(uu,mr,Gr)*conj(BW(uu,mr,Gr)));

  /*
    PR(p1*q) PR(p1*p2) PR(q*k);

    PR(p1) PR(p2) PR(q) PR(k);
  */
  Array<dcomplex,3> Mhad(idx_s1h,idx_s1h,idx_lor);
  
  ubar_ ubar(half,p2);
  u_ u(half,p1);

  for (halfint la1 : {-half,half}) {
    for (halfint la2 : {-half,half}) {
      for (halfint mu(0); mu<4; mu++) {
        Mhad(la1,la2,mu) = 0;
        //------------  spin 1/2+  ------------------------------
        if (N1440) {
          const double mr = Config::get<double>("N1440.mass");
          const double Gr = Config::get<double>("N1440.width");
          const double g0 = Config::get<double>("N1440.g0");
          const double g1 = Config::get<double>("N1440.g1");
          const double l = Config::get<int>("N1440.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoNs;
          const double iso_u = isoNu;
          double fac(1.);
          if (Config::exists("negN1440")) fac = -1.;
          
          Mhad(la1,la2,mu) += 
            fac * F_rho *
            ubar(0,la2) *
            ( iso_s *
              ( (not sch) ? gamma_null
                : ( vertex1hNrho(g1,half,p1+q,mu,-k) *
                    i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                    vertex1hNpi(g0,half,q) ) )
              + iso_u *
              ( (not uch) ? gamma_null
                : ( vertex1hNpi(g0,half,q) *
                    i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                    vertex1hNrho(g1,half,-p1+k,mu,-k) )
                )
              ) * u(0,la1);
        }
        //------------  spin 1/2-  ------------------------------
        if (N1535) {
          const double mr = Config::get<double>("N1535.mass");
          const double Gr = Config::get<double>("N1535.width");
          const double g0 = Config::get<double>("N1535.g0");
          const double g1 = Config::get<double>("N1535.g1");
          const double l = Config::get<int>("N1535.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoNs;
          const double iso_u = isoNu;
          double fac(1.);
          if (Config::exists("negN1535")) fac = -1.;

          Mhad(la1,la2,mu) += 
            fac * F_rho *
            ubar(0,la2) *
            ( iso_s *
              ( (not sch) ? gamma_null
                : ( vertex1hNrho(g1,-half,p1+q,mu,-k) *
                    i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                    vertex1hNpi(g0,-half,q) ) )
              + iso_u *
              ( (not uch) ? gamma_null
                : ( vertex1hNpi(g0,-half,q) *
                    i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                    vertex1hNrho(g1,-half,-p1+k,mu,-k) )
                )
              ) * u(0,la1);
        } 
        if (N1650) {
          const double mr = Config::get<double>("N1650.mass");
          const double Gr = Config::get<double>("N1650.width");
          const double g0 = Config::get<double>("N1650.g0");
          const double g1 = Config::get<double>("N1650.g1");
          const double l = Config::get<int>("N1650.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoNs;
          const double iso_u = isoNu;
          double fac(1.);
          if (Config::exists("negN1650")) fac = -1.;

          Mhad(la1,la2,mu) += 
            fac * F_rho *
            ubar(0,la2) *
            ( iso_s *
              ( (not sch) ? gamma_null
                : ( vertex1hNrho(g1,-half,p1+q,mu,-k) *
                    i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                    vertex1hNpi(g0,-half,q) ) )
              + iso_u *
              ( (not uch) ? gamma_null
                : ( vertex1hNpi(g0,-half,q) *
                    i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                    vertex1hNrho(g1,-half,-p1+k,mu,-k) )
                )
              ) * u(0,la1);
        } 
        if (D1620) {
          const double mr = Config::get<double>("D1620.mass");
          const double Gr = Config::get<double>("D1620.width");
          const double g0 = Config::get<double>("D1620.g0");
          const double g1 = Config::get<double>("D1620.g1");
          const double l = Config::get<int>("D1620.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoDs;
          const double iso_u = isoDu;
          double fac(1.);
          if (Config::exists("negD1620")) fac = -1.;

          Mhad(la1,la2,mu) += 
            fac * F_rho *
            ubar(0,la2) *
            ( iso_s *
              ( (not sch) ? gamma_null
                : ( vertex1hNrho(g1,-half,p1+q,mu,-k) *
                    i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                    vertex1hNpi(g0,-half,q) ) )
              + iso_u *
              ( (not uch) ? gamma_null
                : ( vertex1hNpi(g0,-half,q) *
                    i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                    vertex1hNrho(g1,-half,-p1+k,mu,-k) )
                )
              ) * u(0,la1);
        } 
        //------------  spin 3/2+  ------------------------------
        if (D1232) {
          const double mr = Config::get<double>("D1232.mass");
          const double Gr = Config::get<double>("D1232.width");
          const double g0 = Config::get<double>("D1232.g0");
          const double g1 = Config::get<double>("D1232.g1");
          const double g2 = Config::get<double>("D1232.g2");
          const double g3 = Config::get<double>("D1232.g3");
          const double l = Config::get<int>("D1232.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoDs;
          const double iso_u = isoDu;
          double fac(1.);
          if (Config::exists("negD1232")) fac = -1.;

          for (uint al(0); al<4; al++) {
            for (uint be(0); be<4; be++) {
              Mhad(la1,la2,mu) += 
                fac * F_rho *
                ubar(0,la2) *
                ( iso_s *
                  ( (not sch) ? gamma_null
                    : ( vertex3hNrho(g1,g2,g3,3*half,al,p1+q,mu,-k) *
                        i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                        P3h(p1+q,mr,al,be)*sign_(al)*sign_(be) *
                        vertex3hNpi(g0,3*half,be,-p1-q,q) ) )
                    + iso_u *
                  ( (not uch) ? gamma_null
                    : ( vertex3hNpi(g0,3*half,be,p1-k,q) *
                        i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                        P3h(p1-k,mr,be,al)*sign_(al)*sign_(be) *
                        vertex3hNrho(g1,g2,g3,3*half,al,-p1+k,mu,-k) )
                    )
                  ) * u(0,la1);
            }
          }
        }
        if (D1600) {
          const double mr = Config::get<double>("D1600.mass");
          const double Gr = Config::get<double>("D1600.width");
          const double g0 = Config::get<double>("D1600.g0");
          const double g1 = Config::get<double>("D1600.g1");
          const double g2 = Config::get<double>("D1600.g2");
          const double g3 = Config::get<double>("D1600.g3");
          const double l = Config::get<int>("D1600.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoDs;
          const double iso_u = isoDu;
          double fac(1.);
          if (Config::exists("negD1600")) fac = -1.;

          for (uint al(0); al<4; al++) {
            for (uint be(0); be<4; be++) {
              Mhad(la1,la2,mu) += 
                fac * F_rho *
                ubar(0,la2) *
                ( iso_s *
                  ( (not sch) ? gamma_null
                    : ( vertex3hNrho(g1,g2,g3,3*half,al,p1+q,mu,-k) *
                        i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                        P3h(p1+q,mr,al,be)*sign_(al)*sign_(be) *
                        vertex3hNpi(g0,3*half,be,-p1-q,q) ) )
                    + iso_u *
                  ( (not uch) ? gamma_null
                    : ( vertex3hNpi(g0,3*half,be,p1-k,q) *
                        i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                        P3h(p1-k,mr,be,al)*sign_(al)*sign_(be) *
                        vertex3hNrho(g1,g2,g3,3*half,al,-p1+k,mu,-k) ) )
                  )
                * u(0,la1);
            }
          }
        }
        //------------  spin 3/2-  ------------------------------
        if (N1520) {
          const double mr = Config::get<double>("N1520.mass");
          const double Gr = Config::get<double>("N1520.width");
          const double g0 = Config::get<double>("N1520.g0");
          const double g1 = Config::get<double>("N1520.g1");
          const double g2 = Config::get<double>("N1520.g2");
          const double g3 = Config::get<double>("N1520.g3");
          const double l = Config::get<int>("N1520.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoNs;
          const double iso_u = isoNu;
          double fac(1.);
          if (Config::exists("negN1520")) fac = -1.;

          for (uint al(0); al<4; al++) {
            for (uint be(0); be<4; be++) {
              Mhad(la1,la2,mu) +=
                fac * F_rho * 
                ubar(0,la2) *
                ( iso_s *
                  ( (not sch) ? gamma_null
                    : ( vertex3hNrho(g1,g2,g3,-3*half,al,p1+q,mu,-k) *
                        i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                        P3h(p1+q,mr,al,be)*sign_(al)*sign_(be) *
                        vertex3hNpi(g0,-3*half,be,-p1-q,q) ) )
                  + iso_u *
                  ( (not uch) ? gamma_null
                    : ( vertex3hNpi(g0,-3*half,be,p1-k,q) *
                        i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                        P3h(p1-k,mr,be,al)*sign_(al)*sign_(be) *
                        vertex3hNrho(g1,g2,g3,-3*half,al,-p1+k,mu,-k) ) )
                  )
                * u(0,la1);
            }
          }
        }
        //------------  spin 5/2+  ------------------------------
        if (N1680) {
          const double mr = Config::get<double>("N1680.mass");
          const double Gr = Config::get<double>("N1680.width");
          const double g0 = Config::get<double>("N1680.g0");
          const double g1 = Config::get<double>("N1680.g1");
          const double g2 = Config::get<double>("N1680.g2");
          const double g3 = Config::get<double>("N1680.g3");
          const double l = Config::get<int>("N1680.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoNs;
          const double iso_u = isoNu;
          double fac(1.);
          if (Config::exists("negN1680")) fac = -1.;

          for (uint al1(0); al1<4; al1++) {
            for (uint al2(0); al2<4; al2++) {
              for (uint be1(0); be1<4; be1++) {
                for (uint be2(0); be2<4; be2++) {
                  Mhad(la1,la2,mu) += 
                    fac * F_rho *
                    ubar(0,la2) * 
                    ( iso_s *
                      ( (not sch) ? gamma_null
                        : ( vertex5hNrho(g1,g2,g3,5*half,al1,al2,p1+q,mu,-k) *
                            i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                            P5h(p1+q,mr,al1,al2,be1,be2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) * 
                            vertex5hNpi(g0,5*half,be1,be2,-p1-q,q) ) )
                      + iso_u *
                      ( (not uch) ? gamma_null
                        : ( vertex5hNpi(g0,5*half,be1,be2,p1-k,q) *
                            i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                            P5h(p1-k,mr,be1,be2,al1,al2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) * 
                            vertex5hNrho(g1,g2,g3,5*half,al1,al2,-p1+k,mu,-k) ) )
                      )
                    * u(0,la1);
                }
              }
            }
          }
        }
        //------------  spin 5/2-  ------------------------------
        if (N1675) {
          const double mr = Config::get<double>("N1675.mass");
          const double Gr = Config::get<double>("N1675.width");
          const double g0 = Config::get<double>("N1675.g0");
          const double g1 = Config::get<double>("N1675.g1");
          const double g2 = Config::get<double>("N1675.g2");
          const double g3 = Config::get<double>("N1675.g3");
          const double l = Config::get<int>("N1675.l");
          const dcomplex BWs = BW(s,mr,Gamma_R(s,mr,Gr,l));
          const dcomplex BWu = BW(uu,mr,0);
          const double iso_s = isoNs;
          const double iso_u = isoNu;
          double fac(1.);
          if (Config::exists("negN1675")) fac = -1.;

          for (uint al1(0); al1<4; al1++) {
            for (uint al2(0); al2<4; al2++) {
              for (uint be1(0); be1<4; be1++) {
                for (uint be2(0); be2<4; be2++) {
                  Mhad(la1,la2,mu) += 
                    fac * F_rho *
                    ubar(0,la2) *
                    ( iso_s *
                      ( (not sch) ? gamma_null
                        : ( vertex5hNrho(g1,g2,g3,-5*half,al1,al2,p1+q,mu,-k) *
                            i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                            P5h(p1+q,mr,al1,al2,be1,be2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) *
                            vertex5hNpi(g0,-5*half,be1,be2,-p1-q,q) ) )
                      + iso_u *
                      ( (not uch) ? gamma_null
                        : ( vertex5hNpi(g0,-5*half,be1,be2,p1-k,q) *
                            i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                            P5h(p1-k,mr,be1,be2,al1,al2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) *
                            vertex5hNrho(g1,g2,g3,-5*half,al1,al2,-p1+k,mu,-k) ) )
                      )
                    * u(0,la1);
                }
              }
            }
          }
        }
      }
    }
  }


  //  cerr << "W(mu,nu) = " << endl;
  Array<dcomplex,2> W(idx_lor,idx_lor);
  for (halfint mu(0); mu<4; mu++) {
    for (halfint nu(0); nu<4; nu++) {
      W(mu,nu) = 0;
      for (halfint la1 : {-half,half}) {
        for (halfint la2 : {-half,half}) {
          W(mu,nu) += Mhad(la1,la2,mu)*conj(Mhad(la1,la2,nu));
        }
      }
      //      cerr << W(mu,nu) << "  ";
    }
    //    cerr << endl;
  }

  // Calculate the production density matrix:
  Array<dcomplex,2> result(idx_s1,idx_s1);
  for (halfint la : {-_1,_0,_1}) {
    //    cerr << "epsRe(" << la << ") = " << epsRe(la) << endl;
    //    cerr << "epsIm(" << la << ") = " << epsIm(la) << endl;
    for (halfint lap : {-_1,_0,_1}) {
      result(la,lap) = 0;
      for (halfint mu(0); mu<4; mu++) {
        for (halfint nu(0); nu<4; nu++) {
          result(la,lap) += W(mu,nu)
	    * (epsRe(la)(mu)-i_*epsIm(la)(mu)) * sign_(mu)
	    * (epsRe(lap)(nu)+i_*epsIm(lap)(nu)) * sign_(nu);
        }
      }
    }
  }

  return result;
}


double DiffSigma_piN_Ndilep::B_coeff(double costh) const {
  Array<dcomplex,2> rho = rho_prod(costh);
  double SigmaT = real(rho(-_1,-_1) + rho(_1,_1));
  double SigmaL = 2.*real(rho(_0,_0));
  return (SigmaT - SigmaL)/(SigmaT + SigmaL);
}


double DiffSigma_piN_Ndilep::B_coeff_pipi(double costh) const {
  Array<dcomplex,2> rho = rho_prod(costh);
  double SigmaT = real(rho(-_1,-_1) + rho(_1,_1));
  double SigmaL = 2.*real(rho(_0,_0));
  return (SigmaL - SigmaT)/SigmaT;
}
//*/

/**
   d\sigma / dM * d\cos{\theta_\gamma*}
*/
double DiffSigma_piN_Ndilep::dsig_dMdcosth(double costh) const {
  Array<dcomplex,2> rho = rho_prod(costh);
  double trace_rho_prod = 0;
  for (halfint la : {-_1,_0,_1}) {
    trace_rho_prod += real(rho(la,la));
  }
  double MSQR = 16.*pi_/3. * e*e/M2 * trace_rho_prod;
  const double npol(2.);
  return M/(64.*POW<4>(2.*pi_)*s) * pout_abs/pin_abs * 1./npol * MSQR;
}

double DiffSigma_piN_Ndilep::dsig_pipi_dMdcosth(double costh) const {
  const double grhopipi = Config::get<double>("grho_tilde");
  Array<dcomplex,2> rho = rho_prod(costh);
  double trace_rho_prod = 0;
  for (halfint la : {-_1,_0,_1}) {
    trace_rho_prod += real(rho(la,la));
  }
  double MSQR = 8.*pi_/3. * grhopipi*grhopipi * 2. * (M*M/4.-mpi*mpi) * trace_rho_prod * POW<2>(grho/(e*M2));
  const double npol(2.);
  return M/(64.*POW<4>(2.*pi_)*s) * pout_abs/pin_abs * 1./npol * MSQR;
}


double DiffSigma_piN_Ndilep::dsig_dMdcosth_dcosthe_dphie(double costh, double costh_e, double phi_e) const {
  Array<dcomplex,2> rho = rho_prod(costh);
  double sinth_e = sqrt(1 - costh_e*costh_e);
  dcomplex MM = 
    + (1. + costh_e*costh_e) * (rho(-_1,-_1) + rho(_1,_1))
    + 2.*(1. - costh_e*costh_e) * rho(_0,_0)
    + sqrt(2.) * costh_e*sinth_e
    * (exp(i_*phi_e)*(rho(_0,_1)-rho(-_1,_0)) + exp(-i_*phi_e)*(rho(_1,_0)-rho(_0,-_1)))
    + sinth_e*sinth_e * (exp(2.*i_*phi_e)*rho(-_1,_1) + exp(-2.*i_*phi_e)*rho(_1,-_1));
  if (imag(MM)/real(MM) > 1e-3) { cerr << "warning: M^2 is not real! ( M^2 = " << MM << ")" << endl; }
  double MSQR = e*e/M2 * real(MM);
  const double npol(2.);
  return M/(64.*POW<4>(2.*pi_)*s) * pout_abs/pin_abs * 1./npol * MSQR;
}


double DiffSigma_piN_Ndilep::dsig_dM() const {
  uint nth(20);
  if (Config::exists("nth")) nth = Config::get<int>("nth");
  const double startcosth(-1.);
  const double maxcosth(1.);
  const double dcosth = (maxcosth-startcosth)/nth;
  double sum(0);
  for (double costh(startcosth+dcosth/2.); costh<maxcosth; costh+=dcosth) {
    sum += dsig_dMdcosth(costh);
  }
  return sum/nth * (maxcosth-startcosth);
}

double DiffSigma_piN_Ndilep::dsig_pipi_dM() const {
  uint nth(20);
  if (Config::exists("nth")) nth = Config::get<int>("nth");
  const double startcosth(-1.);
  const double maxcosth(1.);
  const double dcosth = (maxcosth-startcosth)/nth;
  double sum(0);
  for (double costh(startcosth+dcosth/2.); costh<maxcosth; costh+=dcosth) {
    sum += dsig_pipi_dMdcosth(costh);
  }
  return sum/nth * (maxcosth-startcosth);
}

void printHelp(char* name) {
  cerr << "usage: " << name << " [(srt|plab)=<double val>] [nth=<int val>] [sch|uch|D1232|N1440|...] [theta|costh]" << endl;
  cerr << "Order of options is arbitrary." << endl;
  exit(0);
}

int main(int argc, char** argv) {
  int nth(20);
  int nM(20);
  double srt;
  double plab;
  bool in_theta(true);
  bool theta_spect(false);
  bool theta_spect_pipi(false);
  bool mass_spect(false);
  bool dsig_dM(false);
  bool dsig_pipi_dM(false);
  bool event_gen(false);
  bool density_matrix(false);
  bool tab_density_matrix(false);
  bool diff3sig(false);

  Config::load(argc,argv);

  if (Config::exists("srt")) {
    //cerr << "srt given" << endl;
    srt = Config::get<double>("srt");
    double mn = Config::get<double>("Nucleon.mass");
    double mpi = Config::get<double>("pi_pm.mass");
    plab = sqrt(lambda(srt*srt,mn*mn,mpi*mpi))/(2.*mn);
  } else if (Config::exists("plab")) {
    // cerr << "plab given" << endl;
    plab = Config::get<double>("plab");
    double mn = Config::get<double>("Nucleon.mass");
    double mpi = Config::get<double>("pi_pm.mass");
    double Epi = sqrt(mpi*mpi + plab*plab);
    srt = sqrt(mn*mn + mpi*mpi + 2.*mn*Epi);
  } else {
    cerr << "srt or plab NOT given" << endl;
    printHelp(argv[0]);
    exit(0);
  }
  if (Config::exists("dsig_dM")) { dsig_dM=true; }
  if (Config::exists("dsig_pipi_dM")) { dsig_pipi_dM=true; }
  if (Config::exists("theta_spect")) { theta_spect=true; }
  if (Config::exists("theta_spect_pipi")) { theta_spect_pipi=true; }
  if (Config::exists("mass_spect")) { mass_spect=true; }
  if (Config::exists("event_gen")) { event_gen=true; }
  if (Config::exists("density_matrix")) { density_matrix=true; }
  if (Config::exists("tab_density_matrix")) { tab_density_matrix=true; }
  if (Config::exists("diff3sig")) { diff3sig=true; }
  if (not (dsig_dM or dsig_pipi_dM or theta_spect or theta_spect_pipi or mass_spect or
	   event_gen or density_matrix or tab_density_matrix or diff3sig)) {
    cerr << "Specify one of [dsig_dM|dsig_pipi_dM|theta_spect|theta_spect_pipi|mass_spect|event_gen|density_matrix|tab_density_matrix|diff3sig]!" << endl;
    exit(0);
  }

  //  PR(tab_density_matrix);

  if (Config::exists("nth")) nth = Config::get<int>("nth");
  if (Config::exists("nM")) nM = Config::get<int>("nM");
  if (Config::exists("sch")) uch=false;
  if (Config::exists("uch")) sch=false;
  if (sch==false and uch==false) sch=uch=true;

  if (Config::exists("D1232")) D1232=true;
  if (Config::exists("N1440")) N1440=true;
  if (Config::exists("N1535")) N1535=true;
  if (Config::exists("N1520")) N1520=true;
  if (Config::exists("N1650")) N1650=true;
  if (Config::exists("N1675")) N1675=true;
  if (Config::exists("N1680")) N1680=true;
  if (Config::exists("D1600")) D1600=true;
  if (Config::exists("D1620")) D1620=true;

  if (Config::exists("theta")) in_theta=true;
  if (Config::exists("costh")) in_theta=false;

  if ( not (D1232 or N1440 or N1535 or N1520 or N1650 or N1675 or N1680 or D1600 or D1620) ) {
    cerr << "All channels switched off" << endl;
    exit(0);
  }

  //============================================================
  //*
  cout << "# Parameters:" << endl;
  Config::list(cout,"# ");
  cout << "# ======================================================================\n" << endl;
  cout << "#" << endl;;

  //------------------------------------------------------------------------
  if (density_matrix) {
    double M;
    try {
      M = Config::get<double>("M");
    } catch (...) {
      cerr << "dilepton mass M has to be specified!" << endl;
      exit(0);
    }
    double costh;
    if (Config::hasValue("theta")) {
      costh = cos(Config::get<double>("theta"));
    } else if (Config::hasValue("costh")) {
      costh = Config::get<double>("costh");
    } else {
      cerr << "costh or theta must be specified!" << endl;
      exit(0);
    }
    DiffSigma_piN_Ndilep dsigma(srt,M);
    cout << "density matrix:" << endl;
    Array<dcomplex,2> RHO = dsigma.rho_prod(costh);
    for (halfint la1(-_1); la1<=_1; la1++) { 
      cout << "( ";
      for (halfint la2(-_1); la2<=_1; la2++) { 
	cout << setw(20) << RHO(la1,la2) << " ";
      }
      cout << ")" << endl;
    }

    std::ofstream of;
    if (Config::exists("dmfile")) {
      string dmfile = Config::get<string>("dmfile");
      of.open(dmfile,std::ofstream::app);
      if (not of) { of.open(dmfile,std::ofstream::trunc); }
    }
    std::streambuf * buf;
    if(of) {
      buf = of.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream dmout(buf);

    dmout << "\n# Real part tabbed:\n";
    dmout << "#    sqrt(s)    M_{e+e-}  cos(thetaCM)";
    for (halfint la1(_1); la1>=-_1; la1--) {
      for (halfint la2(_1); la2>=-_1; la2--) {
        dmout << "  rho(" << setw(2) << la1 << "," << setw(2) << la2 << ")";
      }
    }
    dmout << endl;
    dmout << "# " << (N1520 ? "N1520" : "") << (N1440 ? (Config::exists("negN1440") ? "-N1440" : "+N1440") : "") << endl;
    dmout << setw(12) << srt << setw(12) << M << setw(14) << costh;
    for (halfint la1(_1); la1>=-_1; la1--) {
      for (halfint la2(_1); la2>=-_1; la2--) {
        double rhoRe = real(RHO(la1,la2));
        double rhoIm = imag(RHO(la1,la2));
        if (fabs(rhoRe) < 1.e-10) { rhoRe = 0; }
        if (fabs(rhoIm) > 1.e-10) {
          cerr << "Warning - RHO(" << la1 << "," << la2 << ") not real: " << RHO(la1,la2) << endl;
          cerr << "    ( "  << (N1520 ? "N1520" : "") << (N1440 ? (Config::exists("negN1440") ? "-N1440" : "+N1440") : " )")
               << "  srt=" << srt << "  M=" << M << "  costh=" << costh << " )" << endl;
        }
        
        dmout << "  " << setw(10) << rhoRe;
      }
    }
    dmout << "\n" << endl;
    
    dcomplex TrRHO(0);
    dcomplex TrRHO2(0); 
    for (halfint la1(-_1); la1<=_1; la1++) {  
      TrRHO += RHO(la1,la1);
      for (halfint la2(-_1); la2<=_1; la2++) {
        TrRHO2 += RHO(la1,la2) * RHO(la2,la1);
      } 
    }
    cout << "trace(rho) = " << TrRHO << endl;
    cout << "trace(rho^2) = " << TrRHO2 << endl;
    cout << "normalized trace(rho^2) = " << TrRHO2/(TrRHO*TrRHO) << endl;
    cout << "\nlambda_theta     = " << (RHO(-_1,-_1)+RHO(_1,_1)-2.*RHO(_0,_0))
      /(RHO(-_1,-_1)+RHO(_1,_1)+2.*RHO(_0,_0)) << endl;
    cout << "\nlambda_theta_phi = " << real(RHO(_0,_1)-RHO(-_1,_0))/(RHO(-_1,-_1)+RHO(_1,_1)+2.*RHO(_0,_0))
         << endl;
    cout << "\nlambda_phi       = " << real(RHO(-_1,_1))/(RHO(-_1,-_1)+RHO(_1,_1)+2.*RHO(_0,_0)) << endl;
  //------------------------------------------------------------------------
  } else if (tab_density_matrix) {
    double dM(10*MeV);
    if(Config::exists("dM")) {
      dM = Config::get<double>("dM");
    }
    const double mn = Config::get<double>("Nucleon.mass");
    double Mmax(srt-mn);
    cout << "#" << setw(5) << "M" << " " << setw(6) << "costh" << " ";
    for (halfint la1 : {_1,_0,-_1}) {
      for (halfint la2 : {_1,_0,-_1}) {
	cout << setw(18) << "rho(" << setw(2) << la1 << "," << setw(2) << la2 << ")           ";
      }
    }
    cout << endl;
    for (double M(dM); M<Mmax; M+=dM) {
      double dcosth = 2./nth;
      for (double costh(-1.); costh <= 1.; costh += dcosth) {
	if (fabs(costh) < 1e-8) { costh=0; }
	DiffSigma_piN_Ndilep dsigma(srt,M);
	cout << setw(6) << M << " " << setw(6) << costh << "  ";
	Array<dcomplex,2> RHO = dsigma.rho_prod(costh);
	for (halfint la1 : {_1,_0,-_1}) {
	  for (halfint la2 : {_1,_0,-_1}) {
	    double RE = real(RHO(la1,la2)) * POW<2>(e/(M*M));
	    double IM = imag(RHO(la1,la2)) * POW<2>(e/(M*M));
	    if (fabs(IM/RE) < 1e-10) { 
	      IM = 0; 
	    } else if (fabs(RE/IM) < 1e-10) { 
	      RE = 0;
	    }
	    cout << "( " << setw(13) << RE << ", " << setw(13) << IM << " )   ";
	  }
	}
	cout << endl;
      }
      cout << endl;
    }
  //------------------------------------------------------------------------
  } else if (mass_spect) {
    double costh;
    if (Config::hasValue("theta")) {
      costh = cos(Config::get<double>("theta"));
    } else if (Config::hasValue("costh")) {
      costh = Config::get<double>("costh");
    } else {
      cerr << "costh or theta must be specified!" << endl;
      exit(0);
    }
    const double mn = Config::get<double>("Nucleon.mass");
    const double startM(0);
    const double maxM(srt-mn);
    const double dM = (maxM-startM)/nM;
    PR(startM) PR(maxM) PR(dM);
    //exit(0);
    for (double M(startM+dM); M<maxM; M+=dM) {
      DiffSigma_piN_Ndilep dsigma(srt,M);
      cout << setw(13) << M;
      cout << " " << setw(15) << mub(dsigma.dsig_dMdcosth(costh));
      cout << " " << setw(15) << dsigma.B_coeff(costh);
      cout << endl;
    }
  //------------------------------------------------------------------------
  } else if (dsig_dM) {
    double M;
    try {
      M = Config::get<double>("M");
    } catch (...) {
      cerr << "dilepton mass M has to be specified!" << endl;
      exit(0);
    }
    DiffSigma_piN_Ndilep dsigma(srt,M);
    cout << "dsigma/dM [mub/GeV] = " << mub(dsigma.dsig_dM()) << endl;
  //------------------------------------------------------------------------
  } else if (diff3sig) {
    double M;
    try {
      M = Config::get<double>("M");
    } catch (...) {
      cerr << "dilepton mass M has to be specified!" << endl;
      exit(0);
    }
    cerr << "# differential cross section dsigma/dM dcosth dcosth_e dphi_e" << endl;
    cerr << "#" << setw(9) << "costh" << setw(10) << "costh_e" << setw(12) << "phi_e" << setw(15) << "diff.sigma" << endl;

    DiffSigma_piN_Ndilep dsigma(srt,M);
    double dcosth = 2./nth;
    double dnphi = 2./nth;

    for (double costh(-1.); costh <= 1.; costh += dcosth) {
      //------------------------------------------------------------------------
      cout << "density matrix for cos(theta) = " << costh << endl;
      Array<dcomplex,2> RHO = dsigma.rho_prod(costh);
      for (halfint la1(-_1); la1<=_1; la1++) { 
	cout << "( ";
	for (halfint la2(-_1); la2<=_1; la2++) { 
	  cout << setw(20) << RHO(la1,la2) << " ";
	}
	cout << ")" << endl;
      }
      //------------------------------------------------------------------------
      for (double costh_e(-1.); costh_e <= 1.; costh_e += dcosth) {
	for (double nphi(0.); nphi <= 2.; nphi += dnphi) {
	  double phi_e = nphi*pi_;
	  cerr << setw(10) << costh << setw(10) << costh_e << setw(7) << nphi << "*pi  " 
	       << setw(15) << dsigma.dsig_dMdcosth_dcosthe_dphie(costh, costh_e, phi_e) << endl;
	}
	cerr << endl;
      }
      cerr << "------------------------------------------------------------------------\n" << endl;
    }
  //------------------------------------------------------------------------
  } else if (dsig_pipi_dM) {
    double M;
    try {
      M = Config::get<double>("M");
    } catch (...) {
      cerr << "dilepton mass M has to be specified!" << endl;
      exit(0);
    }
    DiffSigma_piN_Ndilep dsigma(srt,M);
    cout << "dsigma/dM [mub/GeV] = " << mub(dsigma.dsig_pipi_dM()) << endl;
  //------------------------------------------------------------------------
  } else if (theta_spect) {
    double M;
    try {
      M = Config::get<double>("M");
    } catch (...) {
      cerr << "dilepton mass M has to be specified!" << endl;
      exit(0);
    }
    cout << "# Differential cross section and anisotropy coefficient (B) of dilepton production in pi + N\n#" << endl;
    cout << setw(12) << (mass_spect ? "M (GeV)" : (in_theta ? "theta" : "costh")) << " ";
    cout << setw(15) << "dsig (mub/GeV)" << " " << setw(15) << "B" << endl;
    cout << "#" << endl;
    DiffSigma_piN_Ndilep dsigma(srt,M);
    if (in_theta) {
      const double starttheta(0);
      const double maxtheta(pi_);
      const double dtheta = (maxtheta-starttheta)/nth;
      for (double theta(starttheta); theta<=maxtheta; theta+=dtheta) {
        double costh = cos(theta);
        //        PL(costh) PL(costh_k) PR(phi_k);
        cout << setw(13) << theta;
        cout << " " << setw(15) << mub(dsigma.dsig_dMdcosth(costh));
        cout << " " << setw(15) << dsigma.B_coeff(costh);
        cout << endl;
	// cout << "density matrix:" << endl;
	// Array<dcomplex,2> RHO = dsigma.rho_prod(costh);
	// for (int la1(-1); la1<=1; la1++) { 
	//   cout << "( ";
	//   for (int la2(-1); la2<=1; la2++) { 
	//     cout << setw(20) << RHO(la1,la2) << " ";
	//   }
	//   cout << ")" << endl;
	// }
      } // costh
    } else {
      const double startcosth(-1.);
      const double maxcosth(1.);
      const double dcosth = (maxcosth-startcosth)/nth;
      for (double costh(startcosth); costh<=maxcosth+dcosth/2.; costh+=dcosth) {
        //        PL(costh) PL(costh_k) PR(phi_k);
	if (costh>1) { costh=1.; }
        cout << setw(13) << costh;
        cout << " " << setw(15) << mub(dsigma.dsig_dMdcosth(costh));
        cout << " " << setw(15) << dsigma.B_coeff(costh);
        cout << endl;
      } // costh
    }
  //------------------------------------------------------------------------
  } else if (theta_spect_pipi) {
    double M;
    try {
      M = Config::get<double>("M");
    } catch (...) {
      cerr << "dilepton mass M has to be specified!" << endl;
      exit(0);
    }
    cout << "# Differential cross section and anisotropy coefficient (B) of pion pair production in pi + N\n#" << endl;
    cout << setw(12) << (mass_spect ? "M (GeV)" : (in_theta ? "theta" : "costh")) << " ";
    cout << setw(15) << "dsig (mub/GeV)" << " " << setw(15) << "B" << endl;
    cout << "#" << endl;
    DiffSigma_piN_Ndilep dsigma(srt,M);
    if (in_theta) {
      const double starttheta(0);
      const double maxtheta(pi_);
      const double dtheta = (maxtheta-starttheta)/nth;
      for (double theta(starttheta); theta<=maxtheta; theta+=dtheta) {
        double costh = cos(theta);
        //        PL(costh) PL(costh_k) PR(phi_k);
        cout << setw(13) << theta;
        cout << " " << setw(15) << mub(dsigma.dsig_pipi_dMdcosth(costh));
        cout << " " << setw(15) << dsigma.B_coeff_pipi(costh);
        cout << endl;
	// cout << "density matrix:" << endl;
	// Array<dcomplex,2> RHO = dsigma.rho_prod(costh);
	// for (int la1(-1); la1<=1; la1++) { 
	//   cout << "( ";
	//   for (int la2(-1); la2<=1; la2++) { 
	//     cout << setw(20) << RHO(la1,la2) << " ";
	//   }
	//   cout << ")" << endl;
	// }
      } // costh
    } else {
      const double startcosth(-1.);
      const double maxcosth(1.);
      const double dcosth = (maxcosth-startcosth)/nth;
      for (double costh(startcosth); costh<=maxcosth+dcosth/2.; costh+=dcosth) {
        //        PL(costh) PL(costh_k) PR(phi_k);
	if (costh>1) { costh=1.; }
        cout << setw(13) << costh;
        cout << " " << setw(15) << mub(dsigma.dsig_pipi_dMdcosth(costh));
        cout << " " << setw(15) << dsigma.B_coeff_pipi(costh);
        cout << endl;
      } // costh
    }
  //------------------------------------------------------------------------
  } else if (event_gen) {
    uint Nevent;
    try {
      Nevent = Config::get<int>("Nevent");
    } catch (...) {
      cerr << "No. of generated events, Nevent has to be specified! (e.g. Nevent=1000000)" << endl;
      exit(0);
    }
    if (Config::exists("seed")) {
      uint seed = Config::get<int>("seed");
      rn.setSeed(seed);
    }
    cout << "# " << setw(8) << "plab" << setw(10) << "particle"
	 << setw(14) << "p_x" << setw(14) << "p_y" << setw(14) << "p_z" << setw(14) << "diff. sigma" << endl;
    cout << "# " << setw(8) << "(GeV/c)" << setw(10) << "type  "
	 << setw(14) << "(GeV/c)" << setw(14) << "" << setw(14) << "" << setw(14) << "(mub/GeV)" << endl;
    cout << "#" << endl;
    uint nbins(0);
    if (Config::exists("nbins")) {
      nbins = Config::get<int>("nbins");
    }
    Histogram1d hist_dsig_dcosth(-1.,1,nbins);
    for (uint n(0); n<Nevent; n++) {
      double mn = Config::get<double>("Nucleon.mass");
      double mpi = Config::get<double>("pi_pm.mass");
      double Mmax = srt - mn - 1.*MeV;
      double M = Mmax*rn();
      double costh = 2.*rn() - 1.;
      double sinth = sqrt(1.-costh*costh);
      double phi = 2.*pi_*rn();
      double costh_e = 2.*rn() - 1.;
      double sinth_e = sqrt(1.-costh_e*costh_e);
      double phi_e = 2.*pi_*rn();
      DiffSigma_piN_Ndilep dsigma(srt,M);
      double diffsig = dsigma.dsig_dMdcosth_dcosthe_dphie(costh,costh_e,phi_e);
      double pout_abs = sqrt(lambda(srt*srt,M*M,mn*mn))/(2.*srt);
      double k0 = (srt*srt+M*M-mn*mn)/(2.*srt);
      double EN_out = (srt*srt-M*M+mn*mn)/(2.*srt);

      double pin_abs = sqrt(lambda(srt*srt,mpi*mpi,mn*mn))/(2.*srt);
      double EN_in = (srt*srt-mpi*mpi+mn*mn)/(2.*srt);
      double Epi_in = (srt*srt+mpi*mpi-mn*mn)/(2.*srt);

      double ch = k0/M;
      double sh = -pout_abs/M;
      FourTensor Boost(ch, 0,  0,  sh,
		       0, -1,  0,  0,
		       0,  0, -1,  0,
		      -sh, 0,  0, -ch );

      double chlab = EN_in/mn;
      double shlab = -pin_abs/mn;
      FourTensor Boostlab(chlab, 0,  0, shlab,
		          0,    -1,  0,  0,
		          0,     0, -1,  0,
		          -shlab, 0,  0, -chlab );

      FourTensor Rot( 1,       0,      0,      0,
		      0,  -costh,      0, -sinth,
		      0,       0,     -1,      0,
		      0,   sinth,      0, -costh );

      FourTensor Rot_phi( 1,       0,          0,    0,
		          0,  -cos(phi),  -sin(phi), 0,
		          0,   sin(phi),  -cos(phi), 0,
		          0,       0,          0,   -1 );

      // check boost to lab system:
      // FourVector p1 = Boostlab*FourVector(EN_in,0,0,-pin_abs);
      // FourVector qi = Boostlab*FourVector(Epi_in,0,0,pin_abs);
      FourVector p1 = FourVector(EN_in,0,0,-pin_abs);
      FourVector qi = FourVector(Epi_in,0,0,pin_abs);
      FourVector p1lab = Boostlab*p1;
      FourVector qilab = Boostlab*qi;
      // PR(p1);
      // PR(qi);
      // PR(p1lab);
      // PR(qilab);

      // PR(Boostlab);
      //exit(0);

      ThreeVector ek1(sinth_e*cos(phi_e),sinth_e*sin(phi_e),costh_e);
      FourVector k1_0(M/2.,M/2.*ek1);
      FourVector k2_0(M/2.,-M/2.*ek1);
      FourVector k1 = Boostlab*Rot_phi*Rot*Boost*k1_0;
      FourVector k2 = Boostlab*Rot_phi*Rot*Boost*k2_0;

      ThreeVector kk = pout_abs*ThreeVector(sinth,0,costh);
      FourVector p2 = Boostlab*Rot_phi*FourVector(EN_out,-kk);

      FourVector Ptot = k1+k2+p2;
      // cerr << "total fourmom: " << Ptot << endl;

      //PR(pin_abs);

      FourTensor BoostCM(chlab,  0,  0, -shlab,
		          0,    -1,  0,  0,
		          0,     0, -1,  0,
		         shlab,  0,  0, -chlab );

      FourVector PtotCM = BoostCM*Ptot;
      // cerr << "total CM fourmom:  " << PtotCM << endl;

      cout << setw(10) << plab << setw(10) << "e-"
	   << setw(14) << k1[1] << setw(14) << k1[2] << setw(14) << k1[3] << setw(14) << mub(diffsig) << endl;
      cout << setw(10) << plab << setw(10) << "e+"
	   << setw(14) << k2[1] << setw(14) << k2[2] << setw(14) << k2[3] << setw(14) << mub(diffsig) << endl;
      // cout << setw(10) << plab << setw(10) << "n"
      // 	   << setw(14) << p2[1] << setw(14) << p2[2] << setw(14) << p2[3] << setw(14) << mub(diffsig) << endl;
      // cout << setw(10) << plab << setw(10) << "sum"
      // 	   << setw(14) << k1[1]+k2[1]+p2[1] << setw(14) << k1[2]+k2[2]+p2[2] << setw(14) << k1[3]+k2[3]+p2[3] << setw(14) << mub(diffsig) << endl;

      // histogram of dsigma/dcosth:
      if (nbins>0) {
	FourVector K = k1+k2;
	FourVector KCM = BoostCM*K;
	ThreeVector kk = KCM.spacial();
	ThreeVector ez(0,0,1);
	double costh = kk*ez/sqrt(kk*kk);
	//cerr << "costh = " << costh << endl;
	hist_dsig_dcosth.add(costh,mub(diffsig));
      }
    }

    if (nbins>0) {
      cout << "Histogram of dsigma/dcosth:" << endl;
      hist_dsig_dcosth.print(cout,"cos(theta)","",15,15);
    }
  }
  return 1;
}
