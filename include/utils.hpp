#ifndef UTILS_HPP
#define UTILS_HPP

#include <complex>

using dcomplex = std::complex<double>;

const dcomplex i_(0,1);
bool isnan(const dcomplex&);
double isreal(const dcomplex&, double eps=1e-6);

const double pi_ = 4.0*atan(1.0);

#define DEFPOW(x)                               \
  const double x ## 2(x*x);                     \
  const double x ## 3(x ## 2*x);                \
  const double x ## 4(x ## 3*x);                \
  const double x ## 5(x ## 4*x);                \
  const double x ## 6(x ## 5*x);

#define PR(x) std::cerr << #x << " = " << (x) << std::endl;
#define PL(x) std::cerr << #x << " = " << (x);
#define LOG(x) std::cout << "# " << std::setiosflags(std::ios::left) << setw(10) << #x  \
  << " = " << std::setiosflags(std::ios::right) << (x) << std::endl;
#define LOGC(x,comment) std::cout << "# " << std::setiosflags(std::ios::left) << setw(10) << #x  \
  << " = " << std::setiosflags(std::ios::right) << (x) << "      // " << comment << std::endl;

template <unsigned int n>
double POW(double x) {
  return x*POW<n-1>(x);
}

template <>
double POW<0>(double x);

template <unsigned int n>
dcomplex POW(dcomplex x) {
  return x*POW<n-1>(x);
}

template <>
dcomplex POW<0>(dcomplex x);

double lambda(double x, double y, double z);

template <typename T>
int delta_(T i, T j) {
  return (i==j) ? 1 : 0;
}


#endif // UTILS_HPP
