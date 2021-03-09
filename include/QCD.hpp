#ifndef QCD_HPP
#define QCD_HPP

#include <cmath>

namespace QCD {
  constexpr double NC = 3.0;
  constexpr double TR = 1./2.;
  constexpr double CA = NC;
  constexpr double CF = (NC*NC-1.)/2./NC;
}

using namespace QCD;

class AlphaS
{
private:
  size_t _order;
  double _mz2, _mb2, _mc2;
  double _asmz, _asmb, _asmc;
public:
  AlphaS(const size_t order=1, const double MZ=91.1876,
         const double alphaMZ=0.118, const double mb=4.92,
         const double mc=1.42)
  : _order{order}, _mz2{MZ*MZ}, _mb2{mb*mb}, _mc2{mc*mc}, _asmz{alphaMZ},
    _asmb{(*this)(_mb2)}, _asmc{(*this)(_mc2)}
  {}

  ~AlphaS() {}

  inline static double beta0(const size_t nf) {return 11./6.*CA -2./3.*TR*nf;}
  inline static double beta1(const size_t nf) {return 17./6.*CA*CA - (5./3.*CA+CF)*TR*nf;}
  inline double as0(const double& t) {
    double tref{}, asref{}, b0{};
    if(t >= _mb2){
      tref  = _mz2;
      asref = _asmz;
      b0    = beta0(5)/2./M_PI;
    } else if (t >= _mc2){
      tref  = _mb2;
      asref = _asmb;
      b0    = beta0(4)/2./M_PI;
    } else {
      tref  = _mc2;
      asref = _asmc;
      b0    = beta0(3)/2./M_PI;
    }
    return 1./(1./asref + b0 * log(t/tref));
  }

  inline double as1(const double& t) {
    double tref{}, asref{}, b0{}, b1{};
    if(t >= _mb2){
      tref  = _mz2;
      asref = _asmz;
      b0    = beta0(5)/2./M_PI;
      b1    = beta1(5)/2./M_PI/2./M_PI;
    } else if (t >= _mc2){
      tref  = _mb2;
      asref = _asmb;
      b0    = beta0(4)/2./M_PI;
      b1    = beta1(4)/2./M_PI/2./M_PI;
    } else {
      tref  = _mc2;
      asref = _asmc;
      b0    = beta0(3)/2./M_PI;
      b1    = beta1(3)/2./M_PI/2./M_PI;
    }
    const double w {1. + b0 * asref *log(t/tref)};
    return asref/w * (1. - b1/b0 *asref * log(w)/w);
  }

  inline double operator()(const double& t)
  {
    return (_order == 0) ? as0(t) : as1(t);
  }
};

#endif