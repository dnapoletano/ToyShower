#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cmath>
#include <random>

#include "Particle.hpp"
#include "QCD.hpp"
#include "Random.hpp"

#include "Rivet/Rivet.hh"

struct EWParameters
{
  double mz2, gz2, alpha0, sin2tw, qe, ae;
  EWParameters()
  : mz2{91.1876*91.1876}, gz2{2.4952*2.4952},
    alpha0{1./128.802}, sin2tw{0.22293}, qe{-1},
    ae{-0.5}
  {}
};

struct EventInfo {
  long int EvtNumber;
  double dxs, lome;
  std::vector<Particle> Particles;
  inline friend std::ostream& operator<<(std::ostream& os, EventInfo& evt){
    os << "XS       : " << evt.dxs<<"\n";
    os << "MEWeight : " << evt.lome << "\n";
    for(const auto& p: evt.Particles){
      os << p << " : " << &p << "\n";
    }
    return os;
  }
};

class Matrix
{
private:
  EWParameters _ewparams;
  const double _ecms;
  Random* ran;

public:
  Matrix(const double& ecms, Random* random);
  ~Matrix() {}

  double ME2(const int& flav, const double& s, const double& t) const;
  EventInfo GeneratePoint();
  inline void SetEWParameters(const EWParameters& ewp) {_ewparams = ewp;}
};

#endif