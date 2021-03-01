#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cmath>
#include <random>

#include "QCD.hpp"
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
  Rivet::Particles Particles;
  double dxs;
};

class Matrix
{
private:
  EWParameters _ewparams;
  const double _ecms;

  std::mt19937_64 re;
  std::uniform_real_distribution<double> urng;

public:
  Matrix(const double& ecms);
  ~Matrix() {}

  double ME2(const int& flav, const double& s, const double& t) const;
  EventInfo GeneratePoint();
  inline void SetEWParameters(const EWParameters& ewp) {_ewparams = ewp;}
};

#endif