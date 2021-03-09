#ifndef KERNELS_HPP
#define KERNELS_HPP

#include <vector>

#include "QCD.hpp"
#include "Random.hpp"

class Kernels
{
public:
  std::vector<int> flavs;
  Random* ran;
public:
  Kernels(std::vector<int> fl, Random* random)
  : flavs{fl}, ran{random}
  {}
  virtual ~Kernels() = 0;
  virtual double Value(const double z, const double y)       const = 0;
  virtual double Estimate(const double z)                    const = 0;
  virtual double Integral(const double zm, const double zp)  const = 0;
  virtual double GenerateZ(const double zm, const double zp) const = 0;
  inline friend std::ostream & operator<<(std::ostream & os, Kernels& kern){
    os << "SF : " << kern.flavs[0] << " -> " << kern.flavs[1] << " , " << kern.flavs[2] << "\n";
    return os;
  }
};

inline Kernels::~Kernels() {}

class Pqq : public Kernels
{
public:
  Pqq(const int& fl, Random *random) : Kernels({fl, fl, 21}, random) {}
  ~Pqq() {}
  inline double Value(const double z, const double y) const override
  {
    return QCD::CF * (2./(1.-z*(1.-y))-(1.+z));
  }

  inline double Estimate(const double z) const override
  {
    return QCD::CF * 2./(1.-z);
  }

  inline double Integral(const double zm, const double zp) const override
  {
    return QCD::CF * 2. * log((1. - zm)/(1. - zp));
  }

  inline double GenerateZ(const double zm, const double zp) const override
  {
    return 1. + (zp - 1.) * pow((1. - zm)/(1. - zp), (*ran)());
  }
};

class Pgg : public Kernels
{
public:
  Pgg(Random *random) : Kernels({21,21,21}, random) {}
  ~Pgg() {}
  inline double Value(const double z, const double y) const override
  {
    return QCD::CA /2. * (2./(1.-z*(1.-y))-2.+z*(1.-z));
  }

  inline double Estimate(const double z) const override
  {
    return QCD::CA/(1.-z);
  }

  inline double Integral(const double zm, const double zp) const override
  {
    return QCD::CA * log((1.-zm)/(1.-zp));
  }

  inline double GenerateZ(const double zm, const double zp) const override
  {
    return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), (*ran)());
  }
};

class Pgq : public Kernels
{
public:
  Pgq(const int& fl, Random *random) : Kernels({21,fl,-fl}, random) {}
  ~Pgq() {}
  inline double Value(const double z, const double y) const override
  {
    return QCD::TR/2. *(1.-2.*z*(1.-z));
  }

  inline double Estimate(const double z) const override
  {
    return QCD::TR/2.;
  }

  inline double Integral(const double zm, const double zp) const override
  {
    return QCD::TR/2. * (zp -zm);
  }

  inline double GenerateZ(const double zm, const double zp) const override
  {
    return zm + (zp - zm) * (*ran)();
  }
};

#endif