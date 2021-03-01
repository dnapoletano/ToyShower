#include "QCD.hpp"

class Shower
{
private:
  double _tmax, _tmin;
  QCD::AlphaS _alphaS;
public:
  Shower(const QCD::AlphaS& alphaS, const double t0,
         const double t1);
  ~Shower();
};