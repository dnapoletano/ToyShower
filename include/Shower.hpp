#ifndef SHOWER_HPP
#define SHOWER_HPP

#include <memory>
#include <vector>

#include "Particle.hpp"

struct DipoleInfo{
  typedef std::vector<Particle> Partons;
  Partons::iterator split, spect;
  std::vector<std::unique_ptr<class Kernels> >::iterator selected;
  double m2, zp;
  inline friend std::ostream& operator<<(std::ostream& os, const DipoleInfo& di){
    os << " Splitter  : " << *(di.split) << "\n"
       << " Spectator : " << *(di.spect) << "\n"
       << " Spllitng  : " << di.selected->get() << "\n"
       << " m2 : " << di.m2 << "\n"
       << " zp : " << di.zp;
    return os;
  }
};

class Shower
{
private:
  int _c {0};
  double _tEnd, _tActual, _alphaSMax;
  class Random* _ran;
  class AlphaS* _alphaS;
  std::vector<std::unique_ptr<class Kernels> > _kernels;
  DipoleInfo _dipole;
public:
  typedef std::pair<int,int>  Colour;
  typedef std::vector<Colour> Colours;
  typedef std::vector<Particle> Partons;
  /// pass reference to alphaS class, random class and
  /// shower stopping scale, t0
  Shower(class AlphaS* alphaS, class Random* ran,
         const double t0);
  ~Shower() {}

  Rivet::FourMomenta MakeKinematics(const double& z, const double &y,
                                    const double& phi, const Rivet::FourMomentum& pijt,
                                    const Rivet::FourMomentum& pkt) const;
  Colours MakeColours(const std::vector<int>& flavs, const Colour& colij,
                      const Colour& colk);
  void Run(class EventInfo &evt, const double t);
  void GeneratePoint(class EventInfo& evt);
  void SelectSplitSpect(class EventInfo& evt, double& t);
  static bool CheckEvent (class EventInfo &evt);
  static bool ColourConnected(const Particle& pa, const Particle& pb);
};

#endif