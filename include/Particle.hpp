#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "Rivet/Rivet.hh"

typedef std::pair<int,int> Colour;
typedef std::pair<int, Rivet::FourMomentum> Particle_Info;
typedef std::pair<Particle_Info, Colour> Particle_Data;

class Particle
{
private:
  Particle_Data _pd;
public:
  Particle(const int& fl, const Rivet::FourMomentum& fv,
           const Colour cl = std::make_pair<int,int>(0,0))
  : _pd{{fl,fv},cl}
  { }
  ~Particle() {}
  /// Access members
  inline int GetFlavour()                  const {return _pd.first.first;}
  inline Rivet::FourMomentum GetMomentum() const {return _pd.first.second;}
  inline Colour GetColour()                const {return _pd.second;}
  /// Set members
  inline void SetFlavour(const int& fl)                  {_pd.first.first = fl;}
  inline void SetMomentum(const Rivet::FourMomentum& fv) { _pd.first.second = fv; }
  inline void SetColour(const Colour& cl)                { _pd.second = cl; }

  /// static members
  /// Boost pb wrt to pa
  inline static Rivet::FourMomentum Boost(const Rivet::FourMomentum& pa,
                                          const Rivet::FourMomentum& pb)
  {
    double rsq {pa.mass()};
    double v0 {(pa.E() * pb.E() - pa.px() * pb.px() - pa.py() * pb.py() - pa.pz() * pb.pz() )/rsq};
    double c1 {(pb.E() + v0)/(rsq + pa.E())};
    return Rivet::FourMomentum{v0,
                               pb.px() - c1 * pa.px(),
                               pb.py() - c1 * pa.py(),
                               pb.pz() - c1 * pa.pz()};
  }

  inline static Rivet::FourMomentum BoostBack(const Rivet::FourMomentum& pa,
                                              const Rivet::FourMomentum& pb)
  {
    double rsq {pa.mass()};
    double v0 {(pa.E() * pb.E() + pa.px() * pb.px() + pa.py() * pb.py() + pa.pz() * pb.pz() )/rsq};
    double c1 {(pb.E() + v0)/(rsq + pa.E())};
    return Rivet::FourMomentum{v0,
                               pb.px() + c1 * pa.px(),
                               pb.py() + c1 * pa.py(),
                               pb.pz() + c1 * pa.pz()};
  }

  /// overloaded operators
  inline friend std::ostream &operator<<(std::ostream &os, const Particle &p)
  {
    os << p.GetFlavour() << " : " << p.GetMomentum() << " : " << p.GetColour();
    return os;
  }
};


#endif