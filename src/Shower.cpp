#include "Shower.hpp"

#include "Kernels.hpp"
#include "Matrix.hpp"
#include "Random.hpp"
#include "QCD.hpp"

#include <iomanip>

Shower::Shower(AlphaS* alphaS, Random* ran,
               const double t0)
: _c{0}, _alphaS{alphaS}, _ran{ran}, _tEnd{t0}, _tActual{-1.0}
{
   _alphaSMax = (*_alphaS)(_tEnd);
  /// load Pqq (q -> q g) in kernels
  for(const auto& fl_i : {-5,-4,-3,-2,-1,1,2,3,4,5}){
    _kernels.emplace_back(new Pqq{fl_i, ran});
  }
  /// load Pgq (g -> q bar{q}) in kernels
  for(const auto& fl_i : {1,2,3,4,5}){
    _kernels.emplace_back(new Pgq{fl_i, ran});
  }
  /// load Pgg (g -> gg) kernel
  _kernels.emplace_back(new Pgg{ran});
  _kernels.shrink_to_fit();
}

Rivet::FourMomenta Shower::MakeKinematics(const double& z, const double &y,
                                          const double& phi, const Rivet::FourMomentum& pijt,
                                          const Rivet::FourMomentum& pkt) const
{
  //std::cout << std::setprecision(16);
  Rivet::FourMomentum Q {pijt + pkt};
  //std::cout << "Q : " << Q << std::endl;
  /// kt in the decay frame
  const double rkt {sqrt(Q.mass2() * y * z * (1.-z))};
  //std::cout << "rkt : " << rkt << std::endl;
  Rivet::ThreeVector vkt1 {Rivet::cross(pijt.vector3(),pkt.vector3())};
 // std::cout << " vkt1 : " << vkt1 << std::endl;
  if(vkt1.mod() < 1.e-6){
    vkt1 = Rivet::cross(pijt.vector3(),{1.,0.,0.});
  }
  //std::cout << " vkt1 : " << vkt1 << " , " << cos(phi) << " , "<< vkt1.mod() << std::endl;
  Rivet::FourMomentum kt1 {0.0,vkt1[0],vkt1[1],vkt1[2]};
  kt1 *= (rkt * cos(phi) / kt1.vector3().mod());
 // std::cout << "kt1 : " << kt1 << std::endl;
  /// Boost to CMS
  Rivet::LorentzTransform ToCMS {};
  ToCMS.setBetaVec(Q.betaVec());
  Rivet::ThreeVector vkt2CMS {
    Rivet::cross(ToCMS.transform(pijt).vector3(), kt1.vector3())
  };
  vkt2CMS *= rkt * sin(phi) / vkt2CMS.mod();
  //std::cout << "kt2CMS : " << vkt2CMS << std::endl;
  Rivet::FourMomentum kt2 {
    ToCMS.inverse().transform(
      Rivet::FourMomentum{0.0,vkt2CMS[0],vkt2CMS[1],vkt2CMS[2]}
    )
  };
  //std::cout << "kt2 : " << kt2 << std::endl;
  Rivet::FourMomentum pi {
    z * pijt + (1. - z) * y * pkt + kt1 + kt2
  };
  //std::cout << "pi : " << pi << std::endl;
  Rivet::FourMomentum pj {
    (1. - z) * pijt + z * y * pkt - kt1 - kt2
  };
  //std::cout << "pj : " << pj << std::endl;
  Rivet::FourMomentum pk {
    (1. - y) * pkt
  };
  //std::cout << "pk : " << pk << std::endl;
  return Rivet::FourMomenta{pi,pj,pk};
}

Shower::Colours Shower::MakeColours(const std::vector<int>& flavs,
                                    const Colour& colij,
                                    const Colour& colk)
{
  ++_c;
  /// splitter is quark
  if(flavs[0]!=Rivet::PID::GLUON){
    if(flavs[0] > 0){
      return {{_c,0},{colij.first,_c}};
    } else{
      return {{0,_c},{_c,colij.second}};
    }
  } else{ /// splitter is gluon
    if(flavs[1] == Rivet::PID::GLUON){
      return {{colij.first,0},{0,colij.second}};
    } else {
      return {{0,colij.second},{colij.first,0}};
    }
  }
}

void Shower::Run(class EventInfo &evt, const double t) 
{
  _c = 1;
  _tActual = t;
  while( _tActual > _tEnd){
    GeneratePoint(evt);
  }
}

void Shower::GeneratePoint(EventInfo& evt)
{
  while (_tActual > _tEnd)
  {
    double t {_tEnd};
    /// Select Splitter and Spectator
    SelectSplitSpect(evt,t);
    _tActual = t;
    if(t > _tEnd){
      double z {_dipole.selected->get()->GenerateZ(
        1. - _dipole.zp, _dipole.zp)};
      double y {t/_dipole.m2/z/(1.-z)};
      if(y >= 1) continue;
      const double sf {(1. - y) * (*_alphaS)(t) * _dipole.selected->get()->Value(z,y)};
      const double overestimate{_alphaSMax * _dipole.selected->get()->Estimate(z)};
      if ((*_ran)() < sf / overestimate){
        const double phi {2. * M_PI * (*_ran)()};
        Rivet::FourMomenta moms = MakeKinematics(z,y,phi,_dipole.split->GetMomentum(),
                                                 _dipole.spect->GetMomentum());
        Colours cols {MakeColours(_dipole.selected->get()->flavs,
                                  _dipole.split->GetColour(),_dipole.spect->GetColour())};
        _dipole.split->SetColour(cols[0]);
        _dipole.split->SetFlavour(_dipole.selected->get()->flavs[1]);
        _dipole.split->SetMomentum(moms[0]);

        _dipole.spect->SetMomentum(moms[2]);
        evt.Particles.push_back({
          _dipole.selected->get()->flavs[2], moms[1], cols[1]
        });

        return;
      }
    }
  }
  return;
}

void Shower::SelectSplitSpect(EventInfo& evt, double& t)
{
//  std::cout << evt << std::endl;
  for(Partons::iterator split{evt.Particles.begin()+2};
    split!=evt.Particles.end(); ++split){
    for(Partons::iterator spect{evt.Particles.begin()+2};
      spect!=evt.Particles.end();++spect){
      if(spect == split) continue;
      if(not ColourConnected((*split),(*spect))) continue;
      /// select splitting function
      for(auto kern{_kernels.begin()}; kern!=_kernels.end(); ++kern){
        if(kern->get()->flavs[0] != split->GetFlavour()) continue;
        double m2 {(split->GetMomentum() + spect->GetMomentum()).mass2()};
        if(m2 < 4. * _tEnd) continue;
        double zp {0.5 * (1. + sqrt(1. - 4.*_tEnd/m2))};
        double overestimate {_alphaSMax/(2. * M_PI) * kern->get()->Integral(1.-zp,zp)};
        double tt {_tActual * pow((*_ran)(), 1./overestimate)};
        if(tt > t){
          t = tt;
          _dipole.m2       = m2;
          _dipole.zp       = zp;
          _dipole.split    = split;
          _dipole.spect    = spect;
          _dipole.selected = kern;
        }
      }
    }
  }
}

bool Shower::CheckEvent(EventInfo &evt)
{
  Rivet::FourMomentum TotMom {0.,0.,0.,0.};
  for(const auto& p : evt.Particles){
    TotMom += p.GetMomentum();
  }
  std::cout <<" TotMomentum " << TotMom << std::endl;
  if(abs(TotMom.E()) < 1.e-12 and abs(TotMom.px()) < 1.e-12 
    and abs(TotMom.py()) < 1.e-12 and abs(TotMom.pz()) < 1.e-12) return true;
  else return false;
}

bool Shower::ColourConnected(const Particle& pa, const Particle& pb)
{
  return ((pa.GetColour().first > 0 and pa.GetColour().first == pb.GetColour().second) or
          (pa.GetColour().second > 0 and pa.GetColour().first == pb.GetColour().first));
}



/* /// Test main function, TODO: remove it!
int main()
{
  AlphaS alphaS{1,91.1876,0.118, 4.75,1.3};
  Random ran{0};
  Matrix me{91.2,&ran};
  EventInfo evt{me.GeneratePoint()};
  std::cout << evt << std::endl;
  std::cout << Shower::CheckEvent(evt) << std::endl;
  Shower sh{&alphaS,&ran,1.};
  sh.Run(evt,91.2 * 91.2);
  std::cout << evt << std::endl;
  std::cout << sh.CheckEvent(evt) << std::endl;
  return 0;
} */