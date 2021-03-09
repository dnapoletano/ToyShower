#include "Matrix.hpp"

myMatrix::myMatrix(const double& ecms, Random* random)
: _ecms{ecms}, _ewparams{}, ran{random}
{}

double myMatrix::ME2(const int& flav, const double& s, const double& t) const
{
  const double ve {_ewparams.ae - 2. * _ewparams.qe * _ewparams.sin2tw};
  double qf{}, af{};
  const bool IsUpQuark {
    (abs(flav) == Rivet::PID::UQUARK) or
    (abs(flav) == Rivet::PID::CQUARK)
  };
  if (IsUpQuark) {
    qf = 2./3.;
    af = 0.5;
  } else{
    qf = -1./3.;
    af = -0.5;
  }
  const double vf {af - 2. * qf * _ewparams.sin2tw};
  const double kappa {1./(4.*_ewparams.sin2tw*(1.-_ewparams.sin2tw))};
  const double chi1 {kappa * s * (s - _ewparams.mz2)/((s-_ewparams.mz2)*(s-_ewparams.mz2) + _ewparams.gz2 * _ewparams.mz2)};
  const double chi2 {(kappa * s * kappa * s)/((s-_ewparams.mz2)*(s-_ewparams.mz2) +
    _ewparams.gz2 * _ewparams.mz2)};
  const double term1 {
    (1. + pow(1. + 2.*t/s, 2)) * (pow(qf * _ewparams.qe, 2) +
      2. * (qf * _ewparams.qe * vf * ve)*chi1 +
      (_ewparams.ae * _ewparams.ae + ve * ve) * (af * af + vf * vf) * chi2)
  };
  const double term2 {
    (1. + 2.*t/s) * (4. * _ewparams.qe * qf * _ewparams.ae * af * chi1 +
      8.0 * _ewparams.ae * ve * af * vf * chi2)
  };
  return pow(4.*M_PI*_ewparams.alpha0,2) * 3.0 * (term1 + term2);
}

EventInfo myMatrix::GeneratePoint()
{
  EventInfo evtinfo{};
  evtinfo.Particles.clear();
  const double ct  {2. * (*ran)() - 1.};
  const double st  {sqrt(1. - ct * ct)};
  const double phi {2.* M_PI * (*ran)()};

  const Rivet::FourMomentum pa{_ecms/2.,0.,0.,_ecms/2.};
  const Rivet::FourMomentum pb{_ecms/2.,0.,0.,-_ecms/2.};
  const Rivet::FourMomentum p1{_ecms/2.,
                               _ecms/2. * st * cos(phi),
                               _ecms/2. * st * sin(phi),
                               _ecms/2. * ct};
  const Rivet::FourMomentum p2{_ecms/2.,
                               -_ecms/2. * st * cos(phi),
                               -_ecms/2. * st * sin(phi),
                               -_ecms/2. * ct};

  evtinfo.Particles.push_back({Rivet::PID::POSITRON, -pa});
  evtinfo.Particles.push_back({Rivet::PID::ELECTRON, -pb});
  const int fl {ran->randint()};
  evtinfo.Particles.push_back({fl, p1 ,std::make_pair<int,int>(1,0)});
  evtinfo.Particles.push_back({ -fl, p2,std::make_pair<int,int>(0,1)});

  const double lome {ME2(fl, (pa + pb).mass2(), (pa - p1).mass2())};
  const double dxs {5. * lome * 3.89379656e8 / 8. / M_PI / 2. / _ecms /_ecms};
  evtinfo.dxs = dxs;
  evtinfo.lome = lome;
  return evtinfo;
}