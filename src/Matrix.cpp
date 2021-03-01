#include "Matrix.hpp"

Matrix::Matrix(const double& ecms)
: _ecms{ecms}, _ewparams{}, re{0}, urng{0.,1.}
{}

double Matrix::ME2(const int& flav, const double& s, const double& t) const
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
  const double chi2 {kappa * kappa/((s-_ewparams.mz2)*(s-_ewparams.mz2) +
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

EventInfo Matrix::GeneratePoint()
{
  EventInfo evtinfo{};
  evtinfo.Particles.clear();
  const double ct  {2. * urng(re) - 1.};
  const double st  {sqrt(1. - ct * ct)};
  const double phi {2.* M_PI * urng(re)};

  const Rivet::FourVector pa{_ecms/2.,0.,0.,_ecms/2.};
  const Rivet::FourVector pb{_ecms/2.,0.,0.,-_ecms/2.};
  const Rivet::FourVector p1{_ecms/2.,
                             _ecms/2. * st * cos(phi),
                             _ecms/2. * st * sin(phi),
                             _ecms/2. * ct};
  const Rivet::FourVector p2{_ecms/2.,
                             - _ecms/2. * st * cos(phi),
                             - _ecms/2. * st * sin(phi),
                             - _ecms/2. * ct};

  evtinfo.Particles.push_back({
    Rivet::PID::ELECTRON, -pa});
  evtinfo.Particles.push_back({
    Rivet::PID::POSITRON, -pb});
  const int fl { static_cast<int>(( - 5 * 11 * urng(re)))};
  evtinfo.Particles.push_back({ fl, p1 });
  evtinfo.Particles.push_back({ - fl, p2});

  const double lome {ME2(fl, (pa + pb).mod2(), (pa - p1).mod2())};
  const double dxs {5. * lome * 3.89379656e8 / 8. / M_PI / 2. / _ecms /_ecms};
  evtinfo.dxs = dxs;
  return evtinfo;
}