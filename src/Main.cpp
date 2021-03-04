#include "Shower.hpp"

#include "Kernels.hpp"
#include "Matrix.hpp"
#include "QCD.hpp"
#include "Random.hpp"

#include "Rivet/Rivet.hh"
#include "Rivet/AnalysisHandler.hh"

bool ToHepMCEvent(const EventInfo &evt, HepMC::GenEvent& hepevt)
{
  hepevt.use_units(HepMC::Units::GEV, HepMC::Units::MM);
  hepevt.set_event_number(evt.EvtNumber);
  HepMC::WeightContainer weights{};
  weights["Nominal"] = evt.dxs;
  weights["MEWeight"] = evt.lome;
  hepevt.weights() = weights;
  HepMC::GenVertex * vertex {new HepMC::GenVertex{}};
  std::vector<HepMC::GenParticle* > inparticles;

  for(auto p{evt.Particles.begin()}; p != evt.Particles.begin() + 2; ++p){
    HepMC::FourVector mom {p->GetMomentum().px(), p->GetMomentum().py(),
                           p->GetMomentum().pz(),p->GetMomentum().E()};
    int status{(p->GetColour().first + p->GetColour().second) == 0 ? 4 : 11};
    HepMC::GenParticle *part{
        new HepMC::GenParticle{mom, p->GetFlavour(), status}};
    vertex->add_particle_in(part);
    inparticles.push_back(part);
  }

  for (auto p{evt.Particles.begin()+2}; p != evt.Particles.end(); ++p) {
    HepMC::FourVector mom{p->GetMomentum().px(), p->GetMomentum().py(),
                          p->GetMomentum().pz(), p->GetMomentum().E()};
    HepMC::GenParticle *part{new HepMC::GenParticle{mom, p->GetFlavour(), 1}};
    vertex->add_particle_out(part);
  }
  hepevt.add_vertex(vertex);
  return true;
}

int main()
{
  AlphaS alphaS{1, 91.1876, 0.118, 4.75, 1.3};
  Random ran{0};
  Matrix me{91.2, &ran};
  Shower shower{&alphaS, &ran, 1.};

  Rivet::AnalysisHandler rivet;
  rivet.setIgnoreBeams(true);
  rivet.addAnalyses({"ALEPH_2004_S5765862", "LL_JetRates"});

  long int TotEvents{100000};
  double totalxs{0.0}, sumW{0.0}, sumW2{0.0}, err{0.0};

  for (size_t i{0}; i < TotEvents; ++i){
    EventInfo evt{me.GeneratePoint()};
    const double t {(evt.Particles[0].GetMomentum() + evt.Particles[1].GetMomentum()).mass2()};
    shower.Run(evt, t);
    evt.EvtNumber = i;
    HepMC::GenEvent hepevt;
    ToHepMCEvent(evt,hepevt);
    sumW  += evt.dxs;
    sumW2 += (evt.dxs * evt.dxs);
    totalxs = sumW / TotEvents;
    err = sqrt(abs(sumW2 / TotEvents - totalxs * totalxs) / (TotEvents -1));
    HepMC::GenCrossSection xs;
    xs.set_cross_section(totalxs,err);
    hepevt.set_cross_section(xs);
    rivet.analyze(hepevt);
    if(i % 1000 == 0)
      std::cout << "\rEvent " << i <<  ", \u03c3 = "  << totalxs << " \u00B1 "
                << err << " [pb] (" << 100. * err/totalxs << " %)" << std::flush;
  }
  totalxs = sumW/TotEvents;
  err     = sqrt(abs(sumW2 / TotEvents - totalxs * totalxs) / (TotEvents -1));

  const double wgtfract {rivet.sumW()/rivet.sumW()};
  const double rivetxs {rivet.nominalCrossSection()};
  const double thisxs {rivetxs * wgtfract};
  rivet.setCrossSection(thisxs,0.0,true);
  rivet.finalize();
  rivet.writeData("result.yoda");
  std::cout << std::endl;
  std::cout << "=============================================\n";
  std::cout <<  "  \u03c3 = "  << totalxs << " \u00B1 "
            << err << " [pb] (" << 100. * err/totalxs << " %)" << std::endl;
  std::cout << "=============================================\n";
  std::cout << "Completed!" << std::endl;
  return 0;
}