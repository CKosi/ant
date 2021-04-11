#include "pg_ppi0_pgg.h"
#include "TLorentzVector.h"
#include "analysis/utils/fitter/KinFitter.h"
#include "utils/uncertainties/Interpolated.h"
#include "analysis/utils/uncertainties/FitterSergey.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_kosi_pg_ppi0_pgg::scratch_kosi_pg_ppi0_pgg(const string& name, OptionsPtr opts) :
    Physics(name, opts),
  fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
              utils::UncertaintyModels::Interpolated::Type_t::MC,
              make_shared<utils::UncertaintyModels::FitterSergey>())),
  fitter(nullptr, opts->Get<bool>("FitZVertex", true))
{
    const BinSettings tagger_time_bins(2000, -200, 200);
    const BinSettings theta_time_bins(100, -1, 1);
    const BinSettings numcands_time_bins(10, 0, 10);
    const BinSettings theta_pi0_bins(100, -1, 1);
    const BinSettings gamma_im_bins(1000, 0, 1000);
    const BinSettings gamma_mm_bins(1500, 0, 1500);
    const BinSettings energy_bins(88,1.07352,1.51254);

    fitter.SetZVertexSigma(3.0);

        h_TaggerTime = HistFac.makeTH1D("Tagger Time",     // title
                                        "t [ns]", "#",     // xlabel, ylabel
                                        tagger_time_bins,  // our binnings
                                        "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                        );

        h_NumCands = HistFac.makeTH1D("Number of candidates",
                                       "Candidates", "#",
                                       numcands_time_bins,
                                       "h_NumCands");

        h_ThetaPi0 = HistFac.makeTH1D("Theta of Pi0",
                                      "cos_Theta [deg]", "#",
                                      theta_pi0_bins,
                                      "h_ThetaPi0");

        h_ThetaPi0Cut = HistFac.makeTH1D("Theta of Pi0 after cut",
                                      "cos_Theta [deg]", "#",
                                      theta_pi0_bins,
                                      "h_ThetaPi0Cut");

        h_GammaIM = HistFac.makeTH1D("Invariant mass of decay gammas",
                                      "Mass [MeV]", "#",
                                      gamma_im_bins,
                                      "h_GammaIM");

        h_MM = HistFac.makeTH1D("Missing mass",
                                      "Mass [MeV]", "#",
                                      gamma_mm_bins,
                                      "h_GammaMM");

        h_cs = HistFac.makeTH2D("Cross section",
                                "cos_Theta [deg]","sqrt_S [GeV]",
                                theta_pi0_bins, energy_bins,    // our binnings
                                "h_cs",false);
}

void scratch_kosi_pg_ppi0_pgg::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed().Candidates;
    TrigSim.ProcessEvent(event);
    const double weight = promptrandom.FillWeight();

    fitter.SetUncertaintyModel(fit_model);

    TCandidatePtrList neutral;
    TCandidatePtrList charged;

    TParticleList protons;
    TParticleList photons;
    TParticlePtr proton;

    //Loop over candidates and distinguish between charged and neutral

    for (const auto& c : cands.get_iter()) {
        if (c->VetoEnergy <= vetoEthreshold){ neutral.emplace_back(c);
         photons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon, c));
        }
        else charged.emplace_back(c);
        protons.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Proton, c));
        }


    h_NumCands->Fill(cands.size());


    //Check if decay particles are 2 gammas Pi0->gg
    if(neutral.size() == 2){

        //Create Pi0
        TParticle Pi0;


        //Loop over taggerhits
        for (const auto& taggerhit : event.Reconstructed().TaggerHits) {

            for(const auto& photon : neutral){
                Pi0 += TParticle(ParticleTypeDatabase::Photon, photon);
            }

            h_TaggerTime->Fill(taggerhit.Time,weight);

            //Initial conditions
            auto initial = taggerhit.GetPhotonBeam() + LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
            auto missingMass = initial - Pi0;

//            if (missingMass.M() > 990) continue;
//            if (missingMass.M() < 880) continue;

            h_GammaIM->Fill(Pi0.M(), weight);
            h_MM->Fill(missingMass.M(),weight);
            h_ThetaPi0->Fill(cos(Pi0.Theta()),weight);
            h_cs->Fill(cos(Pi0.Theta()),initial.M()/1000,weight);

            //Cut away everything thats Cos(°) < 0.9

            if(cos(Pi0.Theta()) >= 0.9){

               h_ThetaPi0Cut->Fill(cos(Pi0.Theta()),weight);
            }
        }
    }
}


void scratch_kosi_pg_ppi0_pgg::ShowResult()
{
    ant::canvas(GetName()+": Basic plots")
            << h_TaggerTime
            << h_NumCands
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Pi0 Dist")
            << h_ThetaPi0
            << h_ThetaPi0Cut
            << h_GammaIM
            << h_MM
            << endc;

    ant::canvas(GetName()+": cross section")
            << drawoption("Surf")
            << h_cs
            << endc;

}

void scratch_kosi_pg_ppi0_pgg::Finish()
{
    cout<<"Finished"<<endl;
}

AUTO_REGISTER_PHYSICS(scratch_kosi_pg_ppi0_pgg)
