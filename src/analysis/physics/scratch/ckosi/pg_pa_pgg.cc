#include "pg_pa_pgg.h"
#include "TLorentzVector.h"
#include "analysis/utils/fitter/KinFitter.h"
#include "utils/uncertainties/Interpolated.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_kosi_pg_pa_pgg::scratch_kosi_pg_pa_pgg(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                utils::UncertaintyModels::Interpolated::Type_t::MC,
                make_shared<utils::UncertaintyModels::FitterSergey>())),
    fitter(nullptr, opts->Get<bool>("FitZVertex", true)),
    npfitter(nullptr, opts->Get<bool>("FitZVertex", true))
    //promptrandom(ExpConfig::Setup::Get())
{
    const BinSettings tagger_time_bins(200, -50, 50);
    const BinSettings phi_bins(360, 0, 360);
    const BinSettings numcands_time_bins(10, 0, 10);
    const BinSettings theta_eta_bins(100, -1, 1);
    const BinSettings gamma_im_bins(300, 0, 1000);
    const BinSettings gamma_mm_bins(500, 0, 1500);
    const BinSettings energy_bins(88,1.07352,1.51254);

    promptrandom.AddPromptRange({-8,8});
    promptrandom.AddRandomRange({-50,-12});
    promptrandom.AddRandomRange({12,50});

    fitter.SetZVertexSigma(3.0);
    npfitter.SetZVertexSigma(3.0);

        h_TaggerTime = HistFac.makeTH1D("Tagger Time",     // title
                                        "t [ns]", "#",     // xlabel, ylabel
                                        tagger_time_bins,  // our binnings
                                        "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                        );

        h_NumCands = HistFac.makeTH1D("Number of candidates",
                                       "Candidates", "#",
                                       numcands_time_bins,
                                       "h_NumCands");

        h_Theta = HistFac.makeTH1D("Theta of Meson",
                                   "cos_Theta [deg]", "#",
                                   theta_eta_bins,
                                   "h_Theta");

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
                                theta_eta_bins, energy_bins,    // our binnings
                                "h_cs",false);

        h_Theta_nokin = HistFac.makeTH1D("Theta of particle (no kin)",
                                      "cos_Theta [deg]", "#",
                                      theta_eta_bins,
                                      "h_Theta_nokin");

        h_GammaIM_nocut = HistFac.makeTH1D("Invariant mass of decay gammas (no cuts)",
                                      "Mass [MeV]", "#",
                                      gamma_im_bins,
                                      "h_GammaIM_nocut");

        h_GammaIM_cut1 = HistFac.makeTH1D("Invariant mass of decay gammas (cut on theta)",
                                      "Mass [MeV]", "#",
                                      gamma_im_bins,
                                      "h_GammaIM_cut1");

        h_GammaIM_cut2 = HistFac.makeTH1D("Invariant mass of decay gammas (cut on theta/missing mass/)",
                                      "Mass [MeV]", "#",
                                      gamma_im_bins,
                                      "h_GammaIM_cut2");

        h_GammaIM_np = HistFac.makeTH1D("Invariant mass of decay gammas (cut on theta/missing mass/np fit)",
                                      "Mass [MeV]", "#",
                                      gamma_im_bins,
                                      "h_GammaIM_np");

        h_GammaIM_cut3 = HistFac.makeTH1D("Invariant mass of decay gammas (cut on theta/missing mass/coplanarity)",
                                      "Mass [MeV]", "#",
                                      gamma_im_bins,
                                      "h_GammaIM_cut3");

        h_MM_nokin = HistFac.makeTH1D("Missing mass (no kin)",
                                      "Mass [MeV]", "#",
                                      gamma_mm_bins,
                                      "h_GammaMM_nokin");

        h_coplanar_nokin = HistFac.makeTH1D("Coplanarity (no kin)",
                                      "Phi [deg]", "#",
                                      phi_bins,
                                      "h_coplanar_nokin");
}

void scratch_kosi_pg_pa_pgg::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed().Candidates;
    TrigSim.ProcessEvent(event);

    fitter.SetUncertaintyModel(fit_model);
    npfitter.SetUncertaintyModel(fit_model);

    //Candidates
    TCandidatePtrList neutral;
    TCandidatePtrList charged;

    //Kinfit particles
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


    //Check if decay particles are 2 neutrals and 1 charged
    if(neutral.size() == 2){

        //Loop over taggerhits
        for (const auto& taggerhit : event.Reconstructed().TaggerHits) {
            h_TaggerTime->Fill(taggerhit.Time);

            double correctedTime = TrigSim.GetCorrectedTaggerTime(taggerhit);
            promptrandom.SetTaggerTime(correctedTime);
            if(promptrandom.State() == PromptRandom::Case::Outside) continue;
            const double weight = promptrandom.FillWeight();

            //Create Candidate            
            TParticle Meson_nokin;

            //Before kinfit
            for(const auto& photon : neutral){
                Meson_nokin += TParticle(ParticleTypeDatabase::Photon, photon);
            }
            auto initial_nokin = taggerhit.GetPhotonBeam() + LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
            auto missingMass_nokin = initial_nokin - Meson_nokin;

            h_GammaIM_nocut->Fill(Meson_nokin.M(),weight);

            h_Theta_nokin->Fill(cos(Meson_nokin.Theta()),weight);
            if (cos(Meson_nokin.Theta()) < 0.8) continue;
            h_GammaIM_cut1->Fill(Meson_nokin.M(),weight);

            h_MM_nokin->Fill(missingMass_nokin.M(),weight);
            if (missingMass_nokin.M() > 990) continue;
            if (missingMass_nokin.M() < 880) continue;
            h_GammaIM_cut2->Fill(Meson_nokin.M(),weight);


            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            //h_MM_nokin->Fill(missingMass_nokin.M(),weight);
            //h_ThetaEta_nokin->Fill(cos(Meson_nokin.Theta()),weight);
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            if (charged.size() == 1) {
                proton = protons.at(0);

                auto coplanarity_nokin  = abs(Meson_nokin.Phi() - charged.at(0)->Phi)*radtodeg;

                h_coplanar_nokin->Fill(coplanarity_nokin,weight);
                if (coplanarity_nokin > 190) continue;
                if (coplanarity_nokin < 170) continue;
                h_GammaIM_cut3->Fill(Meson_nokin.M(),weight);

                //Do kinfit
                double best_prob = std_ext::NaN;

                APLCON::Result_t result = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);

                if(result.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(best_prob, result.Probability))
                    continue;

                photons = fitter.GetFittedPhotons();

                auto Meson = *photons.at(0) + *photons.at(1);
                h_GammaIM->Fill(Meson.M(),weight);

//                h_MM->Fill(missingMass.M(),weight);
//                h_ThetaEta->Fill(cos(Meson.Theta()),weight);
//                h_cs->Fill(cos(Meson.Theta()),initial.M()/1000,weight);
//                h_TaggerTime->Fill(correctedTime,weight);

            } else if (charged.size() == 0) {

                vector<LorentzVec> fitphotons;

                //Do kinfit
                double best_prob = std_ext::NaN;

                APLCON::Result_t result = npfitter.DoFit(taggerhit.PhotonEnergy, photons);

                if(result.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(best_prob, result.Probability))
                    continue;

                for(const auto& fph : npfitter.GetFittedPhotons())
                    fitphotons.push_back(*fph);

                LorentzVec Meson;
                for(const auto& ph: fitphotons) Meson += ph;

                h_GammaIM_np->Fill(Meson.M(),weight);
                h_GammaIM->Fill(Meson.M(),weight);

            }
        }
    }
}


void scratch_kosi_pg_pa_pgg::ShowResult()
{
    ant::canvas(GetName()+": Basic plots")
            << drawoption("hist")
            << h_NumCands
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Particle Dist")
            << drawoption("hist")
            << h_Theta
            << h_MM
            << endc;

    ant::canvas(GetName()+": Particle Dist (no kin)")
            << drawoption("hist")               
            << h_TaggerTime
            << h_Theta_nokin
            << h_MM_nokin
            << h_coplanar_nokin
            << endc;

    ant::canvas(GetName()+": cross section")
            << drawoption("Surf")
            << h_cs
            << endc;

    ant::canvas(GetName()+": Invariant masses")
            << drawoption("hist")
            << h_GammaIM_nocut
            << h_GammaIM_cut1
            << h_GammaIM_cut2
            << h_GammaIM_cut3
            << h_GammaIM
            << h_GammaIM_np
            << endc;

}

void scratch_kosi_pg_pa_pgg::Finish()
{
    cout<<"Finished"<<endl;
}

AUTO_REGISTER_PHYSICS(scratch_kosi_pg_pa_pgg)
