#include "pg_peta_pgg.h"
#include "TLorentzVector.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_kosi_pg_peta_pgg::scratch_kosi_pg_peta_pgg(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    const BinSettings tagger_time_bins(2000, -200, 200);
    const BinSettings theta_time_bins(100, -1, 1);
    const BinSettings numcands_time_bins(10, 0, 10);
    const BinSettings theta_eta_bins(100, -1, 1);
    const BinSettings gamma_im_bins(1000, 0, 1000);
    const BinSettings gamma_mm_bins(1500, 0, 1500);
    const BinSettings energy_bins(88,1.07352,1.51254);

        h_TaggerTime = HistFac.makeTH1D("Tagger Time",     // title
                                        "t [ns]", "#",     // xlabel, ylabel
                                        tagger_time_bins,  // our binnings
                                        "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                        );

        h_NumCands = HistFac.makeTH1D("Number of candidates",
                                       "Candidates", "#",
                                       numcands_time_bins,
                                       "h_NumCands");

        h_ThetaEta = HistFac.makeTH1D("Theta of Eta",
                                      "cos_Theta [deg]", "#",
                                      theta_eta_bins,
                                      "h_ThetaEta");

        h_ThetaEtaCut = HistFac.makeTH1D("Theta of Eta after cut",
                                      "cos_Theta [deg]", "#",
                                      theta_eta_bins,
                                      "h_ThetaEtaCut");

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
}

void scratch_kosi_pg_peta_pgg::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed().Candidates;
    TrigSim.ProcessEvent(event);
    const double weight = promptrandom.FillWeight();

    TCandidatePtrList neutral;
    TCandidatePtrList charged;

    //Loop over candidates and distinguish between charged and neutral

    for (const auto& c : cands.get_iter()) {
        if (c->VetoEnergy <= vetoEthreshold) neutral.emplace_back(c);
        else charged.emplace_back(c);
        }

    h_NumCands->Fill(cands.size());


    /*Check if decay particles are 2 gammas
     * Eta->gg
     * Eta->Pi+Pi-Pi0->Pi+Pi-gg
     * Eta->3Pi0->2Pi0gg
     */
    if(neutral.size() == 2 && (charged.size() == 0 || charged.size() == 1)){

        //Loop over taggerhits
        for (const auto& taggerhit : event.Reconstructed().TaggerHits) {

            //Create Candidate
            TParticle Meson;

            for(const auto& photon : neutral){
                Meson += TParticle(ParticleTypeDatabase::Photon, photon);
            }

            //Initial conditions
            auto initial = taggerhit.GetPhotonBeam() + LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
            auto missingMass = initial - Meson;

            //Check if missing mass is in Proton range
            if(     missingMass.M() < ParticleTypeDatabase::Proton.Mass() + 40
                 && missingMass.M() > ParticleTypeDatabase::Proton.Mass() - 40){

                h_GammaIM->Fill(Meson.M(),weight);
                h_MM->Fill(missingMass.M(),weight);
                h_ThetaEta->Fill(cos(Meson.Theta()),weight);
                h_cs->Fill(cos(Meson.Theta()),initial.M()/1000,weight);
                h_TaggerTime->Fill(taggerhit.Time,weight);

                //Cut away everything thats Cos(°) < 0.9
                if(cos(Meson.Theta()) >= 0.9){

                  h_ThetaEtaCut->Fill(cos(Meson.Theta()),weight);
                }
            }
        }
    }
}


void scratch_kosi_pg_peta_pgg::ShowResult()
{
    ant::canvas(GetName()+": Basic plots")
            << h_TaggerTime
            << h_NumCands
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Eta Dist")
            << h_ThetaEta
            << h_ThetaEtaCut
            << h_GammaIM
            << h_MM
            << endc;

    ant::canvas(GetName()+": cross section")
            << drawoption("Surf")
            << h_cs
            << endc;

}

void scratch_kosi_pg_peta_pgg::Finish()
{
    cout<<"Finished"<<endl;
}

AUTO_REGISTER_PHYSICS(scratch_kosi_pg_peta_pgg)
