#include "pg_ppi0_pgg.h"
#include "TLorentzVector.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_kosi_pg_ppi0g_pgg::scratch_kosi_pg_ppi0g_pgg(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    const BinSettings tagger_time_bins(2000, -200, 200);
    const BinSettings theta_time_bins(180, 0, 180);
    const BinSettings numcands_time_bins(10, 0, 10);
    const BinSettings theta_pi0_bins(180, 0, 180);
    const BinSettings theta_pi0Zoom_bins(25, 0, 5);

        h_TaggerTime = HistFac.makeTH1D("Tagger Time",     // title
                                        "t [ns]", "#",     // xlabel, ylabel
                                        tagger_time_bins,  // our binnings
                                        "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                        );

        h_ThetaDist = HistFac.makeTH1D("Theta distribution",
                                       "Angle [°]", "#",
                                       theta_time_bins,
                                       "h_ThetaDist");

        h_NumCands = HistFac.makeTH1D("Number of candidates",
                                       "Candidates", "#",
                                       numcands_time_bins,
                                       "h_NumCands");

        h_ThetaPi0 = HistFac.makeTH1D("Theta of Pi0",
                                      "Angle [°]", "#",
                                      theta_pi0_bins,
                                      "h_ThetaPi0");

        h_ThetaPi0Zoom = HistFac.makeTH1D("Theta of Pi0Zommed",
                                      "Angle [°]", "#",
                                      theta_pi0Zoom_bins,
                                      "h_ThetaPi0Zoom");
}

void scratch_kosi_pg_ppi0g_pgg::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& cands = event.Reconstructed().Candidates;
    TrigSim.ProcessEvent(event);

    TCandidatePtrList gammas;
    TCandidatePtrList protons;

    vector<double> protonCalE;
    vector<double> gammaCalE;

    TParticle pGammas[num_gammas];
    TParticle pProton[num_protons];

    //Loop over candidates

    for (const auto& c : cands.get_iter()) {
        if (c->VetoEnergy <= vetoEthreshold){
            gammas.emplace_back(c);
            gammaCalE.push_back(c->CaloEnergy);
        }else{
            protons.emplace_back(c);
            protonCalE.push_back(c->CaloEnergy);
        }

        h_ThetaDist->Fill(c->Theta*radtodeg);        //Thetas of all candidates
        }

    if(cands.size() == 3){num_3decay += 1;}
    h_NumCands->Fill(cands.size());                 //Distribution of candidate numbers


    //Loop over taggerhits
    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {

        //Check if decay particles are 2 gammas and 1 proton
        if(gammaCalE.size() == num_gammas && protonCalE.size() == num_protons){
            for(int i = 0; i<num_gammas; i++){
                pGammas[i] = TParticle(ParticleTypeDatabase::Photon, gammas[i]);
            }
            for(int i = 0; i<num_protons; i++){
                pProton[i] = TParticle(ParticleTypeDatabase::Proton, protons[0]);
            }

            //Get vector of pi0
            TLorentzVector pi0_decay = pGammas[0] + pGammas[1];

            TParticle pPi0 = TParticle(ParticleTypeDatabase::Pi0,(TLorentzVector)(pi0_decay));
            //Fill with theta of pi0
            h_ThetaPi0->Fill(pPi0.Theta()*radtodeg);
        }

        h_TaggerTime->Fill(taggerhit.Time);                         //Tagger times
        }

}


void scratch_kosi_pg_ppi0g_pgg::ShowResult()
{
    ant::canvas(GetName()+": Basic plots")
            << h_TaggerTime
            << h_ThetaDist
            << h_NumCands
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Pi0 Dist")
            << h_ThetaPi0
            << endc; // actually draws the canvas

}

void scratch_kosi_pg_ppi0g_pgg::Finish()
{
    cout<<"Finished"<<endl;
    cout<<"Number of events with 3 decay products "<<num_3decay<<endl;
}

AUTO_REGISTER_PHYSICS(scratch_kosi_pg_ppi0g_pgg)
