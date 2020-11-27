#include "MC_ppi0_pgg.h"
#include "base/WrapTFile.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "base/ParticleTypeTree.h"
#include "utils/uncertainties/FitterSergey.h"
#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "expconfig/detectors/TAPS.h"

#include <iostream>
#include <memory>

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_damaurer_MC_ppi0_pgg::scratch_damaurer_MC_ppi0_pgg(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    taps = make_shared<expconfig::detector::TAPS_2013_11>(false, false, false);

    const BinSettings tagger_time_bins(1000, -200, 200);
    const BinSettings Veto_Energy_bins(1000, 0, 40);
    const BinSettings theta_bins(180, 0, 180);
    const BinSettings theta_bins_CB(140, 19, 160);
    const BinSettings theta_bins_TAPS(25, 0, 25);
    const BinSettings phi_bins(100, -180, 180);
    const BinSettings E_Photon_bins(100, 0, 1000);
    const BinSettings E_Proton_bins(100, 900, 1600);
    const BinSettings Im_photon_bins(1000, 0, 300);
    const BinSettings Im_proton_bins(1000, 800, 1600);
    const BinSettings initialBeamE_bins(1000,0,0.9);
    const BinSettings statistic_bins(number_of_bins,lower_edge,upper_edge);
    const BinSettings cos_bins(100,-1,1);
    const BinSettings sqrt_S_bins(100,s_square_min,s_square_max);

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    auto hfTagger = new HistogramFactory("Tagger", HistFac, "");
    auto hfTheta = new HistogramFactory("Polar angles", HistFac, "");
    auto hfPhi = new HistogramFactory("Azimuth angles", HistFac, "");
    auto hfIm = new HistogramFactory("Invariant mass", HistFac, "");
    auto hfEnergy = new HistogramFactory("Particle energies", HistFac, "");
    auto hfStatistics = new HistogramFactory("Statistics hist", HistFac, "");
    auto hfCrossCheck = new HistogramFactory("Cross section check", HistFac, "");

    h_TaggerTime = hfTagger->makeTH1D("Tagger Time",     // title
                                    "t [ns]", "#",     // xlabel, ylabel
                                    tagger_time_bins,  // our binnings
                                    "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                    );
    
    h_InitialBeamE = hfTagger->makeTH1D("Photon beam - Energy distribution",     // title
                                   "E [GeV]", "#",     // xlabel, ylabel
                                   initialBeamE_bins,  // our binnings
                                   "h_InitialBeamE"     // ROOT object name, auto-generated if omitted
                                   );

    h_nClusters = hfTagger->makeTH1D("Number of Clusters", "nClusters", "#", BinSettings(15), "h_nClusters");

    h_VetoEnergies = hfTagger->makeTH1D("Veto Energies", "E [MeV]", "#", Veto_Energy_bins, "h_VetoEnergies");


    h_PolarAngles = hfTheta->makeTH1D("Polarangles",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins,  // our binnings
                                       "h_PolarAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_ProtonPolarAngles = hfTheta->makeTH1D("Proton Polarangles",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins,  // our binnings
                                       "h_ProtonPolarAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_PhotonPolarAngles = hfTheta->makeTH1D("Photon Polarangles",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins,  // our binnings
                                       "h_PhotonPolarAngles", false    // ROOT object name, auto-generated if omitted
                                       );

    h_PhotonPolarAngles_CB = hfTheta->makeTH1D("Photon Polarangles CB",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins_CB,  // our binnings
                                       "h_PhotonPolarAngles_CB", false    // ROOT object name, auto-generated if omitted
                                       );

    h_PhotonPolarAngles_TAPS = hfTheta->makeTH1D("Photon Polarangles TAPS",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins_TAPS,  // our binnings
                                       "h_PhotonPolarAngles_TAPS", false    // ROOT object name, auto-generated if omitted
                                       );

    h_AzimuthAngles = hfPhi->makeTH1D("AzimuthAngles",     // title
                                       "#Phi [deg]", "#",     // xlabel, ylabel
                                       phi_bins,  // our binnings
                                       "h_AzimuthAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_PhotonAzimuthAngles = hfPhi->makeTH1D("Photon AzimuthAngles",     // title
                                       "#Phi [deg]", "#",     // xlabel, ylabel
                                       phi_bins,  // our binnings
                                       "h_PhotonAzimuthAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_ProtonAzimuthAngles = hfPhi->makeTH1D("Proton AzimuthAngles",     // title
                                       "#Phi [deg]", "#",     // xlabel, ylabel
                                       phi_bins,  // our binnings
                                       "h_ProtonAzimuthAngles", false     // ROOT object name, auto-generated if omitted
                                       );

    h_PhotonIm = hfIm->makeTH1D("Photons invariant mass",     // title
                                  "M [MeV]", "#",     // xlabel, ylabel
                                  Im_photon_bins,  // our binnings
                                  "h_PhotonIm", false     // ROOT object name, auto-generated if omitted
                                  );

    h_ProtonIm = hfIm->makeTH1D("Protons invariant mass",     // title
                                  "M [MeV]", "#",     // xlabel, ylabel
                                  Im_proton_bins,  // our binnings
                                  "h_ProtonIm", false     // ROOT object name, auto-generated if omitted
                                  );

    h_ProtonPeak = hfIm->makeTH1D("Missing particle reconstructed events around proton mass",     // title
                                  "M [MeV]", "#",     // xlabel, ylabel
                                  Im_proton_bins,  // our binnings
                                  "h_ProtonImAbove80", false     // ROOT object name, auto-generated if omitted
                                  );


    h_PhotonETheta = hfEnergy->makeTH2D("Photon energy vs theta", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins, E_Photon_bins,    // our binnings
                                       "h_PhotonETheta", true    // ROOT object name, auto-generated if omitted
                                       );

    h_PhotonETheta_CB = hfEnergy->makeTH2D("Photon energy vs theta", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins_CB, E_Photon_bins,    // our binnings
                                       "h_PhotonETheta_CB", true    // ROOT object name, auto-generated if omitted
                                       );

    h_PhotonETheta_TAPS = hfEnergy->makeTH2D("Photon energy vs theta", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins_TAPS, E_Photon_bins,    // our binnings
                                       "h_PhotonETheta_TAPS", true    // ROOT object name, auto-generated if omitted
                                       );

    h_ProtonETheta = hfEnergy->makeTH2D("Proton energy vs theta", //title
                                        "Theta [deg]","E [MeV]", // xlabel, ylabel
                                        theta_bins, E_Proton_bins,    // our binnings
                                        "h_ProtonETheta", true    // ROOT object name, auto-generated if omitted
                                        );

    h_PhotonEPhi = hfEnergy->makeTH2D("Photon energy vs phi", //title
                                       "Phi [deg]","E [MeV]",  // xlabel, ylabel
                                       phi_bins, E_Photon_bins,    // our binnings
                                       "h_PhotonEPhi", true    // ROOT object name, auto-generated if omitted
                                       );

    h_ProtonEPhi = hfEnergy->makeTH2D("Proton energy vs phi", //title
                                        "Phi [deg]","E [MeV]", // xlabel, ylabel
                                        phi_bins, E_Proton_bins,    // our binnings
                                        "h_ProtonEPhi", true    // ROOT object name, auto-generated if omitted
                                        );

    h_Reconstructed_Data_Statistics = hfStatistics->makeTH1D("Statistics Hist",     // title
                                     "Cases", "Reconstruction efficiency",     // xlabel, ylabel
                                     statistic_bins,  // our binnings
                                     "h_Reconstructed_Data_Statistics", true    // ROOT object name, auto-generated if omitted
                                     );

    h_doubly_ap_DCS_reconstructed_lab = hfCrossCheck->makeTH2D("Cross section check in lab-frame", //title
                                                                        "cos_Theta [deg]","sqrt_S [GeV]", // xlabel, ylabel
                                                                        cos_bins, sqrt_S_bins,    // our binnings
                                                                        "h_doubly_ap_DCS_reconstructed_lab", true    // ROOT object name, auto-generated if omitted
                                                                        );
    
    h_doubly_ap_DCS_reconstructed_cmFrame = hfCrossCheck->makeTH2D("Cross section check in cm-frame", //title
                                                                        "cos_Theta* [deg]","sqrt_S [GeV]", // xlabel, ylabel
                                                                        cos_bins, sqrt_S_bins,    // our binnings
                                                                        "h_doubly_ap_DCS_reconstructed_cmFrame", true    // ROOT object name, auto-generated if omitted
                                                                        );

    // define some prompt and random windows (in nanoseconds)
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    // create/initialize the tree
    t.CreateBranches(HistFac.makeTTree("tree"));

}

void scratch_damaurer_MC_ppi0_pgg::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event); 

    //Fill e.g. the polar angle distribution into a histogram
    const auto& candidates = event.Reconstructed().Candidates;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;
    vector<double> NeuCanCluSize;
    vector<double> NeuCanCaloE;
    vector<double> ChaCanCluSize;
    vector<double> ChaCanCaloE;
    vector<double> The;
    vector<double> Phi;
    vector<double> PhotonThe;
    vector<double> PhotonTheCB;
    vector<double> PhotonTheTAPS;
    vector<double> PhotonPhi;
    vector<double> PhotonE;

    for(const auto& cand : candidates.get_iter()) {
        h_VetoEnergies->Fill(cand->VetoEnergy);
        The.push_back(cand->Theta*radtodeg);
        Phi.push_back(cand->Phi*radtodeg);
        if(cand->VetoEnergy == 0) {
            neutral.emplace_back(cand);
            NeuCanCluSize.push_back(cand->ClusterSize);
            NeuCanCaloE.push_back(cand->CaloEnergy);
        }       
        else{
            charged.emplace_back(cand);
            ChaCanCluSize.push_back(cand->ClusterSize);
            ChaCanCaloE.push_back(cand->CaloEnergy);
        }

    }

    //getting access to the pi0 decay photons

    double pi0gg_IM=-100;
    TParticle pi0gg2g;
    TParticle g1;
    TParticle g2;
    /*
    vector<bool> g1CB;
    vector<bool> g1BaF2;
    vector<bool> g1PbWO;
    vector<bool> g1undetected;
    vector<bool> g2CB;
    vector<bool> g2BaF2;
    vector<bool> g2PbWO;
    vector<bool> g2undetected;
    */
    bool isg1CB = false;
    bool isg1BaF2 = false;
    bool isg1PbWO = false;
    //bool isg1undetected = false;
    bool isg2CB = false;
    bool isg2BaF2 = false;
    bool isg2PbWO = false;
    //bool isg2undetected = false;

    if(NeuCanCaloE.size() == 2){

        g1 = TParticle(ParticleTypeDatabase::Photon, neutral[0]);
        g2 = TParticle(ParticleTypeDatabase::Photon, neutral[1]);

        auto g1channel = neutral[0]->FindCaloCluster()->CentralElement;
        auto g2channel = neutral[1]->FindCaloCluster()->CentralElement;

        if(neutral[0]->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
        {isg1CB = true;}
        if(neutral[0]->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            if(taps->IsPbWO4(g1channel)){isg1PbWO = true;}else{isg1BaF2 = true;}
        }
        //if(!isg1CB && !isg1BaF2 && !isg1PbWO)
        //{isg1undetected = true;}

        if(neutral[1]->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
        {isg2CB = true;}
        if(neutral[1]->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
            if(taps->IsPbWO4(g2channel)){isg2PbWO = true;}else{isg2BaF2 = true;}
        }
        //if(!isg2CB && !isg2BaF2 && !isg2PbWO)
        //{isg2undetected = true;}

        double pi0gg2g_time = 0;
        for (const auto& photon : neutral) {
            pi0gg2g+=TParticle(ParticleTypeDatabase::Photon, photon);
            pi0gg2g_time+=photon->Time;
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                PhotonTheCB.push_back(photon->Theta*radtodeg);
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                PhotonTheTAPS.push_back(photon->Theta*radtodeg);

            PhotonThe.push_back(photon->Theta*radtodeg);
            PhotonPhi.push_back(photon->Phi*radtodeg);
            PhotonE.push_back(photon->CaloEnergy);
        }
        pi0gg_IM = pi0gg2g.M();
        /*
        g1CB.push_back(isg1CB);
        g1BaF2.push_back(isg1BaF2);
        g1PbWO.push_back(isg1PbWO);
        g1undetected.push_back(isg1undetected);
        g2CB.push_back(isg2CB);
        g2BaF2.push_back(isg2BaF2);
        g2PbWO.push_back(isg2PbWO);
        g2undetected.push_back(isg2undetected);
        */
    }

    //Looping over the taggerhits

    TParticle proton;
    for (const auto& taggerhit : event.Reconstructed().TaggerHits) { // Event loop

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));

        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());

        proton = TParticle(ParticleTypeDatabase::Proton,(TLorentzVector)(InitialPhotonVec+InitialProtonVec-(g1+g2)));

        TLorentzVector Lpion = (TLorentzVector)(g1+g2);
        TLorentzVector Linitial = (TLorentzVector)(InitialPhotonVec+InitialProtonVec);
        TLorentzVector LpionBoosted = Lpion;
        LpionBoosted.Boost(-Linitial.BoostVector());

        //Adding selections & filling histograms

        const double weight = promptrandom.FillWeight();

        h_TaggerTime->Fill(taggerhit.Time, weight);
        h_InitialBeamE->Fill(InitialPhotonVec.E()/1000,weight);

        if(NeuCanCaloE.size() == 2){

            h_ProtonIm->Fill(proton.M(),weight);
            t.ProtonIm = proton.M();

            if(proton.M()<(mp+sigma) && proton.M()>(mp-sigma)){

                h_PhotonIm->Fill(pi0gg_IM,weight);
                t.PhotonIm = pi0gg_IM;

                for (unsigned int i=0; i<The.size(); i++){
                    h_PolarAngles->Fill(The[i],weight);
                }

                for (unsigned int i=0; i<Phi.size(); i++){
                    h_AzimuthAngles->Fill(Phi[i],weight);
                }

                for (unsigned int i=0; i<PhotonThe.size(); i++){
                    h_PhotonPolarAngles->Fill(PhotonThe[i],weight);
                    h_PhotonETheta->Fill(PhotonThe[i],PhotonE[i],weight);
                    //if (PhotonThe[i] < 20){PhotonPolarAngles_TAPS += 1.0;}
                    //if (PhotonThe[i] >= 20){PhotonPolarAngles_CB += 1.0;}
                    h_doubly_ap_DCS_reconstructed_lab->Fill(cos(Lpion.Theta()),Linitial.M(),weight);
                    h_doubly_ap_DCS_reconstructed_cmFrame->Fill(cos(Linitial.Angle(LpionBoosted.Vect())),Linitial.M(),weight);
                }

                for (unsigned int i=0; i<PhotonTheCB.size(); i++){
                    h_PhotonPolarAngles_CB->Fill(PhotonTheCB[i],weight);
                    h_PhotonETheta_CB->Fill(PhotonTheCB[i],PhotonE[i],weight);
                }

                for (unsigned int i=0; i<PhotonTheTAPS.size(); i++){
                    h_PhotonPolarAngles_TAPS->Fill(PhotonTheTAPS[i],weight);
                    h_PhotonETheta_TAPS->Fill(PhotonTheTAPS[i],PhotonE[i],weight);
                }

                h_ProtonPolarAngles->Fill(proton.Theta()*radtodeg,weight);
                h_ProtonETheta->Fill(proton.Theta()*radtodeg,proton.E,weight);
                if (proton.Theta()*radtodeg < 20){ProtonPolarAngles_TAPS += 1.0;}
                if (proton.Theta()*radtodeg >= 20 && proton.Theta()*radtodeg < 160){ProtonPolarAngles_CB += 1.0;}

                for (unsigned int i=0; i<PhotonPhi.size(); i++){
                    h_PhotonAzimuthAngles->Fill(PhotonPhi[i],weight);
                    h_PhotonEPhi->Fill(PhotonPhi[i],PhotonE[i],weight);
                }

                h_ProtonAzimuthAngles->Fill(proton.Phi()*radtodeg,weight);
                h_ProtonEPhi->Fill(proton.Phi()*radtodeg,proton.E,weight);

                h_ProtonPeak->Fill(proton.M(),weight);

                //Filling the statistical histogram

                if (isg1CB && isg2CB){count_first+=promptrandom.FillWeight();}
                if ((isg1CB && isg2BaF2)||(isg2CB && isg1BaF2)){count_second+=promptrandom.FillWeight();}
                if ((isg1CB && isg2PbWO)||(isg2CB && isg1PbWO)){count_third+=promptrandom.FillWeight();}
                if (isg1BaF2 && isg2BaF2){count_fourth+=promptrandom.FillWeight();}
                if ((isg1BaF2 && isg2PbWO)||(isg2BaF2 && isg1PbWO)){count_fifth+=promptrandom.FillWeight();}
                if (isg1PbWO && isg2PbWO){count_sixth+=promptrandom.FillWeight();}
                //if ((isg1PbWO && isg2undetected)||(isg2PbWO && isg1undetected)){count_seventh++;}
                //if ((isg1BaF2 && isg2undetected)||(isg2BaF2 && isg1undetected)){count_eigth++;}
                //if ((isg1CB && isg2undetected)||(isg2CB && isg1undetected)){count_nineth++;}
                //if (isg1undetected && isg2undetected){count_tenth++;}

                weight_res += weight;
            }
        }

        //Filling the tree for further analysis

        t.TaggW = promptrandom.FillWeight();
        t.nClusters = event.Reconstructed().Clusters.size();
        t.PolarAngles = The;
        t.PhotonPolarAngles = PhotonThe;
        t.PhotonPolarAnglesCB = PhotonTheCB;
        t.PhotonPolarAnglesTAPS = PhotonTheTAPS;
        t.ProtonPolarAngles = proton.Theta()*radtodeg;
        t.AzimuthAngles = Phi;
        t.PhotonAzimuthAngles = PhotonPhi;
        t.ProtonAzimuthAngles = proton.Phi()*radtodeg;
        t.ProtonE = proton.E;
        t.PhotonE = PhotonE;      
        t.ProtonIm = proton.M();
        t.PhotonIm = pi0gg_IM;
        t.InitialBeamE = InitialPhotonVec.E();

        t.Tree->Fill();

    }      

    h_nClusters->Fill(event.Reconstructed().Clusters.size());

}

void scratch_damaurer_MC_ppi0_pgg::ShowResult()
{
    h_Reconstructed_Data_Statistics->SetBinContent(1*steps,(float)count_first/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(1*steps+1,(float)count_first/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(2*steps,(float)count_second/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(2*steps+1,(float)count_second/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(3*steps,(float)count_third/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(3*steps+1,(float)count_third/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(4*steps,(float)count_fourth/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(4*steps+1,(float)count_fourth/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(5*steps,(float)count_fifth/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(5*steps+1,(float)count_fifth/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(6*steps,(float)count_sixth/(float)weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(6*steps+1,(float)count_sixth/(float)weight_res);
    /*
    h_Reconstructed_Data_Statistics->SetBinContent(7*steps,weight_res);
    h_Reconstructed_Data_Statistics->SetBinContent(7*steps,count_seventh);
    h_Reconstructed_Data_Statistics->SetBinContent(8*steps,count_eigth);
    h_Reconstructed_Data_Statistics->SetBinContent(9*steps,count_nineth);
    h_Reconstructed_Data_Statistics->SetBinContent(10*steps,count_tenth);
    */

    ant::canvas(GetName()+": tagger stuff")
            << h_TaggerTime
            << h_nClusters
            << h_VetoEnergies
            << h_InitialBeamE
            << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Polarangles whole setup")
            << h_PolarAngles
            << h_PhotonPolarAngles
            << h_ProtonPolarAngles
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Photon Polarangles CB/TAPS")
            << h_PhotonPolarAngles_CB
            << h_PhotonPolarAngles_TAPS
            << h_Reconstructed_Data_Statistics
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Azimuthangles")
            << h_AzimuthAngles
            << h_PhotonAzimuthAngles
            << h_ProtonAzimuthAngles
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Energy vs angles")
            << drawoption("pcolz")
            << h_PhotonETheta
            << h_ProtonETheta
            << h_PhotonEPhi
            << h_ProtonEPhi
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Photon energy vs angles CB/TAPS")
            << drawoption("pcolz")
            << h_PhotonETheta_CB
            << h_PhotonETheta_TAPS
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Particle Invariant masses")
            << h_PhotonIm
            << h_ProtonIm
            << h_ProtonPeak
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Statistics histogram")
            << h_Reconstructed_Data_Statistics
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Cross section checks")
            << drawoption("Surf")
            << h_doubly_ap_DCS_reconstructed_lab
            << h_doubly_ap_DCS_reconstructed_cmFrame
            << endc; // actually draws the canvas

}

void scratch_damaurer_MC_ppi0_pgg::Finish()
{

    int max_particles_detected = h_PolarAngles->GetEntries();
    int Photons_max_detected = h_PhotonPolarAngles->GetEntries();
    int Protons_max_detected = h_ProtonPolarAngles->GetEntries();

    PhotonPolarAngles_TAPS_frac = (h_PhotonPolarAngles_TAPS->GetEntries())/Photons_max_detected;
    PhotonPolarAngles_CB_frac = (h_PhotonPolarAngles_CB->GetEntries())/Photons_max_detected;
    ProtonPolarAngles_TAPS_frac = ProtonPolarAngles_TAPS/Protons_max_detected;
    ProtonPolarAngles_CB_frac = ProtonPolarAngles_CB/Protons_max_detected;

    PhotonPolarAngles_TAPS_err = sqrt((PhotonPolarAngles_TAPS_frac*(1-PhotonPolarAngles_TAPS_frac))/Photons_max_detected);
    PhotonPolarAngles_CB_err = sqrt((PhotonPolarAngles_CB_frac*(1-PhotonPolarAngles_CB_frac))/Photons_max_detected);
    ProtonPolarAngles_TAPS_err = sqrt((ProtonPolarAngles_TAPS_frac*(1-ProtonPolarAngles_TAPS_frac))/Protons_max_detected);
    ProtonPolarAngles_CB_err = sqrt((ProtonPolarAngles_CB_frac*(1-ProtonPolarAngles_CB_frac))/Protons_max_detected);

    cout << "" << endl;
    cout << "Finished processing events, total #events: " << h_nClusters->GetEntries() << endl;
    cout << "Integrated amount of found clusters in total: " << h_nClusters->Integral() << endl;
    cout << "Detection efficiency: " << max_particles_detected/max_particles << endl;
    cout << "Max number of detected photons: " << Photons_max_detected << endl;
    cout << "Max number of detected protons: " << Protons_max_detected << endl;
    cout << "Photon amount falling into TAPS: " << 100*PhotonPolarAngles_TAPS_frac << "% with error: " << 100*PhotonPolarAngles_TAPS_err << "%" << endl;
    cout << "Photon amount falling into CB: " << 100*PhotonPolarAngles_CB_frac << "% with error: " << 100*PhotonPolarAngles_CB_err << "%" << endl;
    cout << "Proton amount falling into TAPS: " << 100*ProtonPolarAngles_TAPS_frac << "% with error: " << 100*ProtonPolarAngles_TAPS_err << "%" << endl;
    cout << "Proton amount falling into CB: " << 100*ProtonPolarAngles_CB_frac << "% with error: " << 100*ProtonPolarAngles_CB_err << "%" << endl;
    cout << "" << endl;

    cout << "1) Both CB: " << (float)100*((float)count_first/(float)weight_res) << "%" << endl;
    cout << "2) CB & BaF2: " << (float)100*((float)count_second/(float)weight_res) << "%" << endl;
    cout << "3) CB & PbWO: " << (float)100*((float)count_third/(float)weight_res) << "%" << endl;
    cout << "4) Both BaF2: " << (float)100*((float)count_fourth/(float)weight_res) << "%" << endl;
    cout << "5) BaF2 & PbWO: " << (float)100*((float)count_fifth/(float)weight_res) << "%" << endl;
    cout << "6) Both PbWO: " << (float)100*((float)count_sixth/(float)weight_res) << "%" << endl;
    //cout << "7) PbWO & not detected: " << (float)100*((float)count_seventh/(float)weight_res) << "%" << endl;
    //cout << "8) BaF2 & not detected: " << (float)100*((float)count_eigth/(float)weight_res) << "%" << endl;
    //cout << "9) CB & not detected: " << (float)100*((float)count_nineth/(float)weight_res) << "%" << endl;
    //cout << "10) Both not detected: " << (float)100*((float)count_tenth/(float)weight_res) << "%" << endl;

}

AUTO_REGISTER_PHYSICS(scratch_damaurer_MC_ppi0_pgg)
