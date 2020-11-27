#include "pw_ppi0g_p3g_kinFit.h"
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
#include "utils/uncertainties/Interpolated.h"

#include <iostream>
#include <memory>

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

APLCON::Fit_Settings_t scratch_damaurer_pw_ppi0g_p3g_kinFit::MakeFitSettings(unsigned max_iterations)
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    return settings;
}

scratch_damaurer_pw_ppi0g_p3g_kinFit::scratch_damaurer_pw_ppi0g_p3g_kinFit(const string& name, OptionsPtr opts) :
    Physics(name, opts),

    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                           utils::UncertaintyModels::Interpolated::Type_t::MC,
                           make_shared<utils::UncertaintyModels::FitterSergey>())),
    fitter(nullptr, opts->Get<bool>("FitZVertex", true))

{

    fitter.SetZVertexSigma(3.0);
    taps = make_shared<expconfig::detector::TAPS_2013_11>(false, false, false);

    const BinSettings time_bins(200, 0,10);
    const BinSettings timeDiffCorTaggCB_bins(200, -5,10);
    const BinSettings timeDiffCorTaggTAPS_bins(200, -5,20);

    const BinSettings Veto_Energy_bins(1000, 0, 10);
    const BinSettings Calo_Energy_bins(1000, 0, 1200);
    const BinSettings initialBeamE_bins(200,0, 1600);
    const BinSettings cos_bins(100,-1,1);
    const BinSettings statistic_bins((nrCuts_total)*10,0,nrCuts_total);
    const BinSettings sqrt_S_bins(100,s_square_min,s_square_max);

    const BinSettings phi_bins(720, -180, 180);
    const BinSettings theta_bins(360, 0, 180);
    const BinSettings theta_bins_CB(140, 19, 160);
    const BinSettings theta_bins_TAPS(25, 0, 25);

    const BinSettings Im_proton_bins(100, 0, 1800);
    const BinSettings Im_omega_bins(100, 0, 1800);
    const BinSettings Im_pi0_bins(100, 0, 1000);

    const BinSettings Ekin_neu_bins(100, 0, 1600);
    const BinSettings Ekin_cha_bins(100, 0, 1000);

    const BinSettings E_photon_bins(100, 0, 1600);
    const BinSettings E_proton_bins(100, 900, 2000);
    const BinSettings E_omega_bins(100, 400, 1800);
    const BinSettings E_pi0_bins(100, 0, 1600);
    //const BinSettings cos_bins_BackToBack(1000,-1,1);
    //const BinSettings wpi0_Q_bins(100, 0, 1500);
    //const BinSettings proton_theta_bins(180, 0, 90);

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    auto hf_missingP_IM = new HistogramFactory("hf_missingP_IM", HistFac, "");
    auto hf_3g_IM = new HistogramFactory("hf_3g_IM", HistFac, "");
    auto hf_2gComb_IM = new HistogramFactory("hf_2gComb_IM", HistFac, "");
    auto hf_2gPi0_IM = new HistogramFactory("hf_2gPi0_IM", HistFac, "");
    auto hf_EvsTheta = new HistogramFactory("hf_EvsTheta", HistFac, "");
    auto hf_EvsPhi = new HistogramFactory("hf_EvsPhi", HistFac, "");
    auto hf_CrossSection = new HistogramFactory("hf_CrossSection", HistFac, "");
    auto hf_CaloEvsVetoE = new HistogramFactory("hf_CaloEvsVetoE", HistFac, "");
    auto hf_RecData_CandStat = new HistogramFactory("hf_RecData_CandStat", HistFac, "");
    auto hf_TimeDiffCorTaggCB = new HistogramFactory("hf_TimeDiffCorTaggCB", HistFac, "");
    auto hf_TimeDiffCorTaggTAPS = new HistogramFactory("hf_TimeDiffCorTaggTAPS", HistFac, "");
    auto hf_Extras = new HistogramFactory("hf_Extras", HistFac, "");

    auto hfTagger = new HistogramFactory("hfTagger", HistFac, "");
    auto hfOnlyEnergy = new HistogramFactory("All_the_candidate_energies", HistFac, "");
    auto hfOnlyAngles = new HistogramFactory("All_the_candidate_angles", HistFac, "");
    auto hf_wpi0g_EvsTheta = new HistogramFactory("hf_wpi0g_EvsTheta", HistFac, "");
    //auto hf_extras = new HistogramFactory("Some_additional_extra_hists", HistFac, "");

    for (unsigned int i=0; i<nrCuts_total; i++){
        /*
        stringstream ss;
        ss << i;
        string myString = ss.str();
        */

        h_AllCaloEvsVetoE_CB[i] = hf_CaloEvsVetoE->makeTH2D("Deposited Calo- VS Veto-energies in CB "+cuts[i],     // title
                                         "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                         Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                         "h_AllCaloEvsVetoE_CB_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_AllCaloEvsVetoE_TAPS[i] = hf_CaloEvsVetoE->makeTH2D("Deposited Calo- VS Veto-energies in TAPS "+cuts[i],     // title
                                         "Calo E [MeV]","Veto E [MeV]",  // xlabel, ylabel
                                         Calo_Energy_bins, Veto_Energy_bins,    // our binnings
                                         "h_AllCaloEvsVetoE_TAPS_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                         );

        h_neuEkinVSTheta[i] = hf_EvsTheta->makeTH2D("Photons Ekin vs Theta "+cuts[i], //title
                                                "#Theta [deg]","E_kin [MeV]",  // xlabel, ylabel
                                                theta_bins, Ekin_neu_bins,    // our binnings
                                                "h_neuEkinVSTheta_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                                );

        h_neuEkinVSPhi[i] = hf_EvsPhi->makeTH2D("Photons Ekin vs Phi "+cuts[i], //title
                                            "#Phi [deg]","E_kin [MeV]",  // xlabel, ylabel
                                            phi_bins, Ekin_neu_bins,    // our binnings
                                            "h_neuEkinVSPhi_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                            );

        h_chaEkinVSTheta[i] = hf_EvsTheta->makeTH2D("Charged Ekin vs Theta "+cuts[i], //title
                                                  "#Theta [deg]","E_kin [MeV]",  // xlabel, ylabel
                                                  theta_bins, Ekin_cha_bins,    // our binnings
                                                  "h_chaEkinVSTheta_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                                  );

        h_chaEkinVSPhi[i] = hf_EvsPhi->makeTH2D("Charged Ekin vs Phi "+cuts[i], //title
                                              "#Phi [deg]","E_kin [MeV]",  // xlabel, ylabel
                                              phi_bins, Ekin_cha_bins,    // our binnings
                                              "h_chaEkinVSPhi_"+cuts[i], true    // ROOT object name, auto-generated if omitted
                                              );

        h_beamE[i] = hf_Extras->makeTH1D("Initial photon beam energy "+cuts[i],     // title
                                                             "E_{photonbeam} [MeV]", "#",     // xlabel, ylabel
                                                             initialBeamE_bins,  // our binnings
                                                             "h_beamE_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_neuTimeDiffCorTaggCB[i] = hf_TimeDiffCorTaggCB->makeTH1D("Neutrals time diff to cor tagg-time in CB "+cuts[i],     // title
                                                             "t_{CB}-t_{Tagger} [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggCB_bins,  // our binnings
                                                             "h_neuTimeDiffCorTaggCB_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_neuTimeDiffCorTaggTAPS[i] = hf_TimeDiffCorTaggTAPS->makeTH1D("Neutrals time diff to cor tagg-time in TAPS "+cuts[i],     // title
                                                             "t_{TAPS}-t_{Tagger} [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggTAPS_bins,  // our binnings
                                                             "h_neuTimeDiffCorTaggTAPS_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_chaTimeDiffCorTaggCB[i] = hf_TimeDiffCorTaggCB->makeTH1D("Charged time diff to cor tagg-time in CB "+cuts[i],     // title
                                                             "(t_{CB}-t_{Tagger}) [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggCB_bins,  // our binnings
                                                             "h_chaTimeDiffCorTaggCB_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_chaTimeDiffCorTaggTAPS[i] = hf_TimeDiffCorTaggTAPS->makeTH1D("Charged time diff to cor tagg-time in TAPS "+cuts[i],     // title
                                                             "t_{TAPS}-t_{Tagger} [ns]", "#",     // xlabel, ylabel
                                                             timeDiffCorTaggTAPS_bins,  // our binnings
                                                             "h_chaTimeDiffCorTaggTAPS_"+cuts[i], true     // ROOT object name, auto-generated if omitted
                                                             );

    }

    for (unsigned int i=0; i<(nrCuts_total-1); i++){
        /*
        stringstream ss;
        ss << i;
        string myString = ss.str();
        */

        h_missingP_IM[i] = hf_missingP_IM->makeTH1D("Im(missingP) "+cuts[i+1],     // title
                                                             "IM(missingP) [MeV]", "#",     // xlabel, ylabel
                                                             Im_proton_bins,  // our binnings
                                                             "h_missingP_IM_"+cuts[i+1], true     // ROOT object name, auto-generated if omitted
                                                             );

        h_3g_IM[i] = hf_3g_IM->makeTH1D("IM(3g) "+cuts[i+1],     // title
                                            "IM(3g) [MeV]", "#",     // xlabel, ylabel
                                            Im_omega_bins,  // our binnings
                                            "h_3g_IM_"+cuts[i+1], true     // ROOT object name, auto-generated if omitted
                                            );

        h_2gComb_IM[i] = hf_2gComb_IM->makeTH1D("IM(2g) combinations "+cuts[i+1],     // title
                                             "IM(2g) [MeV]", "#",     // xlabel, ylabel
                                             Im_pi0_bins,  // our binnings
                                             "h_2gComb_IM_"+cuts[i+1], true     // ROOT object name, auto-generated if omitted
                                             );

        h_doublyDCScm_gp_wp[i] = hf_CrossSection->makeTH2D("gp_wp cross section in cm-frame "+cuts[i+1], //title
                                                                            "cos_Theta*","sqrt_S [MeV]", // xlabel, ylabel
                                                                            cos_bins, sqrt_S_bins,    // our binnings
                                                                            "h_doublyDCScm_gp_wp_"+cuts[i+1], true    // ROOT object name, auto-generated if omitted
                                                                            );

    }

    for (unsigned int i=0; i<nrCuts_pi0; i++){
        stringstream ss;
        ss << i;
        string myString = ss.str();

        h_2gPi0_IM[i] = hf_2gPi0_IM->makeTH1D("IM(2gPi0)",     // title
                                         "IM(2gPi0) [MeV]", "#",     // xlabel, ylabel
                                         Im_pi0_bins,  // our binnings
                                         "h_2gPi0_IM_"+myString, true     // ROOT object name, auto-generated if omitted
                                         );

    }

    h_RecData_Stat = hf_RecData_CandStat->makeTH1D("Amount of events after cuts:",     // title
                                     "Cuts", "#",     // xlabel, ylabel
                                     statistic_bins,  // our binnings
                                     "h_RecData_Stat", true    // ROOT object name, auto-generated if omitted
                                     );

    h_RecData_relStat = hf_RecData_CandStat->makeTH1D("rel. amount of events after cuts:",     // title
                                      "Cuts", "rel. #",     // xlabel, ylabel
                                      statistic_bins,  // our binnings
                                      "h_RecData_relStat", true    // ROOT object name, auto-generated if omitted
                                      );

    h_TaggerTime = hfTagger->makeTH1D("Tagger Time",     // title
                                    "t [ns]", "#",     // xlabel, ylabel
                                    time_bins,  // our binnings
                                    "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                    );
    
    h_InitialBeamE = hfOnlyEnergy->makeTH1D("Photon beam - Energy distribution",     // title
                                   "E [GeV]", "#",     // xlabel, ylabel
                                   initialBeamE_bins,  // our binnings
                                   "h_InitialBeamE"     // ROOT object name, auto-generated if omitted
                                   );

    h_nClusters = hfTagger->makeTH1D("Number of Clusters", "nClusters", "#", BinSettings(15), "h_nClusters");

    h_VetoEnergies = hfTagger->makeTH1D("Veto Energies", "E [MeV]", "#", Veto_Energy_bins, "h_VetoEnergies");


    h_NeuPolarAngles = hfOnlyAngles->makeTH1D("Neutral candidate polarangles",     // title
                                       "#Theta [deg]", "#",     // xlabel, ylabel
                                       theta_bins,  // our binnings
                                       "h_NeuPolarAngles", true     // ROOT object name, auto-generated if omitted
                                       );

    h_NeuAzimuthAngles = hfOnlyAngles->makeTH1D("Neutral candidate azimuthangles",     // title
                                       "#Phi [deg]", "#",     // xlabel, ylabel
                                       phi_bins,  // our binnings
                                       "h_NeuAzimuthAngles", true     // ROOT object name, auto-generated if omitted
                                       );

    h_ChaPolarAngles = hfOnlyAngles->makeTH1D("Charged particles polarangles",     // title
                                                    "#Theta [deg]", "#",     // xlabel, ylabel
                                                    theta_bins,  // our binnings
                                                    "h_ChaPolarAngles", true     // ROOT object name, auto-generated if omitted
                                                    );

    h_ChaAzimuthAngles = hfOnlyAngles->makeTH1D("Charged particles azimuthangles",     // title
                                                    "#Phi [deg]", "#",     // xlabel, ylabel
                                                    phi_bins,  // our binnings
                                                    "h_ChaAzimuthAngles", true     // ROOT object name, auto-generated if omitted
                                                    );

    h_3gPolarAngles = hfOnlyAngles->makeTH1D("Rec. Photons polarangles",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins,  // our binnings
                                             "h_3gPolarAngles", true     // ROOT object name, auto-generated if omitted
                                             );

    h_3gPolarAnglesCB = hfOnlyAngles->makeTH1D("Rec. Photons polarangles only in CB",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins_CB,  // our binnings
                                             "h_3gPolarAnglesCB", true     // ROOT object name, auto-generated if omitted
                                             );

    h_3gPolarAnglesTAPS = hfOnlyAngles->makeTH1D("Rec. Photons polarangles only in TAPS",     // title
                                             "#Theta [deg]", "#",     // xlabel, ylabel
                                             theta_bins_TAPS,  // our binnings
                                             "h_3gPolarAnglesTAPS", true     // ROOT object name, auto-generated if omitted
                                             );

    h_3gAzimuthAngles = hfOnlyAngles->makeTH1D("Rec. Photons azimuthangles",     // title
                                               "#Phi [deg]", "#",     // xlabel, ylabel
                                               phi_bins,  // our binnings
                                               "h_3gAzimuthAngles", true     // ROOT object name, auto-generated if omitted
                                               );

    h_3g_EvTheta_CB = hf_wpi0g_EvsTheta->makeTH2D("Photon energy vs theta only in CB", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins_CB, E_photon_bins,    // our binnings
                                       "h_3g_EvTheta_CB", true    // ROOT object name, auto-generated if omitted
                                       );

    h_3g_EvTheta_TAPS = hf_wpi0g_EvsTheta->makeTH2D("Photon energy vs theta only in TAPS", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins_TAPS, E_photon_bins,    // our binnings
                                       "h_3g_EvTheta_TAPS", true    // ROOT object name, auto-generated if omitted
                                       );

    h_p_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. proton energy vs theta", //title
                                                 "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                                 theta_bins, E_proton_bins,    // our binnings
                                                 "h_p_EvTheta", true    // ROOT object name, auto-generated if omitted
                                                 );

    h_w_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. omega meson energy vs theta", //title
                                                 "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                                 theta_bins, E_omega_bins,    // our binnings
                                                 "h_w_EvTheta", true    // ROOT object name, auto-generated if omitted
                                                 );

    h_wg_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. omega decay-photon energy vs theta", //title
                                       "Theta [deg]","E [MeV]",  // xlabel, ylabel
                                       theta_bins, E_photon_bins,    // our binnings
                                       "h_wg_EvTheta", true    // ROOT object name, auto-generated if omitted
                                       );

    h_wpi0_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. Pi0 energy vs theta", //title
                                        "Theta [deg]","E [MeV]", // xlabel, ylabel
                                        theta_bins, E_pi0_bins,    // our binnings
                                        "h_wpi0_EvTheta", true    // ROOT object name, auto-generated if omitted
                                        );

    h_wpi02g_EvTheta = hf_wpi0g_EvsTheta->makeTH2D("Rec. Pi0 decay-photons energy vs theta", //title
                                        "Theta [deg]","E [MeV]", // xlabel, ylabel
                                        theta_bins, E_pi0_bins,    // our binnings
                                        "h_wpi02g_EvTheta", true    // ROOT object name, auto-generated if omitted
                                        );

    // define some prompt and random windows (in nanoseconds)
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    // create/initialize the tree
    t.CreateBranches(HistFac.makeTTree("tree"));

}

void scratch_damaurer_pw_ppi0g_p3g_kinFit::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event); 

    // set fitter uncertainty models
    fitter.SetUncertaintyModel(fit_model);

    const auto& data = event.Reconstructed();
    const auto& candidates = data.Candidates;

    //Fill e.g. the polar angle distribution into a histogram
    TCandidatePtrList all;
    TCandidatePtrList neutral;
    TCandidatePtrList charged;

    vector<bool> chaCanInCB;
    vector<bool> chaCanInTAPS;
    vector<bool> neuCanInCB;
    vector<bool> neuCanInTAPS;

    vector<double> neuCanCluSize;
    vector<double> neuCanCaloE;
    vector<double> neuCanVetoE;
    vector<double> neuCanTime;
    vector<double> neuThe;
    vector<double> neuPhi;
    vector<double> chaCanCluSize;
    vector<double> chaCanCaloE;
    vector<double> chaCanVetoE;
    vector<double> chaCanTime;
    vector<double> chaThe;
    vector<double> chaPhi;

    for(const auto& cand : candidates.get_iter()) {
        //h_VetoEnergies->Fill(cand->VetoEnergy);
        all.emplace_back(cand);

        if(cand->VetoEnergy <= vetoEthreshold) {
            neutral.emplace_back(cand);
            neuThe.push_back(cand->Theta*radtodeg);
            neuPhi.push_back(cand->Phi*radtodeg);
            neuCanCluSize.push_back(cand->ClusterSize);
            neuCanCaloE.push_back(cand->CaloEnergy);
            neuCanVetoE.push_back(cand->VetoEnergy);
            neuCanTime.push_back(cand->Time);

            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                neuCanInCB.push_back(true);
            }
            else{
                neuCanInCB.push_back(false);
            }
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                neuCanInTAPS.push_back(true);
            }
            else{
                neuCanInTAPS.push_back(false);
            }
        }       
        else{

            charged.emplace_back(cand);
            chaThe.push_back(cand->Theta*radtodeg);
            chaPhi.push_back(cand->Phi*radtodeg);
            chaCanCluSize.push_back(cand->ClusterSize);
            chaCanCaloE.push_back(cand->CaloEnergy);
            chaCanVetoE.push_back(cand->VetoEnergy);
            chaCanTime.push_back(cand->Time);

            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                chaCanInCB.push_back(true);
            }
            else{
                chaCanInCB.push_back(false);
            }
            if(cand->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS){
                chaCanInTAPS.push_back(true);
            }
            else{
                chaCanInTAPS.push_back(false);
            }
        }

    }

    //getting access to the pi0 decay photons

    double wpi0gOnly3g_IM=-100;
    TParticle g[neu_nrSel];
    TParticle proton_tmp;

    vector<double> PhotonThe;
    vector<double> PhotonTheCB;
    vector<double> PhotonTheTAPS;
    vector<double> PhotonPhi;
    vector<double> PhotonE;
    vector<double> wpi0g_times;

    if(neuCanCaloE.size() == neu_nrSel && chaCanCaloE.size() == cha_nrSel){

        for(int i = 0; i<neu_nrSel; i++){
            g[i] = TParticle(ParticleTypeDatabase::Photon, neutral[i]);
        }

        proton_tmp = TParticle(ParticleTypeDatabase::Proton, charged[0]);

        TParticle wpi0gOnly3g;
        double wpi0gOnly3g_time = 0;
        for (const auto& photon : neutral) {
            wpi0gOnly3g += TParticle(ParticleTypeDatabase::Photon, photon);
            wpi0gOnly3g_time+=photon->Time;
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB)
                PhotonTheCB.push_back(photon->Theta*radtodeg);
            if(photon->FindCaloCluster()->DetectorType == Detector_t::Type_t::TAPS)
                PhotonTheTAPS.push_back(photon->Theta*radtodeg);
            PhotonThe.push_back(photon->Theta*radtodeg);
            PhotonPhi.push_back(photon->Phi*radtodeg);
            PhotonE.push_back(photon->CaloEnergy);
        }
        wpi0gOnly3g_IM = wpi0gOnly3g.M();
        wpi0gOnly3g_time = wpi0gOnly3g_time/3.;
        wpi0g_times.push_back(wpi0gOnly3g_time);

    }

    //Looping over the taggerhits

    TParticle proton;
    TParticle omega_tmp;
    double corTaggTime;

    for (const auto& taggerhit : data.TaggerHits) { // Event loop

        corTaggTime = triggersimu.GetCorrectedTaggerTime(taggerhit);
        promptrandom.SetTaggerTime(corTaggTime);

        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
        TLorentzVector Linitial = (TLorentzVector)(InitialPhotonVec+InitialProtonVec);

        //Adding selections & filling histograms

        const double weight = promptrandom.FillWeight();

        stat[0]+=weight;

        h_TaggerTime->Fill(taggerhit.Time, weight);
        h_beamE[0]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[0]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuEkinVSTheta[0]->Fill(neuThe[i],neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllCaloEvsVetoE_CB[0]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggCB[0]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
            if(neuCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[0]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggTAPS[0]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[0]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaEkinVSTheta[0]->Fill(chaThe[i],chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllCaloEvsVetoE_CB[0]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggCB[0]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
            if(chaCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[0]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggTAPS[0]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
        }

        if(!(neuCanCaloE.size() == neu_nrSel && chaCanCaloE.size() == cha_nrSel))
            continue;

        stat[1]+=weight;

        TLorentzVector L3g;

        for(int i = 0; i<neu_nrSel; i++){
            L3g += (TLorentzVector)(g[i]);
        }

        omega_tmp = TParticle(ParticleTypeDatabase::Omega,(TLorentzVector)(L3g));

        TLorentzVector Lw_tmp = (TLorentzVector)(omega_tmp);
        TLorentzVector Lp_tmp = (TLorentzVector)(proton_tmp);

        TLorentzVector LmissingProton = Linitial-Lw_tmp;

        TLorentzVector Lp = Lp_tmp;
        TLorentzVector Lw = Lw_tmp;
        TLorentzVector Lw_boosted = Lw;
        Lw_boosted.Boost(-Linitial.BoostVector());

        h_missingP_IM[0]->Fill(LmissingProton.M(),weight);
        h_3g_IM[0]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[0]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[0]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[0]->Fill((g[1]+g[2]).M(),weight);

        h_doublyDCScm_gp_wp[0]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);
        h_beamE[1]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[1]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuEkinVSTheta[1]->Fill(neuThe[i],neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllCaloEvsVetoE_CB[1]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggCB[1]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
            if(neuCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[1]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggTAPS[1]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[1]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaEkinVSTheta[1]->Fill(chaThe[i],chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllCaloEvsVetoE_CB[1]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggCB[1]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
            if(chaCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[1]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggTAPS[1]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
        }

        if(!(LmissingProton.M() > (mpFit-1.5*sigmaMissingP) && LmissingProton.M() < (mpFit+1.5*sigmaMissingP)))
            continue;

        stat[2]+=weight;

        h_missingP_IM[1]->Fill(LmissingProton.M(),weight);
        h_3g_IM[1]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[1]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[1]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[1]->Fill((g[1]+g[2]).M(),weight);

        h_doublyDCScm_gp_wp[1]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);
        h_beamE[2]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[2]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuEkinVSTheta[2]->Fill(neuThe[i],neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllCaloEvsVetoE_CB[2]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggCB[2]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
            if(neuCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[2]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggTAPS[2]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[2]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaEkinVSTheta[2]->Fill(chaThe[i],chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllCaloEvsVetoE_CB[2]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggCB[2]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
            if(chaCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[2]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggTAPS[2]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
        }

        if(!(InitialPhotonVec.E()>=Omega_Ethreshold))
            continue;

        stat[3]+=weight;

        h_missingP_IM[2]->Fill(LmissingProton.M(),weight);
        h_3g_IM[2]->Fill(omega_tmp.M(),weight);

        h_2gComb_IM[2]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[2]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[2]->Fill((g[1]+g[2]).M(),weight);

        h_doublyDCScm_gp_wp[2]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);
        h_beamE[3]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[3]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuEkinVSTheta[3]->Fill(neuThe[i],neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllCaloEvsVetoE_CB[3]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggCB[3]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
            if(neuCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[3]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggTAPS[3]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[3]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaEkinVSTheta[3]->Fill(chaThe[i],chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllCaloEvsVetoE_CB[3]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggCB[3]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
            if(chaCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[3]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggTAPS[3]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
        }

        //LorentzBoost-Stuff with opening angle selections:

        /*
        TParticle wpi0_tmp1 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g1+g2));
        TParticle wpi0_tmp2 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g1+g3));
        TParticle wpi0_tmp3 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g2+g3));

        TParticle wg_tmp1 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g3));
        TParticle wg_tmp2 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g2));
        TParticle wg_tmp3 = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g1));

        TLorentzVector wpi0_tmp1_boosted = (TLorentzVector)wpi0_tmp1;
        TLorentzVector wpi0_tmp2_boosted = (TLorentzVector)wpi0_tmp2;
        TLorentzVector wpi0_tmp3_boosted = (TLorentzVector)wpi0_tmp3;

        TLorentzVector wg_tmp1_boosted = (TLorentzVector)wg_tmp1;
        TLorentzVector wg_tmp2_boosted = (TLorentzVector)wg_tmp2;
        TLorentzVector wg_tmp3_boosted = (TLorentzVector)wg_tmp3;

        wpi0_tmp1_boosted.Boost(-Lw.BoostVector());
        wpi0_tmp2_boosted.Boost(-Lw.BoostVector());
        wpi0_tmp3_boosted.Boost(-Lw.BoostVector());

        wg_tmp1_boosted.Boost(-Lw.BoostVector());
        wg_tmp2_boosted.Boost(-Lw.BoostVector());
        wg_tmp3_boosted.Boost(-Lw.BoostVector());

        double cmAngleComb1 = cos(wpi0_tmp1_boosted.Angle(wg_tmp1_boosted.Vect()));
        double cmAngleComb2 = cos(wpi0_tmp2_boosted.Angle(wg_tmp2_boosted.Vect()));
        double cmAngleComb3 = cos(wpi0_tmp3_boosted.Angle(wg_tmp3_boosted.Vect()));

        double array[] = {1.2,1.1,1.4,1.5,0.9};
        double min_i;
        double min = 10;
        int arraySize = sizeof(array)/sizeof(array[0]);

        for (int i=0; i<arraySize; i++){
            //array[i] = i;
            if(array[i] <= min){min = array[i]; min_i = i;}
            cout << "Arrayelement = " << array[i] << endl;
        }

        */

        TLorentzVector w_decayComb[neu_nrSel] = {(g[0] + g[1]),(g[0] + g[2]),(g[1] + g[2])};
        double inv_mass_comb[neu_nrSel] = {sqrt(Lw.M2()-w_decayComb[1].M2()-w_decayComb[2].M2()),sqrt(Lw.M2()-w_decayComb[0].M2()-w_decayComb[2].M2()),sqrt(Lw.M2()-w_decayComb[0].M2()-w_decayComb[1].M2())};
        double pi0_absMassDiff[neu_nrSel];

        TLorentzVector wpi0g[nr_pi0g];
        TLorentzVector wpi0;
        TLorentzVector wg;

        TParticle wpi0_tmp[neu_nrSel];
        TParticle wg_tmp[neu_nrSel];

        for(int i=0; i<neu_nrSel;i++){
            wpi0_tmp[i] = TParticle(ParticleTypeDatabase::Pi0,(TLorentzVector)(w_decayComb[i]));
            wg_tmp[i] = TParticle(ParticleTypeDatabase::Photon,(TLorentzVector)(g[neu_nrSel-1-i]));
            pi0_absMassDiff[i] = abs(inv_mass_comb[i]-mpi0);
        }

        if(pi0_absMassDiff[0] <= pi0_absMassDiff[1] && pi0_absMassDiff[0] <= pi0_absMassDiff[2]){wpi0 = wpi0_tmp[0]; wg = wg_tmp[0]; wpi0g[0] = g[0]; wpi0g[1] = g[1];}
        else if(pi0_absMassDiff[1] <= pi0_absMassDiff[0] && pi0_absMassDiff[1] <= pi0_absMassDiff[2]){wpi0 = wpi0_tmp[1]; wg = wg_tmp[1]; wpi0g[0] = g[0]; wpi0g[1] = g[2];}
        else{wpi0 = wpi0_tmp[2]; wg = wg_tmp[2]; wpi0g[0] = g[1]; wpi0g[1] = g[2];}

        h_2gPi0_IM[0]->Fill(wpi0.M(),weight);

        if(!(wpi0.M()>(mpi0Fit-3*sigmaPi0IM) && wpi0.M()<(mpi0Fit+3*sigmaPi0IM)))
            continue;

        stat[4]+=weight;

        h_2gPi0_IM[1]->Fill(wpi0.M(),weight);

        h_missingP_IM[3]->Fill(LmissingProton.M(),weight);
        h_3g_IM[3]->Fill(omega_tmp.M(),weight);
        h_2gComb_IM[3]->Fill((g[0]+g[1]).M(),weight);
        h_2gComb_IM[3]->Fill((g[0]+g[2]).M(),weight);
        h_2gComb_IM[3]->Fill((g[1]+g[2]).M(),weight);

        h_doublyDCScm_gp_wp[3]->Fill(cos(Linitial.Angle(Lw_boosted.Vect())),Linitial.M(),weight);
        h_beamE[4]->Fill(InitialPhotonVec.E(),weight);

        for (unsigned int i=0; i<neutral.size(); i++){
            h_neuEkinVSPhi[4]->Fill(neuPhi[i],neuCanCaloE[i],weight);
            h_neuEkinVSTheta[4]->Fill(neuThe[i],neuCanCaloE[i],weight);
            if(neuCanInCB[i]){
                h_AllCaloEvsVetoE_CB[4]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggCB[4]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
            if(neuCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[4]->Fill(neuCanCaloE[i],neuCanVetoE[i],weight);
                h_neuTimeDiffCorTaggTAPS[4]->Fill(neuCanTime[i]-corTaggTime,weight);
            }
        }
        for (unsigned int i=0; i<charged.size(); i++){
            h_chaEkinVSPhi[4]->Fill(chaPhi[i],chaCanCaloE[i],weight);
            h_chaEkinVSTheta[4]->Fill(chaThe[i],chaCanCaloE[i],weight);
            if(chaCanInCB[i]){
                h_AllCaloEvsVetoE_CB[4]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggCB[4]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
            if(chaCanInTAPS[i]){
                h_AllCaloEvsVetoE_TAPS[4]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
                h_chaTimeDiffCorTaggTAPS[4]->Fill(chaCanTime[i]-corTaggTime,weight);
            }
        }

        /*
        if(!(charged[0]->VetoEnergy>0.9))
                continue;

        if(!(proton_tmp.Theta()*radtodeg <= 80))
            continue;

        for (unsigned int i=0; i<chaThe.size(); i++){
            if(chaThe[i]>=20 && chaThe[i]<=160)
                h_ChaCaloVSVetoEnergies_CB[3]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
            if(chaThe[i]<20 && chaThe[i]>=2)
                h_ChaCaloVSVetoEnergies_TAPS[3]->Fill(chaCanCaloE[i],chaCanVetoE[i],weight);
        }

        proton = TParticle(ParticleTypeDatabase::Proton,(TLorentzVector)(proton_tmp));

        h_wOnly3g_Im[4]->Fill(omega_tmp.M(),weight);
        h_missingProton_Im[4]->Fill(LmissingProton.M(),weight);
        h_2gMassComb[4]->Fill((g[0]+g[1]).M(),weight);
        h_2gMassComb[4]->Fill((g[0]+g[2]).M(),weight);
        h_2gMassComb[4]->Fill((g[1]+g[2]).M(),weight);

        */

        h_p_EvTheta->Fill(Lp.Theta()*radtodeg,Lp.E(),weight);
        h_w_EvTheta->Fill(Lw.Theta()*radtodeg,Lw.E(),weight);
        h_wg_EvTheta->Fill(wg.Theta()*radtodeg,wg.E(),weight);
        h_wpi0_EvTheta->Fill(wpi0.Theta()*radtodeg,wpi0.E(),weight);

        for(int i=0; i<nr_pi0g;i++){
            h_wpi02g_EvTheta->Fill(wpi0g[i].Theta()*radtodeg,wpi0g[i].E(),weight);
        }

        //weight_res += weight;

        //Filling the tree for further analysis

        t.TaggW = promptrandom.FillWeight();
        t.nClusters = data.Clusters.size();
        t.PhotonAzimuthAngles = PhotonPhi;
        t.PhotonPolarAngles = PhotonThe;
        t.PhotonPolarAnglesCB = PhotonTheCB;
        t.PhotonPolarAnglesTAPS = PhotonTheTAPS;
        t.Tree->Fill();     
    }
    h_nClusters->Fill(data.Clusters.size());

}

void scratch_damaurer_pw_ppi0g_p3g_kinFit::ShowResult()
{    

    int lower_edge = 0;
    int upper_edge = nrCuts_total;
    int number_of_bins = nrCuts_total*10;
    int steps = (int)(number_of_bins/(upper_edge-lower_edge));

    for (unsigned int i=0; i<nrCuts_total; i++){
        cout << "Amount events after " << i << " applied cuts: " << stat[i] << endl;
        h_RecData_Stat->SetBinContent(i*steps+1, stat[i]);
        h_RecData_Stat->GetXaxis()->SetBinLabel(i*steps+1,cuts[i].c_str());
    }

    for (unsigned int i=0; i<(nrCuts_total-1); i++){
        cout << "Relative amount of events after " << (i+1) << " applied cuts: " << stat[i+1]/stat[0] << endl;
        h_RecData_relStat->SetBinContent((i+1)*steps+1, stat[i+1]/stat[0]);
        h_RecData_relStat->GetXaxis()->SetBinLabel((i+1)*steps+1,cuts[i+1].c_str());
    }

    ant::canvas(GetName()+": Event-Statistics after applied cuts:")
            << h_RecData_Stat
            << h_RecData_relStat
            << endc; // actually draws the canvas

    ant::canvas c1(GetName()+": IM(missingP)");
    for (unsigned int i=0; i<(nrCuts_total-1); i++){
            c1 << h_missingP_IM[i];
    }
            c1 << endc; // actually draws the canvas

    ant::canvas c2(GetName()+": IM(3neu)");
    for (unsigned int i=0; i<(nrCuts_total-1); i++){
            c2 << h_3g_IM[i];
    }
            c2 << endc; // actually draws the canvas


    ant::canvas c3(GetName()+": IM(2neu) combinations");
    for (unsigned int i=0; i<(nrCuts_total-1); i++){
            c3 << h_2gComb_IM[i];
    }
            c3 << endc; // actually draws the canvas

    ant::canvas c4(GetName()+": gp_wp cross section in cm-frame");
            c4 << drawoption("Surf");
    for (unsigned int i=0; i<(nrCuts_total-1); i++){
            c4 << h_doublyDCScm_gp_wp[i];
    }
            c4 << endc; // actually draws the canvas

    ant::canvas c5(GetName()+": IM(2neuPi0)");
    for (unsigned int i=0; i<nrCuts_pi0; i++){
            c5 << h_2gPi0_IM[i];
    }
            c5 << endc; // actually draws the canvas

    ant::canvas c6(GetName()+": Deposited Calo- VS Veto-energies in CB");
            c6 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c6 << h_AllCaloEvsVetoE_CB[i];
    }
            c6 << endc; // actually draws the canvas

    ant::canvas c7(GetName()+": Deposited Calo- VS Veto-energies in TAPS");
            c7 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c7 << h_AllCaloEvsVetoE_TAPS[i];
    }
            c7 << endc; // actually draws the canvas

    ant::canvas c8(GetName()+": Initial photonbeam energy");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c8 << h_beamE[i];
    }
            c8 << endc; // actually draws the canvas

/*

    ant::canvas c9(GetName()+": Photons Ekin vs Theta");
            c9 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c9 << h_neuEkinVSTheta[i];
    }
            c9 << endc; // actually draws the canvas

    ant::canvas c10(GetName()+": charged Ekin vs Theta");
            c10 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c10 << h_chaEkinVSTheta[i];
    }
            c10 << endc; // actually draws the canvas

    ant::canvas c11(GetName()+": Photons Ekin vs Phi");
            c11 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c11 << h_neuEkinVSPhi[i];
    }
            c11 << endc; // actually draws the canvas

    ant::canvas c12(GetName()+": charged Ekin vs Phi");
            c12 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c12 << h_chaEkinVSPhi[i];
    }
            c12 << endc; // actually draws the canvas

    */

    ant::canvas c17(GetName()+": charged time-diff to cor. Tagg time in CB");
            c17 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c17 << h_chaTimeDiffCorTaggCB[i];
    }
            c17 << endc; // actually draws the canvas

    ant::canvas c18(GetName()+": charged time-diff to cor. Tagg time in TAPS");
            c18 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c18 << h_chaTimeDiffCorTaggTAPS[i];
    }
            c18 << endc; // actually draws the canvas

    ant::canvas c19(GetName()+": neutrals time-diff to cor. Tagg time in CB");
            c19 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c19 << h_neuTimeDiffCorTaggCB[i];
    }
            c19 << endc; // actually draws the canvas

    ant::canvas c20(GetName()+": neutrals time-diff to cor. Tagg time in TAPS");
            c20 << drawoption("pcolz");
    for (unsigned int i=0; i<nrCuts_total; i++){
            c20 << h_neuTimeDiffCorTaggTAPS[i];
    }
            c20 << endc; // actually draws the canvas

   /*
   ant::canvas(GetName()+": All candidates deposited Calo vs Veto energy")
           << drawoption("pcolz")
           << h_AllCaloVSVetoEnergies_CB
           << h_AllCaloVSVetoEnergies_TAPS
           << endc; // actually draws the canvas

   ant::canvas(GetName()+": Deposited Calo vs Veto energy of expected uncharged candidates.")
           << drawoption("pcolz")
           << h_NeuCaloVSVetoEnergies_CB
           << h_NeuCaloVSVetoEnergies_TAPS
           << endc; // actually draws the canvas

   for (unsigned int i=0; i<nrCuts_IM; i++){

       stringstream ss;
       ss << i;
       string myString = ss.str();
       ant::canvas(GetName()+": Deposited Calo vs Veto energy of expected charged candidates after "+myString+ " applied cuts")
               << drawoption("pcolz")
               << h_ChaCaloVSVetoEnergies_CB[i]
               << h_ChaCaloVSVetoEnergies_TAPS[i]
               << endc; // actually draws the canvas
   }

    for (unsigned int i=0; i<nrCuts_pi0; i++){

        stringstream ss;
        ss << i;
        string myString = ss.str();

        ant::canvas(GetName()+": Rec. pi0: Invariant mass after "+myString+ " applied cuts")
                << h_Pi0Only2g_Im[i]
                << endc; // actually draws the canvas

    }

    ant::canvas(GetName()+": Rec. Omega & protons: EvsTheta distribution")
            << drawoption("pcolz")
            << h_w_EvTheta
            << h_p_EvTheta
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Rec. wpi0g decay: EvsTheta distribution")
            << drawoption("pcolz")
            << h_wg_EvTheta
            << h_wpi0_EvTheta
            << h_wpi02g_EvTheta
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Cross section checks")
            << drawoption("Surf")
            << h_doubly_wp_DCS_reconstructed_lab
            << h_doubly_wp_DCS_reconstructed_cmFrame
            << endc; // actually draws the canvas
*/


}

void scratch_damaurer_pw_ppi0g_p3g_kinFit::Finish()
{
    cout << "please work!" << endl;
/*
    int max_particles_detected = h_ALL_PolarAngles->GetEntries();
    int Photons_max_detected = h_3gPolarAngles->GetEntries();
    int MissingCHA_max_detected = h_missingChaPolarAngles->GetEntries();

    PhotonPolarAngles_TAPS_frac = (h_3gPolarAnglesTAPS->GetEntries())/Photons_max_detected;
    PhotonPolarAngles_CB_frac = (h_3gPolarAnglesCB->GetEntries())/Photons_max_detected;
    //ProtonPolarAngles_TAPS_frac = ProtonPolarAngles_TAPS/MissingCHA_max_detected;
    //ProtonPolarAngles_CB_frac = ProtonPolarAngles_CB/MissingCHA_max_detected;

    PhotonPolarAngles_TAPS_err = sqrt((PhotonPolarAngles_TAPS_frac*(1-PhotonPolarAngles_TAPS_frac))/Photons_max_detected);
    PhotonPolarAngles_CB_err = sqrt((PhotonPolarAngles_CB_frac*(1-PhotonPolarAngles_CB_frac))/Photons_max_detected);
    //ProtonPolarAngles_TAPS_err = sqrt((ProtonPolarAngles_TAPS_frac*(1-ProtonPolarAngles_TAPS_frac))/MissingCHA_max_detected);
    //ProtonPolarAngles_CB_err = sqrt((ProtonPolarAngles_CB_frac*(1-ProtonPolarAngles_CB_frac))/MissingCHA_max_detected);

    cout << "" << endl;
    cout << "Finished processing events, total #events: " << h_nClusters->GetEntries() << endl;
    cout << "Integrated amount of found clusters in total: " << h_nClusters->Integral() << endl;
    cout << "Detection efficiency: " << max_particles_detected/max_particles << endl;
    cout << "Max number of detected photons: " << Photons_max_detected << endl;
*/
    /*
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
    */

}

AUTO_REGISTER_PHYSICS(scratch_damaurer_pw_ppi0g_p3g_kinFit)
