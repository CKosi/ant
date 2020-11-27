#include "Tutorial_kinFit.h"
#include "utils/uncertainties/Interpolated.h"
#include "utils/uncertainties/FitterSergey.h"

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

APLCON::Fit_Settings_t scratch_damaurer_Tutorial_kinFit::MakeFitSettings(unsigned max_iterations)
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    return settings;
}

scratch_damaurer_Tutorial_kinFit::scratch_damaurer_Tutorial_kinFit(const string& name, OptionsPtr opts) :
    Physics(name, opts),

    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                             utils::UncertaintyModels::Interpolated::Type_t::MC,
                             make_shared<utils::UncertaintyModels::FitterSergey>())),
    fitter(nullptr, opts->Get<bool>("FitZVertex", true))
{

    fitter.SetZVertexSigma(3.0);

    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    pid_detector = ExpConfig::Setup::GetDetector<expconfig::detector::PID>();
    veto_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPSVeto>();

    const BinSettings bins_nClusters(20);
    const BinSettings tagger_time_bins(2000, -200, 200);

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    h_nClusters = HistFac.makeTH1D("Number of Clusters", // title
                                   "nClusters","#",      // xlabel, ylabel
                                   bins_nClusters,       // our binnings, may write directly BinSettings(10) here
                                   "h_nClusters"         // ROOT object name, auto-generated if omitted
                                   );

    h_TaggerTime = HistFac.makeTH1D("Tagger Time",     // title
                                        "t [ns]", "#",     // xlabel, ylabel
                                        tagger_time_bins,  // our binnings
                                        "h_TaggerTime"     // ROOT object name, auto-generated if omitted
                                        );

    //hist = HistFac.makeTH1D(" Accepted Events", "step", "#", BinSettings(10), "steps");

    // define some prompt and random windows (in nanoseconds)
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});

    // create/initialize the tree
    t.CreateBranches(HistFac.makeTTree("t"));  

}

void scratch_damaurer_Tutorial_kinFit::ProcessEvent(const TEvent& event, manager_t&)
{

    triggersimu.ProcessEvent(event);

    // set fitter uncertainty models
    fitter.SetUncertaintyModel(fit_model);

    const auto& data = event.Reconstructed();
    const auto& candidates = data.Candidates;

    //Set up particle combinations for the kinFit
    utils::ProtonPhotonCombs proton_photons(candidates);
    particle_combs_t protphotcombs = proton_photons();

    h_nClusters->Fill(data.Clusters.size());

    double best_probability = std_ext::NaN;

    for(auto& taggerhit : data.TaggerHits) {

        //Get corrected Taggertime
        double cortagtime = triggersimu.GetCorrectedTaggerTime(taggerhit);
        promptrandom.SetTaggerTime(cortagtime);

        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;


        const double TaggWeight = promptrandom.FillWeight();

        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();
        TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());

        h_TaggerTime->Fill(taggerhit.Time, TaggWeight);

        t.TaggW = TaggWeight;
        t.nClusters = data.Clusters.size();
        t.Tree->Fill();

        //Particles for the later kinFit
        TParticlePtr proton;
        TParticleList photons;

        //Photon-proton comb:

        //Performing the kinFit

        APLCON::Result_t fitresult = fitter.DoFit(taggerhit.PhotonEnergy, proton, photons);

        // check if the fit converged
        if (fitresult.Status != APLCON::Result_Status_t::Success)
            continue;

        // check if we found a better probability for this fit and copy it if true, continue otherwise
        if (!std_ext::copy_if_greater(best_probability, fitresult.Probability))
            continue;

        // retrieve the fitted photon and proton information as well as the number of iterations
        auto fitted_proton = fitter.GetFittedProton();
        auto fitted_photons = fitter.GetFittedPhotons();
        auto iterations = fitresult.NIterations;


    } 

}

void scratch_damaurer_Tutorial_kinFit::ShowResult()
{

    // ShowResult is called after processing of events has finished,
    // and interactive mode (aka non-batchmode) is chosen

    // ant::canvas nice wrapper around TCanvas

    ant::canvas(GetName()+": Basic plots")
            << h_nClusters
            << h_TaggerTime
            << TTree_drawable(t.Tree, "nClusters >> (20,0,20)", "TaggW")
            << endc; // actually draws the canvas
}

void scratch_damaurer_Tutorial_kinFit::Finish()
{

    cout << "Finished processing events, total #events: " << h_nClusters->GetEntries() << endl;
    cout << "Integrated amount of found clusters in total: " << h_nClusters->Integral() << endl;

}

// use the classes name to register the physics class inside Ant
// this is black macro magic what's used here...but it works :)

AUTO_REGISTER_PHYSICS(scratch_damaurer_Tutorial_kinFit)
