#include "Tutorial.h"

#include "base/Logger.h"

#include "plot/root_draw.h"

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Tutorial::Tutorial(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    BinSettings bins_nClusters(20);

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    h_nClusters = HistFac.makeTH1D("Number of Clusters", // title
                                   "nClusters","#",      // xlabel, ylabel
                                   bins_nClusters,       // our binnings, may write directly BinSettings(10) here
                                   "h_nClusters"         // ROOT object name, auto-generated if omitted
                                   );

    h_nClusters_pr = HistFac.makeTH1D("Number of Clusters - prompt-random",
                                      "nClusters","#",
                                      bins_nClusters,
                                      "h_nClusters_pr"
                                      );

    // define some prompt and random windows (in nanoseconds)
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});
}

void Tutorial::ProcessEvent(const TEvent& event, manager_t&)
{
    for(auto& taggerhit : event.Reconstructed().TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        h_nClusters_pr->Fill(event.Reconstructed().Clusters.size(), promptrandom.FillWeight());
    }
    h_nClusters->Fill(event.Reconstructed().Clusters.size());
}

void Tutorial::ShowResult()
{
    // ShowResult is called after processing of events has finished,
    // and interactive mode (aka non-batchmode) is chosen

    // ant::canvas nice wrapper around TCanvas
    ant::canvas(GetName()+": Basic plots")
            << h_nClusters
            << h_nClusters_pr
            << endc; // actually draws the canvas
}

// use the classes name to register the physics class inside Ant
// this is black macro magic what's used here...but it works :)
AUTO_REGISTER_PHYSICS(Tutorial)