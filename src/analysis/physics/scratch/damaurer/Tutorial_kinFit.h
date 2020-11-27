// better than using the standard ifndef include guards...
#pragma once

// physics classes need to derive from this interface
#include "physics/Physics.h"
#include "TLorentzVector.h"

// Ant provides many utility classes,
// such as support for prompt-random handling,
// or trigger simulation
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"

// handle ROOT's TTree
#include "base/WrapTTree.h"

//Including detecor properties
#include "expconfig/detectors/Tagger.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/TAPSVeto.h"
#include "expconfig/detectors/CB.h"

//Include Kinematic fitter
#include "analysis/utils/fitter/KinFitter.h"
#include "analysis/utils/Uncertainties.h"
#include "analysis/utils/ProtonPhotonCombs.h"

// the physics classes reside in this nested namespace
namespace ant {
namespace analysis {
namespace physics {

// choose a nice name for your class
class scratch_damaurer_Tutorial_kinFit : public Physics {
public:

    // use "struct" to avoid all those public keywords...
    struct tree_t : WrapTTree {
        // define two branches with ADD_BRANCH_T (again, this is black macro magic)
        ADD_BRANCH_T(double,   TaggW)
        ADD_BRANCH_T(unsigned, nClusters)
        // you may even use more complex types with non-default ctor, for example
//        ADD_BRANCH_T(std::vector<double>, Numbers, 3) // vector Numbers is three items large now
    };

protected:

    // use this instance to handle prompt-random weighting
    PromptRandom::Switch promptrandom;
    // this encapsulates trigger timings and other related things (CBESum)
    utils::TriggerSimulation triggersimu;

    using particle_comb_t = utils::ProtonPhotonCombs::comb_t;
    using particle_combs_t = utils::ProtonPhotonCombs::Combinations_t;

    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

    utils::UncertaintyModelPtr fit_model;
    utils::KinFitter fitter;

    std::shared_ptr<expconfig::detector::Tagger> tagger_detector;
    std::shared_ptr<expconfig::detector::CB> cb_detector;
    std::shared_ptr<expconfig::detector::PID> pid_detector;
    std::shared_ptr<expconfig::detector::TAPSVeto> veto_detector;


private:

    TH1D* h_nClusters;   
    TH1D* h_TaggerTime;
    //TH1D* hist;

    tree_t t;

public:
    // physics need to implement this public constructor...
    scratch_damaurer_Tutorial_kinFit(const std::string& name, OptionsPtr opts);

    // ...and the following method (use override keyword)
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    // optional method (called in interactive mode to show histograms)
    virtual void ShowResult() override;
    virtual void Finish() override;

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);
};


}}} // end of ant::analysis::physics
