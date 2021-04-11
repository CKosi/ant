#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"
#include "utils/fitter/KinFitter.h"
#include "analysis/utils/Uncertainties.h"

namespace ant {
namespace analysis {
namespace physics {

class scratch_kosi_pg_ppi0_pgg : public Physics {
public:
    scratch_kosi_pg_ppi0_pgg (const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;

protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    utils::UncertaintyModelPtr fit_model;
    utils::KinFitter fitter;

private:
    TH1D* h_nClusters;
    TH1D* h_TaggerTime;
    TH1D* h_NumCands;
    TH1D* h_ThetaPi0;
    TH1D* h_ThetaPi0Cut;
    TH1D* h_GammaIM;
    TH1D* h_MM;
    TH2D* h_cs;

    double vetoEthreshold = 0.1;

    utils::TriggerSimulation TrigSim;
    PromptRandom::Switch promptrandom;

};
}}}
