#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"

namespace ant {
namespace analysis {
namespace physics {

class scratch_kosi_pg_ppi0g_pgg : public Physics {
public:
    scratch_kosi_pg_ppi0g_pgg (const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
    virtual void Finish() override;

protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
private:
    TH1D* h_nClusters;
    TH1D* h_TaggerTime;
    TH1D* h_ThetaDist;
    TH1D* h_NumCands;
    TH1D* h_ThetaPi0;
    TH1D* h_ThetaPi0Zoom;

    double vetoEthreshold = 0.1;
    static const int num_gammas = 2;
    static const int num_protons = 1;

    unsigned int num_3decay;

    utils::TriggerSimulation TrigSim;
    PromptRandom::Switch promptrandom;
};
}}}
