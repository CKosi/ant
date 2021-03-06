#pragma once

#include "analysis/physics/Physics.h"

#include <string>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class ProtonCheck : public Physics {
protected:

    TH2D* tof;
    TH2D* tof_trueE;
    TH2D* dEE;

    TH2D* e_recov;

    TH1D* cand_mult;

    TH1D* theta;
    TH2D* theta_corr;

public:
    ProtonCheck(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
