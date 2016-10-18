#pragma once

#include <vector>
#include <list>
#include <type_traits>
#include <random>

#include "analysis/physics/Physics.h"
#include "analysis/utils/Fitter.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/plot/PromptRandomHist.h"
#include "base/WrapTTree.h"

#include "root-addons/cbtaps_display/TH2CB.h"

#include "TLorentzVector.h"

namespace ant {
namespace analysis {
namespace physics {

class EtapDalitz : public Physics {

public:
    struct Tree_t : WrapTTree {
        Tree_t();

        ADD_BRANCH_T(unsigned,                    nCands)

        ADD_BRANCH_T(std::vector<TLorentzVector>, photons, 3)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons_kinfitted, 3)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons_kinfit_freeZ, 3)
        ADD_BRANCH_T(std::vector<TLorentzVector>, photons_treefitted, 3)
        ADD_BRANCH_T(std::vector<double>,         photons_Time, 3)
        ADD_BRANCH_T(std::vector<TVector2>,       photons_PSA, 3)
        ADD_BRANCH_T(std::vector<int>,            photons_detector, 3)
        ADD_BRANCH_T(std::vector<int>,            photons_clusterSize, 3)
        ADD_BRANCH_T(std::vector<int>,            photons_centralElem, 3)
        ADD_BRANCH_T(std::vector<double>,         photons_vetoE, 3)
        ADD_BRANCH_T(std::vector<int>,            photons_vetoChannel, 3)

        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_E_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_theta_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_phi_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_freeZ_E_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_freeZ_theta_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_kinfit_freeZ_phi_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_treefit_E_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_treefit_theta_pulls, 3)
        ADD_BRANCH_T(std::vector<double>,         photon_treefit_phi_pulls, 3)

        ADD_BRANCH_T(TLorentzVector,              p)
        ADD_BRANCH_T(TLorentzVector,              p_kinfitted)
        ADD_BRANCH_T(TLorentzVector,              p_kinfit_freeZ)
        ADD_BRANCH_T(TLorentzVector,              p_treefitted)
        ADD_BRANCH_T(double,                      p_Time)
        ADD_BRANCH_T(TVector2,                    p_PSA)
        ADD_BRANCH_T(int,                         p_detector)
        ADD_BRANCH_T(int,                         p_clusterSize, 3)
        ADD_BRANCH_T(int,                         p_centralElem, 3)
        ADD_BRANCH_T(double,                      p_vetoE)
        ADD_BRANCH_T(int,                         p_vetoChannel, 3)

        ADD_BRANCH_T(double,                      p_kinfit_theta_pull)
        ADD_BRANCH_T(double,                      p_kinfit_phi_pull)
        ADD_BRANCH_T(double,                      p_kinfit_freeZ_theta_pull)
        ADD_BRANCH_T(double,                      p_kinfit_freeZ_phi_pull)
        ADD_BRANCH_T(double,                      p_treefit_theta_pull)
        ADD_BRANCH_T(double,                      p_treefit_phi_pull)

        ADD_BRANCH_T(double,                      TaggW)
        ADD_BRANCH_T(double,                      TaggW_wide)
        ADD_BRANCH_T(double,                      TaggE)
        ADD_BRANCH_T(double,                      TaggT)
        ADD_BRANCH_T(unsigned,                    TaggCh)

        ADD_BRANCH_T(double,                      beam_E_kinfitted)
        ADD_BRANCH_T(double,                      beam_E_kinfit_freeZ)
        ADD_BRANCH_T(double,                      beam_E_treefitted)
        ADD_BRANCH_T(double,                      beam_kinfit_E_pull)
        ADD_BRANCH_T(double,                      beam_kinfit_freeZ_E_pull)
        ADD_BRANCH_T(double,                      beam_treefit_E_pull)
        ADD_BRANCH_T(double,                      kinfit_ZVertex)
        ADD_BRANCH_T(double,                      kinfit_ZVertex_pull)
        ADD_BRANCH_T(double,                      kinfit_freeZ_ZVertex)
        ADD_BRANCH_T(double,                      kinfit_freeZ_ZVertex_pull)
        ADD_BRANCH_T(double,                      treefit_ZVertex)
        ADD_BRANCH_T(double,                      treefit_ZVertex_pull)

        ADD_BRANCH_T(double,                      kinfit_chi2)
        ADD_BRANCH_T(double,                      kinfit_probability)
        ADD_BRANCH_T(unsigned,                    kinfit_iterations)
        ADD_BRANCH_T(unsigned,                    kinfit_DoF)
        ADD_BRANCH_T(double,                      kinfit_freeZ_chi2)
        ADD_BRANCH_T(double,                      kinfit_freeZ_probability)
        ADD_BRANCH_T(unsigned,                    kinfit_freeZ_iterations)
        ADD_BRANCH_T(unsigned,                    kinfit_freeZ_DoF)
        ADD_BRANCH_T(double,                      treefit_chi2)
        ADD_BRANCH_T(double,                      treefit_probability)
        ADD_BRANCH_T(unsigned,                    treefit_iterations)
        ADD_BRANCH_T(unsigned,                    treefit_DoF)

        ADD_BRANCH_T(double,                      CBSumE)
        ADD_BRANCH_T(double,                      CBAvgTime)

        ADD_BRANCH_T(unsigned,                    channel)
        ADD_BRANCH_T(bool,                        MCtrue)
        ADD_BRANCH_T(double,                      trueZVertex)

        ADD_BRANCH_T(TLorentzVector,              etap)
        ADD_BRANCH_T(TLorentzVector,              etap_kinfit)
        ADD_BRANCH_T(TLorentzVector,              etap_kinfit_freeZ)
        ADD_BRANCH_T(TLorentzVector,              etap_treefit)
        ADD_BRANCH_T(TLorentzVector,              mm)
        ADD_BRANCH_T(double,                      copl)
    };

protected:
    TH1D* h_tagger_time = nullptr;
    TH1D* h_tagger_time_CBavg = nullptr;

    TH1D* h_pTheta = nullptr;
    TH1D* h_protonVeto = nullptr;
    TH1D* h_etapIM_final = nullptr;
    TH2D* h_IM2d = nullptr;
    TH2* h_etap = nullptr;
    TH2* h_proton = nullptr;

    TH1D* h_counts = nullptr;
    TH1D* h_nCands = nullptr;
    TH1D* h_cluster_CB = nullptr;
    TH1D* h_cluster_TAPS = nullptr;
    TH1D* missed_channels = nullptr;
    TH1D* found_channels  = nullptr;

    static constexpr unsigned N_FINAL_STATE = 4;
    static constexpr double ETAP_IM = 957.78;
    static constexpr double ETAP_SIGMA = 50.;
    // threshold to check if double value should be treated as zero
    static constexpr double EPSILON = 2*std::numeric_limits<double>::epsilon();
    // threshold for cluster energies
    static constexpr double CLUSTER_TRESH = 25.;
    // cuts
    static constexpr bool Q2_PRESELECTION = false;
    static constexpr bool PROBABILITY_CUT = false;
    static constexpr double PROBABILITY = .02;
    static constexpr bool ANTI_PI0_CUT = false;
    static constexpr double ANTI_PI0_LOW = 102.;
    static constexpr double ANTI_PI0_HIGH = 170.;
    // which fit should be used to determine best candidate combination?
    static constexpr bool USE_TREEFIT = false;

    struct PerChannel_t {
        std::string title;
        std::string name;
        TH1D* steps = nullptr;
        TH1D* etapIM = nullptr;
        TH1D* etapIM_kinfit = nullptr;
        TH1D* etapIM_kinfit_freeZ = nullptr;
        TH1D* etapIM_treefit = nullptr;
        TH1D* etapIM_cand = nullptr;
        TH1D* etapIM_final = nullptr;
        TH2D* IM2d = nullptr;
        TH1D* MM = nullptr;
        TH1D* hCopl = nullptr;
        TH1D* hCopl_final = nullptr;
        TH1D* treefitChi2 = nullptr;
        TH1D* treefitProb = nullptr;
        TH1D* treefitIter = nullptr;
        TH1D* kinfitChi2 = nullptr;
        TH1D* kinfitProb = nullptr;
        TH1D* kinfitIter = nullptr;
        TH1D* kinfit_freeZ_chi2 = nullptr;
        TH1D* kinfit_freeZ_prob = nullptr;
        TH1D* kinfit_freeZ_iter = nullptr;
        TH1D* kinfit_ZVertex = nullptr;
        TH1D* kinfit_freeZ_ZVertex = nullptr;
        TH1D* treefit_ZVertex = nullptr;
        TH1D* effect_rad = nullptr;
        TH2D* effect_rad_E = nullptr;
        TH1D* cluster_size = nullptr;
        TH2D* cluster_size_E = nullptr;

        TH2* proton_E_theta = nullptr;

        PerChannel_t(const std::string& Name, const std::string& Title, HistogramFactory& hf);

        void Show();
        void Fill(const TEventData& d);
    };

    std::map<std::string, PerChannel_t> channels;
    std::map<std::string, HistogramFactory&> productions;

    Tree_t t;
    PromptRandom::Switch promptrandom;
    using uncertainty_model_t = utils::UncertaintyModels::Optimized_Oli1;
    utils::UncertaintyModelPtr model;
    utils::KinFitter kinfit;
    utils::KinFitter kinfit_freeZ;
    utils::TreeFitter treefitter_etap;

    std::shared_ptr<ant::Detector_t> cb;

    template<typename T>
    void shift_right(std::vector<T>&);

    template <typename iter>
    LorentzVec sumlv(iter start, iter end);

    void remove_char(std::string&, char);
    void remove_chars(std::string&, std::initializer_list<char>);

    double calc_effective_radius(const TCandidatePtr);

    ParticleTypeTree base_tree();
    ParticleTypeTree etap_3g();

    void count_clusters(const TCandidateList&);
    bool q2_preselection(const TEventData&, const double) const;

public:

    EtapDalitz(const std::string& name, OptionsPtr opts);

    static APLCON::Fit_Settings_t MakeFitSettings(unsigned);

    bool doFit_checkProb(const TTaggerHit& taggerhit,
                         const TParticlePtr proton,
                         const TParticleList photons,
                         PerChannel_t& h,
                         Tree_t& t,
                         double& best_prob_fit);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

    using decaytree_t = ant::Tree<const ParticleTypeDatabase::Type&>;

    struct ReactionChannel_t {
        std::string name = "";
        std::shared_ptr<decaytree_t> tree = nullptr;
        int color = kBlack;

        ReactionChannel_t() = default;
        ReactionChannel_t(const std::string& n);
        ReactionChannel_t(const std::shared_ptr<decaytree_t>& t, const int c);
        ReactionChannel_t(const std::shared_ptr<decaytree_t>& t, const std::string& n, const int c);
        ~ReactionChannel_t();
    };

    struct ReactionChannelList_t {
        static const unsigned other_index;
        std::map<unsigned, ReactionChannel_t> channels;
        unsigned identify(const TParticleTree_t &tree) const;
    };

    static ReactionChannelList_t makeChannels();
    static const ReactionChannelList_t reaction_channels;
};

}}} // namespace ant::analysis::physics