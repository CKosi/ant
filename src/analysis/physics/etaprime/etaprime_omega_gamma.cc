#include "etaprime_omega_gamma.h"

#include "plot/root_draw.h"
#include "plot/TH1Dcut.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/combinatorics.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"

#include <TTree.h>

#include <limits>


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

const ParticleTypeTree EtapOmegaG::ptreeSignal = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g);
const ParticleTypeTree EtapOmegaG::ptreeReference = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g);


EtapOmegaG::EtapOmegaG(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    kinfitter_2("kinfitter_2",2),
    kinfitter_4("kinfitter_4",4)
{
    const interval<double> prompt_range{-2.5,1.5};
    promptrandom.AddPromptRange(prompt_range); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    promptrandom_tight.AddPromptRange(prompt_range); // slight offset due to CBAvgTime reference
    promptrandom_tight.AddRandomRange({-30,-10});  // just ensure to be way off prompt peak
    promptrandom_tight.AddRandomRange({  10,30});


    h_CommonCuts = HistFac.makeTH1D("Common Cuts", "", "#", BinSettings(15),"h_TotalEvents");
    h_MissedBkg = HistFac.makeTH1D("Missed Background", "", "#", BinSettings(25),"h_MissedBkg");


    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw runtime_error("EtapOmegaG needs a setup");
    kinfitter_2.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
    kinfitter_4.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

    t.CreateBranches(HistFac.makeTTree("treeCommon"));

    Sig.t.CreateBranches(HistFac.makeTTree("treeSig"));
    SigFitted.t.CreateBranches(HistFac.makeTTree("treeSigFitted"));
    Ref.t.CreateBranches(HistFac.makeTTree("treeRef"));
    RefFitted.t.CreateBranches(HistFac.makeTTree("treeRefFitted"));

}

void EtapOmegaG::ProcessEvent(const TEvent& event, manager_t&)
{
    // we start with some general candidate handling,
    // later we split into ref/sig analysis according to
    // number of photons

    TEventData& data = *event.Reconstructed;

    h_CommonCuts->Fill("Seen",1.0);

    auto& particletree = event.MCTrue->ParticleTree;

    h_CommonCuts->Fill("MCTrue #eta'", 0); // ensure the bin is there...
    if(particletree) {
        // note: this might also match to g p -> eta' eta' p,
        // but this is kinematically forbidden
        if(utils::ParticleTools::FindParticle(ParticleTypeDatabase::EtaPrime, particletree, 1)) {
            h_CommonCuts->Fill("MCTrue #eta'", 1);
        }
    }

    if(data.Trigger.CBEnergySum<=550)
        return;
    h_CommonCuts->Fill("CBEnergySum>550",1.0);
    t.CBSumE = data.Trigger.CBEnergySum;

    t.CBAvgTime = event.Reconstructed->Trigger.CBTiming;
    if(!isfinite(t.CBAvgTime))
        return;
    h_CommonCuts->Fill("CBAvgTime ok",1.0);

    // use struct to gather particles
    Particles_t particles;

    // identify the proton here as slowest cluster in TAPS
    // using CBAvgTime does not help here, since it's constant
    // over each TAPS clusters
    /// \todo think about using beta here as in EtapProton?
    t.ProtonTime = std_ext::NaN;
    for(const TCandidatePtr& cand : data.Candidates) {
        if(cand->Detector & Detector_t::Type_t::TAPS) {
            if(!isfinite(t.ProtonTime) || t.ProtonTime < cand->Time) {
                t.ProtonTime = cand->Time;
                particles.Proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand);
            }
        }
    }

    if(!particles.Proton)
        return;

    h_CommonCuts->Fill("p in TAPS", 1.0);

    // remaining candidates are photons
    const auto nPhotons = data.Candidates.size() - 1;
    if(nPhotons != 2 && nPhotons != 4)
        return;
    h_CommonCuts->Fill("nPhotons==2|4", 1.0);

    particles.PhotonSum.SetPxPyPzE(0,0,0,0);
    t.nPhotonsCB = 0;
    t.nPhotonsTAPS = 0;
    t.CBSumVetoE = 0;
    for(const TCandidatePtr& cand : data.Candidates) {
         if(cand == particles.Proton->Candidate)
             continue;
         if(cand->Detector & Detector_t::Type_t::CB) {
             t.nPhotonsCB++;
             t.CBSumVetoE += cand->VetoEnergy;
         }
         if(cand->Detector & Detector_t::Type_t::TAPS)
             t.nPhotonsTAPS++;
         auto photon = make_shared<TParticle>(ParticleTypeDatabase::Photon, cand);
         particles.PhotonSum += *photon;
         particles.Photons.emplace_back(move(photon));
    }
    assert(particles.Photons.size() == nPhotons);
    assert(nPhotons == t.nPhotonsCB + t.nPhotonsTAPS);

    // don't bother with events where proton coplanarity is not ok
    // we use some rather large window here...
    t.ProtonCopl = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(particles.Proton->Phi() - particles.PhotonSum.Phi() - M_PI ));
    const interval<double> ProtonCopl_cut(-30, 30);
    if(!ProtonCopl_cut.Contains(t.ProtonCopl))
        return;
    h_CommonCuts->Fill("ProtonCopl ok", 1.0);

    // sum up the PID energy
    // (might be different to matched CB/PID Veto energy)
    t.PIDSumE = 0;
    for(const TClusterPtr& cl : data.Clusters) {
        if(cl->DetectorType == Detector_t::Type_t::PID) {
            t.PIDSumE += cl->Energy;
        }
    }



    // select signal or reference according to nPhotons

    t.IsSignal = nPhotons == 4; // else nPhotons==2, so reference
    utils::KinFitter& fitter = t.IsSignal ? kinfitter_4 : kinfitter_2;



    // do some MCTrue identification (if available)
    t.MCTrue = 0; // indicate data by default
    TParticleTree_t ptree_sigref = nullptr; // used by Sig_t::Process to check matching
    if(particletree) {
        auto get_mctrue = [] (ParticleTypeTree expected, TParticleTree_t particletree) -> unsigned {
            if(particletree->IsEqual(expected, utils::ParticleTools::MatchByParticleName))
                return 1;
            unsigned n = 10;
            for(const auto& ptreeBkg : ptreeBackgrounds) {
                if(particletree->IsEqual(ptreeBkg.Tree, utils::ParticleTools::MatchByParticleName))
                    return n;
                n++;
            }
            return 9; // missed background
        };

        // 1=Signal, 2=Reference, 9=MissedBkg, >=10 found in ptreeBackgrounds
        if(t.IsSignal) {
            t.MCTrue = get_mctrue(ptreeSignal, particletree);
            if(t.MCTrue==1)
                ptree_sigref = particletree;
        }
        else {
            t.MCTrue = get_mctrue(ptreeReference, particletree);
            if(t.MCTrue==1) {
                ptree_sigref = particletree;
                t.MCTrue++;
            }
        }

        if(t.MCTrue==9) {
            const auto& decaystr = utils::ParticleTools::GetDecayString(particletree);
            h_MissedBkg->Fill(decaystr.c_str(), 1.0);
        }
    }

    Sig.ResetBranches();
    Ref.ResetBranches();
    if(t.IsSignal)
        Sig.Process(particles, ptree_sigref);
    else
        Ref.Process(particles);

    // loop over tagger hits, do KinFit
    bool kinfit_ok = false;
    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        promptrandom_tight.SetTaggerHit(taggerhit.Time - t.CBAvgTime);
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // missing mass
        const TLorentzVector beam_target = taggerhit.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        t.MissingMass = (beam_target - particles.PhotonSum).M();

        t.TaggW = promptrandom.FillWeight();
        t.TaggW_tight = promptrandom_tight.FillWeight();
        t.TaggE = taggerhit.PhotonEnergy;
        t.TaggT = taggerhit.Time;
        t.TaggCh = taggerhit.Channel;

        // do kinfit
        fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
        fitter.SetProton(particles.Proton);
        fitter.SetPhotons(particles.Photons);
        const auto& fit_result = fitter.DoFit();

        SigFitted.ResetBranches();
        RefFitted.ResetBranches();

        t.KinFitChi2 = std_ext::NaN;
        t.KinFitIterations = 0;
        if(fit_result.Status == APLCON::Result_Status_t::Success) {
            kinfit_ok = true;
            t.KinFitChi2 = fit_result.ChiSquare;
            t.KinFitIterations = fit_result.NIterations;
            Particles_t fitted_particles;
            fitted_particles.Proton = fitter.GetFittedProton();
            fitted_particles.Photons = fitter.GetFittedPhotons();
            fitted_particles.PhotonSum.SetPxPyPzE(0,0,0,0);

            for(const TParticlePtr& p : fitted_particles.Photons)
                fitted_particles.PhotonSum += *p;

            if(t.IsSignal)
                SigFitted.Process(fitted_particles, ptree_sigref);
            else
                RefFitted.Process(fitted_particles);
        }

        t.Tree->Fill();
        Sig.t.Tree->Fill();
        SigFitted.t.Tree->Fill();
        Ref.t.Tree->Fill();
        RefFitted.t.Tree->Fill();

    }

    if(kinfit_ok)
        h_CommonCuts->Fill("KinFit OK", 1.0);
    if(nPhotons==2)
        h_CommonCuts->Fill("nPhotons==2", 1.0);
    if(nPhotons==4)
        h_CommonCuts->Fill("nPhotons==4", 1.0);

}



EtapOmegaG::Sig_t::Sig_t() :
    treefitter("sig_treefitter",
               utils::ParticleTools::GetProducedParticle(EtapOmegaG::ptreeSignal)),
    fitted_EtaPrime(treefitter.GetTreeNode(ParticleTypeDatabase::EtaPrime)),
    fitted_Omega(treefitter.GetTreeNode(ParticleTypeDatabase::Omega)),
    fitted_Pi0(treefitter.GetTreeNode(ParticleTypeDatabase::Pi0))
{
    t.ggg().resize(4);
    t.gg_gg1().resize(3);
    t.gg_gg2().resize(3);

    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    if(!setup)
        throw runtime_error("EtapOmegaG needs a setup");
    treefitter.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");

    // search dependent gammas and remember the tree nodes in the fitter

    auto find_photons = [] (utils::TreeFitter::tree_t fitted) {
        std::vector<utils::TreeFitter::tree_t> photons;
        for(const auto& d : fitted->Daughters())
            if(d->Get().TypeTree->Get() == ParticleTypeDatabase::Photon)
                photons.emplace_back(d);
        return photons;
    };

    fitted_g_EtaPrime = find_photons(fitted_EtaPrime).at(0);
    fitted_g_Omega = find_photons(fitted_Omega).at(0);
    fitted_g1_Pi0 = find_photons(fitted_Pi0).at(0);
    fitted_g2_Pi0 = find_photons(fitted_Pi0).at(1);
}

void EtapOmegaG::Sig_t::ResetBranches()
{
    t.TreeFitChi2 = std_ext::NaN;
    t.TreeFitIterations = 0;

    std::fill(t.ggg().begin(), t.ggg().end(), std_ext::NaN);
    std::fill(t.gg_gg1().begin(), t.gg_gg1().end(), std_ext::NaN);
    std::fill(t.gg_gg2().begin(), t.gg_gg2().end(), std_ext::NaN);

    t.IM_EtaPrime_fitted  = std_ext::NaN;
    t.IM_Omega_fitted = std_ext::NaN;
    t.IM_Pi0_fitted = std_ext::NaN;

    t.IM_EtaPrime_best = std_ext::NaN;
    t.IM_Omega_best = std_ext::NaN;
    t.IM_Pi0_best = std_ext::NaN;

    t.Bachelor_best_best = std_ext::NaN;
    t.Bachelor_best_fit = std_ext::NaN;
    t.Bachelor_fit_best = std_ext::NaN;
    t.Bachelor_fit_fit = std_ext::NaN;

    t.MCTrueMatch = 0;
}

void EtapOmegaG::Sig_t::Process(const EtapOmegaG::Particles_t& particles, TParticleTree_t ptree_sigref)
{
    ResetBranches();

    assert(particles.Photons.size() == 4);

    //  ggg combinatorics
    auto it_ggg = t.ggg().begin();
    for( auto comb = utils::makeCombination(particles.Photons,3); !comb.Done(); ++comb ) {
        *it_ggg = (*comb.at(0) + *comb.at(1) + *comb.at(2)).M();
        ++it_ggg;
    }

    // gg/gg "Goldhaber" combinatorics
    const auto goldhaber_comb = vector<vector<unsigned>>({{0,1,2,3},{0,2,1,3},{0,3,1,2}});
    auto it_gg_gg1 = t.gg_gg1().begin();
    auto it_gg_gg2 = t.gg_gg2().begin();
    for(auto i : goldhaber_comb) {
        const auto& p = particles.Photons;
        *it_gg_gg1 = (*p[i[0]] + *p[i[1]]).M();
        *it_gg_gg2 = (*p[i[2]] + *p[i[3]]).M();
        ++it_gg_gg1;
        ++it_gg_gg2;
    }


    // g_Omega to check against MCTrue
    TParticlePtr g_Omega_best;
    // the EtaPrime bachelor photon is most important to us...
    TParticlePtr g_EtaPrime_best;
    TLorentzVector g_EtaPrime_fitted;
    const TLorentzVector& EtaPrime_best = particles.PhotonSum; // does not change with permutation
    t.IM_EtaPrime_best = EtaPrime_best.M();
    TLorentzVector EtaPrime_fitted;

    treefitter.SetLeaves(particles.Photons);
    APLCON::Result_t r;

    while(treefitter.NextFit(r)) {
        if(r.Status != APLCON::Result_Status_t::Success)
            continue;
        if(isfinite(t.TreeFitChi2) && r.ChiSquare>t.TreeFitChi2)
            continue;
        // found fit with better chi2
        t.TreeFitChi2 = r.ChiSquare;
        t.TreeFitIterations = r.NIterations;

        EtaPrime_fitted = fitted_EtaPrime->Get().Particle;

        // IM fitted expected to be delta peaks since they were fitted...
        t.IM_EtaPrime_fitted = EtaPrime_fitted.M();
        t.IM_Omega_fitted = fitted_Omega->Get().Particle.M();
        t.IM_Pi0_fitted = fitted_Pi0->Get().Particle.M();

        // have a look at the assigned gammas to Pi0/Omega
        const TLorentzVector& Pi0_best = *fitted_g1_Pi0->Get().SetParticle + *fitted_g2_Pi0->Get().SetParticle;
        t.IM_Pi0_best = Pi0_best.M();

        g_Omega_best = fitted_g_Omega->Get().SetParticle;
        const TLorentzVector& Omega_best = *g_Omega_best + Pi0_best;
        t.IM_Omega_best = Omega_best.M();

        // have a look at the EtaPrime bachelor photon
        g_EtaPrime_fitted = fitted_g_EtaPrime->Get().Particle;
        g_EtaPrime_best = fitted_g_EtaPrime->Get().SetParticle;
    }

    if(isfinite(t.TreeFitChi2)) {

        // do that Bachelor photon boosting
        auto do_boost = [] (const TLorentzVector& bachelor, const TLorentzVector& etaprime) {
            TLorentzVector boosted(bachelor);
            boosted.Boost(-etaprime.BoostVector());
            return boosted;
        };

        t.Bachelor_best_best = do_boost(*g_EtaPrime_best, EtaPrime_best).E();
        t.Bachelor_best_fit  = do_boost(*g_EtaPrime_best, EtaPrime_fitted).E();
        t.Bachelor_fit_best  = do_boost(g_EtaPrime_fitted, EtaPrime_best).E();
        t.Bachelor_fit_fit   = do_boost(g_EtaPrime_fitted, EtaPrime_fitted).E();


        // check MC matching
        if(ptree_sigref) {

            auto true_photons = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, ptree_sigref);
            assert(true_photons.size() == 4);
            auto match_bycandidate = [] (const TParticlePtr& mctrue, const TParticlePtr& recon) {
                return mctrue->Angle(*recon->Candidate); // TCandidate converts into TVector3
            };
            auto matched = utils::match1to1(true_photons, particles.Photons,
                                            match_bycandidate,IntervalD(0.0, 15.0*TMath::DegToRad()));
            if(matched.size() == 4) {
                // do that tedious photon determination (rewriting the matcher somehow would be nice....)
                auto select_daughter = [] (TParticleTree_t tree, const ParticleTypeDatabase::Type& type) {
                    auto d = tree->Daughters().front()->Get()->Type() == type ?
                                 tree->Daughters().front() : tree->Daughters().back();
                    assert(d->Get()->Type() == type);
                    return d;
                };

                auto etap = select_daughter(ptree_sigref, ParticleTypeDatabase::EtaPrime);
                auto g_Etap = select_daughter(etap, ParticleTypeDatabase::Photon);
                auto omega = select_daughter(etap, ParticleTypeDatabase::Omega);
                auto g_Omega = select_daughter(omega, ParticleTypeDatabase::Photon);

                auto g_Etap_matched = utils::FindMatched(matched, g_Etap->Get());
                auto g_Omega_matched = utils::FindMatched(matched, g_Omega->Get());
                if(g_Etap_matched == g_EtaPrime_best)
                    t.MCTrueMatch += 1;
                if(g_Omega_matched == g_Omega_best)
                    t.MCTrueMatch += 2;

            }
        }
    }


}

void EtapOmegaG::Ref_t::ResetBranches()
{
    t.IM_2g = std_ext::NaN;
}

void EtapOmegaG::Ref_t::Process(const EtapOmegaG::Particles_t& particles)
{
    assert(particles.Photons.size() == 2);
    t.IM_2g = particles.PhotonSum.M();
}

void EtapOmegaG::ShowResult()
{
    canvas("Overview") << h_CommonCuts << h_MissedBkg << endc;
}


const std::vector<EtapOmegaG::Background_t> EtapOmegaG::ptreeBackgrounds = {
    {"1Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g)},
    {"2Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)},
    {"Pi0Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)},
    {"3Pi0", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g)},
    {"OmegaPi0g", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g)},
    {"OmegaPi0PiPPiM", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_Pi0PiPPiM_2g)},
    {"EtaP2Pi0Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g)},
    {"2Pi0Dalitz", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_2ggEpEm)},
    {"3Pi0Dalitz", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_4ggEpEm)},
    {"1Eta", ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Eta_2g)},
};

EtapOmegaG_MC::EtapOmegaG_MC(const std::string& name, OptionsPtr opts) : Physics(name, opts),
    sig_HistFac("Sig",HistFac),
    ref_HistFac("Ref",HistFac),
    sig_hists(sig_HistFac),
    ref_hists(ref_HistFac),
    sig_TTree(sig_HistFac.makeTTree("tree")),
    ref_TTree(ref_HistFac.makeTTree("tree"))
{
    h_TotalEvents = HistFac.makeTH1D("Total Events", "", "#", BinSettings(5),"h_TotalEvents");
    h_TotalEvents->Fill("Total",0);
    h_TotalEvents->Fill("#eta'", 0);
    h_TotalEvents->Fill("Signal",0);
    h_TotalEvents->Fill("Reference",0);

    sig_TTree.SetBranches();
    ref_TTree.SetBranches();

    treeSignal = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g);
    treeReference = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g);

    sig_perDecayHists.emplace_back(
                "Signal",
                sig_HistFac,
                treeSignal
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_2pi0",
                sig_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_pi0eta",
                sig_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_3pi0",
                sig_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_6g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_OmegaPi0g",
                sig_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gPi0_3g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_1pi0",
                sig_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_Etap_Eta2Pi0",
                sig_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_2pi0_1Dalitz",
                sig_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_2ggEpEm)
                );
    sig_perDecayHists.emplace_back(
                "Sig_Bkg_3pi0_1Dalitz",
                sig_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::ThreePi0_4ggEpEm)
                );
    sig_perDecayHists.emplace_back("Signal_Bkg_Other", sig_HistFac);


    ref_perDecayHists.emplace_back(
                "Reference",
                ref_HistFac,
                treeReference
                );
    ref_perDecayHists.emplace_back(
                "Ref_Bkg_2pi0",
                ref_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g)
                );
    ref_perDecayHists.emplace_back(
                "Ref_Bkg_pi0eta",
                ref_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g)
                );
    ref_perDecayHists.emplace_back(
                "Ref_Bkg_1pi0",
                ref_HistFac,
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g)
                );
    ref_perDecayHists.emplace_back("Ref_Bkg_Other", ref_HistFac);
}

EtapOmegaG_MC::histogram_t::histogram_t(SmartHistFactory HistFac)
{
    Steps = HistFac.makeTH1D("Steps", "", "#", BinSettings(10),"steps");
    MissedBkg = HistFac.makeTH1D("Missed Bkg channels", "", "#", BinSettings(20),"missed_bkg");
}

EtapOmegaG_MC::sig_perDecayHists_t::sig_perDecayHists_t(SmartHistFactory HistFac)
{
    Steps = HistFac.makeTH1D("Steps", "", "#", BinSettings(10),"steps");

    const BinSettings bins_im(1200);

    gggg = HistFac.makeTH1D("4#gamma IM","4#gamma IM / MeV","events",bins_im,"gggg");
    ggg = HistFac.makeTH1D("3#gamma IM","3#gamma IM / MeV","events",bins_im,"ggg");
    gg = HistFac.makeTH1D("2#gamma IM","2#gamma IM / MeV","events",bins_im,"gg");

    const BinSettings bins_goldhaber(400, 0, 800);
    const string axislabel_goldhaber("2#gamma IM / MeV");

    IM_gg_gg = HistFac.makeTH2D("IM 2#gamma vs. 2#gamma (Goldhaber plot)",
                                axislabel_goldhaber, axislabel_goldhaber,
                                bins_goldhaber, bins_goldhaber,
                                "IM_gg_gg"
                                );
    IM_gg_gg_cut = HistFac.makeTH2D("IM 2#gamma vs. 2#gamma after cut",
                                axislabel_goldhaber, axislabel_goldhaber,
                                bins_goldhaber, bins_goldhaber,
                                "IM_gg_gg_cut"
                                );

    Proton_Copl = HistFac.makeTH1D<TH1Dcut>("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"Proton_Coplanarity");

    IM_etap_omega = HistFac.makeTH2D("#eta' vs. #omega IM",
                               "#eta' IM / MeV",
                               "#omega IM / MeV",
                               BinSettings(400, 600, 1200),
                               BinSettings(400, 400, 1000),
                               "IM_etap_omega");
    IM_pi0 = HistFac.makeTH1D("#pi^{0}","#pi^{0} IM / MeV","events",BinSettings(300, 0, 400),"IM_pi0");

    BinSettings bins_mm(300, 600, 1300);
    MM_gggg = HistFac.makeTH1D("M_{miss}","M_{miss} / MeV","events",bins_mm,"MM_gggg");
    MM_etap = HistFac.makeTH1D("M_{miss}^{cut}","M_{miss} / MeV","events",bins_mm,"MM_etap");

    Chi2_All = HistFac.makeTH1D("#chi^{2} all","#chi^{2}","",BinSettings(300,0,100),"Chi2_All");
    Chi2_Best = HistFac.makeTH1D("#chi^{2} minimal","#chi^{2}","",BinSettings(300,0,12),"Chi2_Min");

    g_EtaPrime_E = HistFac.makeTH1D("#gamma^{#eta'} E in #eta'","E / MeV","#",BinSettings(500,0,300),"g_EtaPrime_E");
}

EtapOmegaG_MC::ref_perDecayHists_t::ref_perDecayHists_t(SmartHistFactory HistFac)
{
    Steps = HistFac.makeTH1D("Steps", "", "#", BinSettings(10),"steps");

    BinSettings bins_im(1200);

    gg = HistFac.makeTH1D("2#gamma IM","2#gamma IM / MeV","events",bins_im,"gg");

    Proton_Copl = HistFac.makeTH1D("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"Proton_Coplanarity");

    IM_etap = HistFac.makeTH1D("#eta' IM","#eta' IM / MeV","events",BinSettings(300, 700, 1100),"IM_etap");

    BinSettings bins_mm(300, 600, 1300);
    MM_etap = HistFac.makeTH1D("M_{miss}","M_{miss} / MeV","events",bins_mm,"MM_etap");


    Proton_ThetaPhi = HistFac.makeTH2D("p #delta(#theta-#phi)",
                                       "#delta#theta / degree",
                                       "#delta#phi / degree",
                                       BinSettings(100, -10, 10),
                                       BinSettings(100, -30, 30),
                                       "Proton_ThetaPhi"
                                       );
    Proton_Energy = HistFac.makeTH1D("p #delta(E)","#deltaE / MeV","",BinSettings(400,-50,350),"Proton_Energy");
}

void EtapOmegaG_MC::sig_TTree_t::SetBranches()
{
    Tree->Branch("MCTrueIndex", &MCTrueIndex);

    Proton.SetBranches(Tree, "Proton");
    ProtonTrue.SetBranches(Tree, "ProtonTrue");

    Tree->Branch("Chi2", &Chi2);

    g_Pi0_0.SetBranches(Tree, "g_Pi0_0");
    g_Pi0_1.SetBranches(Tree, "g_Pi0_1");
    g_Omega.SetBranches(Tree, "g_Omega");
    g_EtaPrime.SetBranches(Tree, "g_EtaPrime");
    g_EtaPrime_Boosted.SetBranches(Tree, "g_EtaPrime_Boosted");

    Pi0.SetBranches(Tree,"Pi0");
    Omega.SetBranches(Tree, "Omega");
    EtaPrime.SetBranches(Tree, "EtaPrime");
}

void EtapOmegaG_MC::ref_TTree_t::SetBranches()
{
    Tree->Branch("MCTrueIndex", &MCTrueIndex);

    Proton.SetBranches(Tree, "Proton");
}

void EtapOmegaG_MC::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& data = *event.Reconstructed;

    auto& particletree = event.MCTrue->ParticleTree;

    ProcessRef(particletree, data);
    ProcessSig(particletree, data);

    h_TotalEvents->Fill("Total",1);
    if(particletree) {
        // note: this might also match to g p -> eta' eta' p,
        // but this is kinematically forbidden...
        if(utils::ParticleTools::FindParticle(ParticleTypeDatabase::EtaPrime, particletree, 1)) {

            h_TotalEvents->Fill("#eta'", 1);

            if(particletree->IsEqual(treeSignal, utils::ParticleTools::MatchByParticleName)) {
                h_TotalEvents->Fill("Signal",1);
            }
            else if(particletree->IsEqual(treeReference, utils::ParticleTools::MatchByParticleName)) {
                h_TotalEvents->Fill("Reference",1);
            }
        }



    }
}

void EtapOmegaG_MC::ProcessSig(const TParticleTree_t& particletree,
                            const TEventData& data)
{
    TH1D* steps = sig_hists.Steps;

    const auto nParticles = data.Particles.GetAll().size();

    steps->Fill("Seen",1);

    if(data.Trigger.CBEnergySum <= 550)
        return;
    steps->Fill("CBESum>550MeV",1);

    if(nParticles != 5)
        return;
    steps->Fill("nParticles==5",1);

    const auto& photons = data.Particles.Get(ParticleTypeDatabase::Photon);
    const auto& protons = data.Particles.Get(ParticleTypeDatabase::Proton);

    const auto nPhotons = photons.size();
    const auto nProtons = protons.size();

    if(nPhotons != 4)
        return;
    steps->Fill("nPhotons==4",1);

    if(nProtons != 1)
        return;
    steps->Fill("nProtons==1",1);
    const TParticlePtr& proton = protons.front();

    const sig_perDecayHists_t& h = getHistogram(particletree, sig_perDecayHists, sig_TTree.MCTrueIndex);
    const bool other_channel = (unsigned)sig_TTree.MCTrueIndex == sig_perDecayHists.size()-1;
    steps = h.Steps;
    steps->Fill("Seen",1);

    if(!(proton->Candidate->Detector & Detector_t::Type_t::TAPS))
        return;
    steps->Fill("p in TAPS", 1);

    // gamma combinatorics
    assert(photons.size() == 4);
    utils::ParticleTools::FillIMCombinations(h.gg,   2, photons);
    utils::ParticleTools::FillIMCombinations(h.ggg,  3, photons);

    TLorentzVector photon_sum(0,0,0,0);
    for(const auto& p : photons) {
        photon_sum += *p;
    }
    h.gggg->Fill(photon_sum.M());

    // use tagged photon
    for(const auto& th : data.TaggerHits) {
        /// \todo make prompt/random cut
        const TLorentzVector beam_target = th.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TLorentzVector mm = beam_target - photon_sum;
        h.MM_gggg->Fill(mm.M());
    }

    // proton coplanarity
    const double d_phi = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
    h.Proton_Copl->Fill(d_phi);

    const interval<double> Proton_Copl_cut(-19, 19);
    if(!Proton_Copl_cut.Contains(d_phi))
        return;
    steps->Fill("Copl p in 2#sigma",1);

    // fill Goldhaber plot and make cut
    bool is_pi0pi0 = false;
    bool is_pi0eta = false;
    const auto IM_pi0_cut = Pi0.makeCutInterval();
    const auto IM_eta_cut = Eta.makeCutInterval();

    auto goldhaber_comb = std::vector<std::vector<unsigned>>({{0,1,2,3},{0,2,1,3},{0,3,1,2}});
    for(auto i : goldhaber_comb) {
        auto& p = photons;
        const TLorentzVector pair1(*p[i[0]]+*p[i[1]]);
        const TLorentzVector pair2(*p[i[2]]+*p[i[3]]);
        h.IM_gg_gg->Fill(pair1.M(), pair2.M());
        if(   IM_pi0_cut.Stop() > pair1.M()
           && IM_pi0_cut.Stop() > pair2.M())
            is_pi0pi0 = true;
        if(   IM_eta_cut.Contains(pair1.M())
           && IM_pi0_cut.Contains(pair2.M()))
            is_pi0eta = true;
        if(   IM_pi0_cut.Contains(pair1.M())
           && IM_eta_cut.Contains(pair2.M()))
            is_pi0eta = true;
    }
    if(is_pi0pi0)
        return;
    steps->Fill("Not <#pi^{0}#pi^{0}",1);

    if(is_pi0eta)
        return;
    steps->Fill("Not #pi^{0}#eta",1);

    for(auto i : goldhaber_comb) {
        auto& p = photons;
        const TLorentzVector pair1(*p[i[0]]+*p[i[1]]);
        const TLorentzVector pair2(*p[i[2]]+*p[i[3]]);
        h.IM_gg_gg_cut->Fill(pair1.M(), pair2.M());
    }


    // bottom-up assignment of photons using Chi2

    struct result_t {
        double Chi2 = std::numeric_limits<double>::infinity();
        TParticlePtr g_pi0_0;
        TParticlePtr g_pi0_1;
        TParticlePtr g_omega;
        TParticlePtr g_etap;
        TLorentzVector Pi0;
        TLorentzVector Omega;
        TLorentzVector EtaPrime;
    };

    result_t result; // best chi2

    assert(photons.size() == 4);
    vector<unsigned> indices = {0,1,2,3};
    auto comparer = [] (unsigned a, unsigned b) {
        // make 0 and 1 look equal
        if(a==0 && b==1)
            return false;
        return a<b;
    };

    do {
        VLOG(9) << indices;
        // the i vector tells us what particle
        // should be used as daughter particle
        // 0,1 : from pi0
        // 2   : from omega
        // 3   : from eta'

        result_t tmp;

        for(unsigned i=0;i<indices.size();i++) {
            switch(indices[i]) {
            case 0: tmp.g_pi0_0 = photons[i]; break;
            case 1: tmp.g_pi0_1 = photons[i]; break;
            case 2: tmp.g_omega = photons[i]; break;
            case 3: tmp.g_etap  = photons[i]; break;
            }
        }

        tmp.Pi0 = *(tmp.g_pi0_0) + *(tmp.g_pi0_1);

        tmp.Omega = tmp.Pi0 + *tmp.g_omega;
        tmp.EtaPrime = tmp.Omega + *tmp.g_etap;

        const auto Chi2_Pi0 =  std_ext::sqr((tmp.Pi0.M() - Pi0.Mean) / Pi0.Sigma);
        const auto Chi2_Omega = std_ext::sqr((tmp.Omega.M() - Omega.Mean) / Omega.Sigma);
        const auto Chi2_EtaPrime = std_ext::sqr((tmp.EtaPrime.M() - EtaPrime_sig.Mean) / EtaPrime_sig.Sigma);

        tmp.Chi2 = Chi2_Pi0 + Chi2_Omega + Chi2_EtaPrime;

        h.Chi2_All->Fill(tmp.Chi2);

        if(tmp.Chi2<result.Chi2)
            result = tmp;
    }
    while(next_permutation(indices.begin(), indices.end(), comparer));

    // photon assignment successful?
    h.Chi2_Best->Fill(result.Chi2);
    const double Chi2cut = 6;
    if(result.Chi2>=Chi2cut)
        return;
    const string Chi2cut_desc(std_ext::formatter() << "MinChi2<" << Chi2cut);
    steps->Fill(Chi2cut_desc.c_str(),1);

    // boost the eta' photon into rest system of the eta'
    // should be monoenergetic there
    TLorentzVector g_etap_boosted(*result.g_etap);
    g_etap_boosted.Boost(-photon_sum.BoostVector());

    h.g_EtaPrime_E->Fill(g_etap_boosted.E());

    h.IM_etap_omega->Fill(result.EtaPrime.M(), result.Omega.M());
    h.IM_pi0->Fill(result.Pi0.M());

    // was this some unidentified channel?
    if(other_channel) {
        sig_hists.MissedBkg->Fill(utils::ParticleTools::GetDecayString(particletree).c_str(),1);
    }

    // use tagged photon
    for(const auto& th : data.TaggerHits) {
        /// \todo make prompt/random cut
        const TLorentzVector beam_target = th.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TLorentzVector mm = beam_target - result.EtaPrime;
        h.MM_etap->Fill(mm.M());
    }

    // fill tree
    sig_TTree.Proton = *proton;
    if(auto protonTrue = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton, particletree, 1)) {
        sig_TTree.ProtonTrue = *protonTrue;
    }
    else {
        sig_TTree.ProtonTrue.Clear();
    }

    sig_TTree.Chi2 = result.Chi2;
    sig_TTree.g_Pi0_0 = *result.g_pi0_0;
    sig_TTree.g_Pi0_1 = *result.g_pi0_1;
    sig_TTree.g_Omega = *result.g_omega;
    sig_TTree.g_EtaPrime = *result.g_etap;
    sig_TTree.g_EtaPrime_Boosted = utils::ParticleVars(g_etap_boosted, ParticleTypeDatabase::Photon);


    sig_TTree.Pi0 = utils::ParticleVars(result.Pi0, ParticleTypeDatabase::Pi0);
    sig_TTree.Omega = utils::ParticleVars(result.Omega, ParticleTypeDatabase::Omega);
    sig_TTree.EtaPrime = utils::ParticleVars(result.EtaPrime, ParticleTypeDatabase::EtaPrime);

    sig_TTree.Tree->Fill();
}

void EtapOmegaG_MC::ProcessRef(const TParticleTree_t& particletree,
                            const TEventData& data)
{
    TH1D* steps = ref_hists.Steps;

    const auto nParticles = data.Particles.GetAll().size();

    steps->Fill("Seen",1);

    if(data.Trigger.CBEnergySum <= 550)
        return;
    steps->Fill("CBESum>550MeV",1);

    if(nParticles != 3)
        return;
    steps->Fill("nParticles==3",1);

    const auto& photons = data.Particles.Get(ParticleTypeDatabase::Photon);
    const auto& protons = data.Particles.Get(ParticleTypeDatabase::Proton);

    const auto nPhotons = photons.size();
    const auto nProtons = protons.size();

    if(nPhotons != 2)
        return;
    steps->Fill("nPhotons==2",1);

    if(nProtons != 1)
        return;
    steps->Fill("nProtons==1",1);
    const TParticlePtr& proton = protons.front();


    const ref_perDecayHists_t& h = getHistogram(particletree, ref_perDecayHists, ref_TTree.MCTrueIndex);
    const bool other_channel = (unsigned)ref_TTree.MCTrueIndex == ref_perDecayHists.size()-1;
    steps = h.Steps;
    steps->Fill("Seen",1);

    if(!(proton->Candidate->Detector & Detector_t::Type_t::TAPS))
        return;
    steps->Fill("p in TAPS", 1);

    // gamma combinatorics
    assert(photons.size() == 2);

    TLorentzVector photon_sum(0,0,0,0);
    for(const auto& p : photons) {
        photon_sum += *p;
    }
    const double gg_im = photon_sum.M();
    h.gg->Fill(gg_im);

    // proton coplanarity
    const double d_phi = std_ext::radian_to_degree(TVector2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
    h.Proton_Copl->Fill(d_phi);
    const interval<double> Proton_Copl_cut(-19, 19);
    if(!Proton_Copl_cut.Contains(d_phi))
        return;
    steps->Fill("Copl p in 2#sigma",1);

    auto IM_etap_cut = EtaPrime_ref.makeCutInterval();

    h.IM_etap->Fill(photon_sum.M());
    if(!IM_etap_cut.Contains(gg_im))
        return;
    steps->Fill("IM #eta' in 2#sigma",1);


    // was this some unidentified channel?
    if(other_channel) {
        ref_hists.MissedBkg->Fill(utils::ParticleTools::GetDecayString(particletree).c_str(),1);
    }

    // use tagged photon
    for(const auto& th : data.TaggerHits) {
        /// \todo make prompt/random cut
        const TLorentzVector beam_target = th.GetPhotonBeam() + TLorentzVector(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
        const TLorentzVector mm = beam_target - photon_sum;
        h.MM_etap->Fill(mm.M());
        // don't assume we always have a proton...
        if(proton) {
            const double diff_Theta = mm.Theta() - proton->Theta();
            const double diff_Phi = mm.Phi() - proton->Phi();
            h.Proton_ThetaPhi->Fill(std_ext::radian_to_degree(diff_Theta), std_ext::radian_to_degree(diff_Phi));
            h.Proton_Energy->Fill(mm.E() - proton->E());
        }
    }

    // fill tree
    ref_TTree.Proton = *proton;

    ref_TTree.Tree->Fill();
}


void EtapOmegaG_MC::Finish()
{
}

void EtapOmegaG_MC::ShowResult()
{
    canvas c_steps(GetName()+": Overview");
    c_steps << h_TotalEvents
            << sig_hists.Steps << sig_hists.MissedBkg
            << h_TotalEvents // draw it 2x to maintain layout...
            << ref_hists.Steps << ref_hists.MissedBkg
            << endc;

    for(const auto& it : sig_perDecayHists) {
        const sig_perDecayHists_t& h = it.PerDecayHists;
        if(h.IM_etap_omega->GetEntries()==0)
            continue;
        canvas c(GetName()+": "+it.ShortName);
        c << h.Steps
          << h.gg << h.ggg << h.gggg
          << h.MM_gggg
          << h.Proton_Copl
          << drawoption("colz") << h.IM_gg_gg
          << drawoption("colz") << h.IM_gg_gg_cut
          << h.g_EtaPrime_E
          << h.Chi2_All << h.Chi2_Best
          << h.IM_pi0
          << drawoption("colz") << h.IM_etap_omega
          << h.MM_etap
          << endc;
    }

    for(const auto& it : ref_perDecayHists) {
        const ref_perDecayHists_t& h = it.PerDecayHists;
        if(h.IM_etap->GetEntries()==0)
            continue;
        canvas c(GetName()+": "+it.ShortName);
        c << h.Steps
          << h.gg
          << h.Proton_Copl << h.IM_etap
          << h.MM_etap
          << drawoption("colz") << h.Proton_ThetaPhi
          << h.Proton_Energy
          << endc;
    }
}


AUTO_REGISTER_PHYSICS(EtapOmegaG)
AUTO_REGISTER_PHYSICS(EtapOmegaG_MC)
