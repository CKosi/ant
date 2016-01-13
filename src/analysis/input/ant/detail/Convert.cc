#include "Convert.h"
#include "data/Event.h"

#include "base/Detector_t.h"

#include "tree/TEvent.h"
#include "tree/TCluster.h"
#include "tree/TCandidate.h"
#include "tree/TTaggerHit.h"
#include "tree/TDataRecord.h"
#include "data/Event.h"

#include "base/Detector_t.h"

#include <memory>
#include <limits>

using namespace ant;
using namespace ant::analysis::input;
using namespace ant::analysis::data;
using namespace std;

template <class Cont1, class Cont2>
void Copy(const Cont1& from, Cont2& to) {
    to.reserve(from.size());
    for(auto& from_element : from) {
        to.emplace_back(Converter::Convert(from_element));
    }
}

Event Converter::Convert(const TEvent& event)
{
    Event antevent;

    Copy(event.Candidates,  antevent.Reconstructed.Candidates);
    Copy(event.Tagger.Hits, antevent.Reconstructed.TaggerHits);
    antevent.Reconstructed.AllClusters = event.AllClusters;

    // calcuclate some trigger stuff

    /// @todo The multiplicity is a much harder business, see acqu/root/src/TA2BasePhysics.cc
    /// the code there might only apply to the old trigger system before 2012
    /// so it might be best to implement such algorithms with some nicely designed interface into the
    /// pseudo-detector Trigger in expconfig/detectors

    double Esum = 0.0;
    double TimeEsum = 0.0;
    double TimeE = 0.0;
    for(const TCluster& cluster : antevent.Reconstructed.AllClusters) {
        if(cluster.GetDetectorType() == Detector_t::Type_t::CB) {
            if(isfinite(cluster.Energy)) {
                Esum += cluster.Energy;
                if(isfinite(cluster.Time)) {
                    TimeEsum += cluster.Energy;
                    TimeE += cluster.Energy*cluster.Time;
                }
            }
        }
    }

    auto& triggerinfos = antevent.Reconstructed.Trigger;
    triggerinfos.EventID = event.ID;
    triggerinfos.CBEnergySum = Esum;
    triggerinfos.CBTiming = TimeE/TimeEsum;

    return antevent;
}


CandidatePtr Converter::Convert(const TCandidate& candidate)
{
    return std::make_shared<TCandidate>(candidate);
}

TaggerHit Converter::Convert(const TTaggerHit& taggerhit)
{
    return TaggerHit(
                taggerhit.Electrons.front().Key, /// @bug only first tagger electron is considered
                taggerhit.PhotonEnergy,
                taggerhit.Time
                );
}


//Cluster Converter::Convert(const TCluster& cluster)
//{

//    Cluster cl(
//                cluster.Energy,
//                0.0,
//                cluster.Time,
//                cluster.GetDetectorType(),
//                cluster.CentralElement,
//                cluster.Position
//                );

//    if(cluster.HasFlag(TCluster::Flags_t::Split)) cl.flags.Set(Cluster::Flag::Split);
//    if(cluster.HasFlag(TCluster::Flags_t::TouchesHole)) cl.flags.Set(Cluster::Flag::TouchesHole);

//    double eshort = 0.0;
//    for(const auto& hit : cluster.Hits) {
//       Cluster::Hit anthit;
//       anthit.Channel = hit.Channel;
//       for(const auto& datum : hit.Data) {
//           anthit.Data.emplace_back(static_cast<Channel_t::Type_t>(datum.Type), datum.Value);

//           // sum up short energy
//           ///@todo implement short Energy in TCluster?
//           if(static_cast<Channel_t::Type_t>(datum.Type) == Channel_t::Type_t::IntegralShort) {
//               eshort += datum.Value;
//           }

//       }
//       cl.Hits.emplace_back(move(anthit));
//    }

//    cl.ShortEnergy = eshort;

//    return cl;
//}

