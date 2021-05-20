// class to produce 2 pat::MuonCollections
// one matched to the Park triggers
// another fitered wrt the Park triggers


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TLorentzVector.h>
#include "helper.h"

using namespace std;

constexpr bool debug = false;

class MuonTriggerSelector : public edm::EDProducer {

public:

  explicit MuonTriggerSelector(const edm::ParameterSet &iConfig);
  ~MuonTriggerSelector() override {};

  // Unprescaled triggers
  enum ParkingTriggers_t {
    kHLT_Mu7_IP4,
    kHLT_Mu9_IP5,
    kHLT_Mu9_IP6,
    kHLT_Mu12_IP6
  };
  static const std::vector<ParkingTriggers_t> parkingTriggers_;
  static const std::map<ParkingTriggers_t, std::string> triggerStrings_;

  enum ParkingL1Seeds_t {
    kL1_Mu7er1p5,
    kL1_Mu8er1p5,
    kL1_Mu9er1p5,
    kL1_Mu10er1p5,
    kL1_Mu12er1p5
  };
  static const std::vector<ParkingL1Seeds_t> parkingL1seeds_;
  static const std::map<ParkingL1Seeds_t, std::string> l1Strings_;
  static const std::map<ParkingL1Seeds_t, double> l1SeedPtCuts_; // Determine pass/fail by emulating cuts on L1 objects
  static const std::map<ParkingL1Seeds_t, double> l1SeedEtaCuts_;


private:

  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

  //for trigger match
  const double maxdR_;
  const double l1dR_;

  //for filter wrt trigger
  const double dzTrg_cleaning_; // selects primary vertex

  const double ptMin_;          // min pT in all muons for B candidates
  const double absEtaMax_;      //max eta ""
  const bool softMuonsOnly_;    //cuts muons without soft ID

  edm::EDGetTokenT<l1t::MuonBxCollection> l1Muon_;
};

const std::vector<MuonTriggerSelector::ParkingTriggers_t> MuonTriggerSelector::parkingTriggers_ = {MuonTriggerSelector::kHLT_Mu7_IP4, MuonTriggerSelector::kHLT_Mu9_IP5, MuonTriggerSelector::kHLT_Mu9_IP6, MuonTriggerSelector::kHLT_Mu12_IP6};
const std::map<MuonTriggerSelector::ParkingTriggers_t, std::string> MuonTriggerSelector::triggerStrings_ = {
  {MuonTriggerSelector::kHLT_Mu7_IP4, "HLT_Mu7_IP4"},  
  {MuonTriggerSelector::kHLT_Mu9_IP5, "HLT_Mu9_IP5"},  
  {MuonTriggerSelector::kHLT_Mu9_IP6, "HLT_Mu9_IP6"},  
  {MuonTriggerSelector::kHLT_Mu12_IP6, "HLT_Mu12_IP6"}};

const std::vector<MuonTriggerSelector::ParkingL1Seeds_t> MuonTriggerSelector::parkingL1seeds_ = {
  MuonTriggerSelector::kL1_Mu7er1p5, 
  MuonTriggerSelector::kL1_Mu8er1p5, 
  MuonTriggerSelector::kL1_Mu9er1p5, 
  MuonTriggerSelector::kL1_Mu10er1p5, 
  MuonTriggerSelector::kL1_Mu12er1p5};

const std::map<MuonTriggerSelector::ParkingL1Seeds_t, std::string> MuonTriggerSelector::l1Strings_ = {
  {MuonTriggerSelector::kL1_Mu7er1p5, "L1_Mu7er1p5"},
  {MuonTriggerSelector::kL1_Mu8er1p5, "L1_Mu8er1p5"},
  {MuonTriggerSelector::kL1_Mu9er1p5, "L1_Mu9er1p5"},
  {MuonTriggerSelector::kL1_Mu10er1p5, "L1_Mu10er1p5"},
  {MuonTriggerSelector::kL1_Mu12er1p5, "L1_Mu12er1p5"}
};

const std::map<MuonTriggerSelector::ParkingL1Seeds_t, double> MuonTriggerSelector::l1SeedPtCuts_ = {
  {MuonTriggerSelector::kL1_Mu7er1p5, 7.},
  {MuonTriggerSelector::kL1_Mu8er1p5, 8.},
  {MuonTriggerSelector::kL1_Mu9er1p5, 9.},
  {MuonTriggerSelector::kL1_Mu10er1p5, 10.},
  {MuonTriggerSelector::kL1_Mu12er1p5, 12.}
};

const std::map<MuonTriggerSelector::ParkingL1Seeds_t, double> MuonTriggerSelector::l1SeedEtaCuts_ = {
  {MuonTriggerSelector::kL1_Mu7er1p5, 1.5},
  {MuonTriggerSelector::kL1_Mu8er1p5, 1.5},
  {MuonTriggerSelector::kL1_Mu9er1p5, 1.5},
  {MuonTriggerSelector::kL1_Mu10er1p5, 1.5},
  {MuonTriggerSelector::kL1_Mu12er1p5, 1.5}
};

MuonTriggerSelector::MuonTriggerSelector(const edm::ParameterSet &iConfig):
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
  maxdR_(iConfig.getParameter<double>("maxdR_matching")),
  l1dR_(iConfig.getParameter<double>("l1dR_matching")),
  dzTrg_cleaning_(iConfig.getParameter<double>("dzForCleaning_wrtTrgMuon")),
  ptMin_(iConfig.getParameter<double>("ptMin")),
  absEtaMax_(iConfig.getParameter<double>("absEtaMax")),
  softMuonsOnly_(iConfig.getParameter<bool>("softMuonsOnly")), 
  l1Muon_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1Muon")))
{
  // produce 2 collections: trgMuons (tags) and SelectedMuons (probes & tags if survive preselection cuts)
  produces<pat::MuonCollection>("trgMuons");
  produces<pat::MuonCollection>("SelectedMuons");
  produces<TransientTrackCollection>("SelectedTransientMuons");
}



void MuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexSrc_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();

  if (debug) {
    std::cout << " MuonTriggerSelector::produce " << std::endl;
  }

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  std::vector<pat::TriggerObjectStandAlone> triggeringMuons;
  std::map<unsigned int, std::map<ParkingTriggers_t, bool>> whichTrigger;

  //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  if (debug) {
    std::cout << "\n TRIGGER OBJECTS " << std::endl;
  }

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
    obj.unpackFilterLabels(iEvent, *triggerBits);
    obj.unpackPathNames(names);

    bool isTriggerMuon = false;
    for (unsigned h = 0; h < obj.filterIds().size(); ++h) {
      if (obj.filterIds()[h] == 83) {
        isTriggerMuon = true;
        if (debug) {
          std::cout << "\t   Type IDs:   " << 83 << std::endl;;  //83 = muon
          std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
        }
        break;
      }
    }
    if (!isTriggerMuon) continue;

    bool isParkingTriggerMuon = false;
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {
      std::string filterName = obj.filterLabels()[h];
      if (filterName.find("hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos) {
        isParkingTriggerMuon = true;
        if (debug) {
          std::cout << std::endl << "[debug] Object used in L3 parking filter = " << filterName << std::endl;
          for (unsigned h2 = 0; h2 < obj.filterLabels().size(); ++h2) {
            std::cout << "\t " << obj.filterLabels()[h2] << std::endl;
          }
          for (unsigned h2 = 0; h2 < obj.pathNames().size(); ++h2) {
            std::cout << "\t " << obj.pathNames()[h2] << " / isLF=" << obj.hasPathName(obj.pathNames()[h2], true, false) << " / isL3 " << obj.hasPathName(obj.pathNames()[h2], false, true) <<  std::endl;
          }
          std::cout << "\t   Filters:   " << filterName << std::endl;
        }
        break;
      }
    }

    if (!isParkingTriggerMuon) continue;
    triggeringMuons.push_back(obj);

    // Find which triggers fired for this muon
    unsigned int iTrgMuon(triggeringMuons.size() - 1);

    for (auto& it_trig : parkingTriggers_) {
      whichTrigger[iTrgMuon][it_trig] = false;
      for (auto& it_pathName : obj.pathNames()) {
        if (it_pathName.compare(0, triggerStrings_.at(it_trig).length(), triggerStrings_.at(it_trig)) == 0) {
          whichTrigger[iTrgMuon][it_trig] = true;
        }
      }
    }

    //if(debug){
    //  std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
    // Print trigger object collection and type
    //     std::cout << "\t   Collection: " << obj.collection() << std::endl;
    //}
  }//trigger objects

  if (debug) {
    std::cout << "\n total n of triggering muons = " << triggeringMuons.size() << std::endl;
    for (auto ij : triggeringMuons) {
      std::cout << " >>> components (pt, eta, phi) = " << ij.pt() << " " << ij.eta() << " " << ij.phi() << std::endl;
    }
  }


  std::unique_ptr<pat::MuonCollection>      trgmuons_out   ( new pat::MuonCollection );
  std::unique_ptr<pat::MuonCollection>      muons_out      ( new pat::MuonCollection );
  std::unique_ptr<TransientTrackCollection> trans_muons_out( new TransientTrackCollection );


  //now check for reco muons matched to triggering muons
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonSrc_, muons);

  std::vector<int> muonIsTrigger(muons->size(), 0);
  std::map<int, std::map<ParkingTriggers_t, bool>> recoMuonWhichTrigger;

  for (const pat::Muon & muon : *muons) {
    //this is for triggering muon not really need to be configurable
    unsigned int iMuo(&muon - & (muons->at(0)) );
    if (!(muon.isLooseMuon() && muon.isSoftMuon(PV))) continue;

    float dRMuonMatching = -1.;
    int recoMuonMatching_index = -1;
    int trgMuonMatching_index = -1;
    for (auto& it_trig : parkingTriggers_) {
      recoMuonWhichTrigger[iMuo][it_trig] = false;
    }
    for (unsigned int iTrg = 0; iTrg < triggeringMuons.size(); ++iTrg) {
      float dR = reco::deltaR(triggeringMuons[iTrg], muon);
      if ((dR < dRMuonMatching || dRMuonMatching == -1) && dR < maxdR_) {
        dRMuonMatching = dR;
        recoMuonMatching_index = iMuo;
        trgMuonMatching_index = iTrg;
        if (debug) std::cout << " dR = " << dR
                               << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << " "
                               << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " " << triggeringMuons[iTrg].phi()
                               << std::endl;
        for (auto& it_trig : parkingTriggers_) {
          recoMuonWhichTrigger[iMuo][it_trig] = whichTrigger[iTrg][it_trig];
        }
      }
    }

    //save reco muon
    if (recoMuonMatching_index != -1) {
      pat::Muon recoTriggerMuonCand (muon);
      recoTriggerMuonCand.addUserInt("trgMuonIndex", trgMuonMatching_index);
      trgmuons_out->emplace_back(recoTriggerMuonCand);

      //keep track of original muon index for SelectedMuons collection
      muonIsTrigger[iMuo] = 1;
    }
  }

  // Match reco muons to L1 objects
  std::map<int, std::map<ParkingL1Seeds_t, bool>> recoMuonWhichL1Seed;
  edm::Handle<l1t::MuonBxCollection> l1Muon;
  iEvent.getByToken(l1Muon_, l1Muon);
  std::vector<pair<pat::TriggerObjectStandAlone,int>> l1Objects;

  /*
  for(l1t::MuonBxCollection::const_iterator it=l1Muon->begin(0); it!=l1Muon->end(0); it++){
    pat::TriggerObjectStandAlone l1obj(it->p4());
    l1obj.setCollection("L1Mu");
    l1obj.addTriggerObjectType(trigger::TriggerL1Mu);
    l1obj.setCharge(it->charge());
    l1Objects.emplace_back(l1obj,it->hwIso());
  }
  */
  for (unsigned int muIdx = 0; muIdx < muons->size(); ++muIdx) {
    const pat::Muon& reco_muon = (*muons)[muIdx];

    // Find L1 seed matches in dR cone
    std::vector<l1t::MuonBxCollection::const_iterator> l1_matches;
    for (l1t::MuonBxCollection::const_iterator l1_muon = l1Muon->begin(0); l1_muon != l1Muon->end(0); ++l1_muon) {
      if (l1_muon->hwQual() < 12) {
        continue;
      }
      double reco_eta = reco_muon.eta();
      double reco_phi = reco_muon.phi();
      double l1_eta = l1_muon->etaAtVtx();
      double l1_phi = l1_muon->phiAtVtx();
      if (pow(reco_eta - l1_eta, 2) + pow(reco_phi - l1_phi, 2) < l1dR_ * l1dR_) {
        l1_matches.push_back(l1_muon);
      }
    }

    // Loop over matched L1 seeds, looking for a seed that passes cuts
    for (auto& it_l1_seed : parkingL1seeds_) {
      recoMuonWhichL1Seed[muIdx][it_l1_seed] = false;
      for (auto& it_l1_muon : l1_matches) {
        if ((it_l1_muon->pt() > l1SeedPtCuts_.at(it_l1_seed)) && (fabs(it_l1_muon->eta()) < l1SeedEtaCuts_.at(it_l1_seed))) {
          recoMuonWhichL1Seed[muIdx][it_l1_seed] = true;
          break; 
        }
      }
    }
  }


  // now produce output for analysis (code simplified loop of trg inside)
  // trigger muon + all compatible in dz with any tag
  for (unsigned int muIdx = 0; muIdx < muons->size(); ++muIdx) {
    const pat::Muon& mu = (*muons)[muIdx];
    //selection cuts
    if (mu.pt() < ptMin_) continue;
    if (fabs(mu.eta()) > absEtaMax_) continue;
    //following ID is needed for trigger muons not here
    // anyway it is off in the configuration
    if (softMuonsOnly_ && !mu.isSoftMuon(PV)) continue;

    // same PV as the tag muon, both tag and probe only dz selection
    //bool SkipMuon = true;
    //for (const pat::Muon & trgmu : *trgmuons_out) {
    //  if ( fabs(mu.vz() - trgmu.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0 )
    //    continue;
    //  SkipMuon = false;
    //}
    //// needs decision: what about events without trg muon? now we SKIP them
    //if (SkipMuon)  continue;


    // build transient track
    const reco::TransientTrack muonTT((*(mu.bestTrack())), &(*bFieldHandle)); //sara: check, why not using inner track for muons?
    if (!muonTT.isValid()) continue;

    muons_out->emplace_back(mu);
    muons_out->back().addUserInt("isTriggering", muonIsTrigger[muIdx]);
    for (auto& it_trig : parkingTriggers_) {
      muons_out->back().addUserInt("isTriggering_" + triggerStrings_.at(it_trig), (recoMuonWhichTrigger[muIdx][it_trig] ? 1 : 0));
    }
    for (auto& it_l1_seed : parkingL1seeds_) {
      muons_out->back().addUserInt("isTriggering_" + l1Strings_.at(it_l1_seed), (recoMuonWhichL1Seed[muIdx][it_l1_seed] ? 1 : 0));
    }
    trans_muons_out->emplace_back(muonTT);
  }

  iEvent.put(std::move(trgmuons_out),    "trgMuons");
  iEvent.put(std::move(muons_out),       "SelectedMuons");
  iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
}



DEFINE_FWK_MODULE(MuonTriggerSelector);
