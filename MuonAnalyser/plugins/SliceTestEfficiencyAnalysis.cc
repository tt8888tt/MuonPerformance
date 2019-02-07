// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/GEMDigi/interface/GEMAMCdataCollection.h"
#include "DataFormats/GEMDigi/interface/GEMGEBdataCollection.h"
#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigiCollection.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "EventFilter/Utilities/interface/json.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class SliceTestEfficiencyAnalysis : public edm::EDAnalyzer {
public:
  explicit SliceTestEfficiencyAnalysis(const edm::ParameterSet&);
  ~SliceTestEfficiencyAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  const GEMEtaPartition* findEtaPartition(const GEMChamber*& chamber, GlobalPoint& tsosGP);
  bool checkEtaPartitionGood(const GEMEtaPartition* part, int amcBx);
  bool getInvariantMass(const reco::Muon mu);
  const GEMRecHit* findMatchedHit(const GEMRecHitCollection* gemRecHits, GEMDetId gemid, LocalPoint locPos);

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  //edm::EDGetTokenT<reco::TrackCollection> muons_;
  edm::EDGetTokenT<GEMAMCdataCollection> amcData_;
  edm::EDGetTokenT<GEMGEBdataCollection> gebStatusCol_;
  edm::EDGetTokenT<GEMVfatStatusDigiCollection> vfatStatusCol_;
  double minPt; bool fiducialCut; bool amcBxCut;

  edm::Handle<GEMAMCdataCollection> amcData;
  edm::Handle<GEMGEBdataCollection> gebStatusCol;
  edm::Handle<GEMVfatStatusDigiCollection> vfatStatusCol;
  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_; 

  static const int MAXCHAMBERS = 36;
  static const int MAXLAYERS = 2;
  static const int MAXROLL = 8;

  int nGEMHitInMuontrack;

  TH1D* h_nGEMHitInMuontrack;
  TH2D* h_inMap;
  TH2D* h_hitMap;

  TH2D* h_hitLumiMap[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_stripLumiMap[MAXCHAMBERS][MAXLAYERS];

  TH1D* h_muon_pt[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_inRoll[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_inStrip[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inVfat[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inStripPerRoll[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inPos[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_inPhi[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_inErrX[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_inErrY[MAXCHAMBERS][MAXLAYERS][MAXROLL];

  TH1D* h_hitErrX[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_hitErrY[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_hitRoll[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_hitStrip[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_hitVfat[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_hitStripPerRoll[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_hitPos[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_hitNstrip[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_resX_etaPart[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_resY_etaPart[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_resPhi_etaPart[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_pullX_etaPart[MAXCHAMBERS][MAXLAYERS][MAXROLL];
  TH1D* h_pullY_etaPart[MAXCHAMBERS][MAXLAYERS][MAXROLL];

  TH1D* h_resX[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_resY[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_resZ[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_resPhi[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_resXvsNstrip[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_pullX[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_pullY[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_muon_pt_matched[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_inStrip_matched[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inPos_matched[MAXCHAMBERS][MAXLAYERS];
  TH1D* h_inPhi_matched[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inPhiVsHitPhi[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inXVsHitX[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inStripVsHitStrip[MAXCHAMBERS][MAXLAYERS];
  TH2D* h_inPhiVsHitStrip[MAXCHAMBERS][MAXLAYERS];

  int nEvents;
  int b_event, b_run, b_lumi;
  int b_nMuons, b_nMuonsInGEMRegion, b_nGEMHits;
  int b_latency;

  bool physicsEvent;

  TTree *t_event;
  edm::Handle<View<reco::Muon> > muons;
  //edm::Handle<reco::TrackCollection> muons;
};

SliceTestEfficiencyAnalysis::SliceTestEfficiencyAnalysis(const edm::ParameterSet& iConfig) :
  nEvents(0)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  //muons_ = consumes<reco::TrackCollection>(iConfig.getParameter<InputTag>("muons"));
  amcData_ = consumes<GEMAMCdataCollection>(iConfig.getParameter<edm::InputTag>("amcData"));
  gebStatusCol_ = consumes<GEMGEBdataCollection>(iConfig.getParameter<edm::InputTag>("gebStatusCol"));
  vfatStatusCol_ = consumes<GEMVfatStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("vfatStatusCol"));

  minPt = iConfig.getParameter<double>("minPt");
  fiducialCut = iConfig.getParameter<bool>("fiducialCut");
  amcBxCut = iConfig.getParameter<bool>("amcBxCut");

  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);

  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("nMuonsInGEMRegion", &b_nMuonsInGEMRegion, "nMuonsInGEMRegion/I");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("event", &b_event, "event/I");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");
  t_event->Branch("latency", &b_latency, "latency/I");

  h_nGEMHitInMuontrack = fs->make<TH1D>("nGEMHitInMuontrack","nGEMHitInMuontrack",4,0,4);
  h_inMap = fs->make<TH2D>("inMap","inMap",18,26.5,31,18,0,9);
  h_hitMap = fs->make<TH2D>("hitMap","hitMap",18,26.5,31,18,0,9);

  Double_t ptbins[] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,200.};
  for (int ichamber=27; ichamber<=30;++ichamber) {
    for (int ilayer=0; ilayer<MAXLAYERS;++ilayer) {
     h_muon_pt[ichamber][ilayer] = fs->make<TH1D>(Form("pt ch %i lay %i",ichamber, ilayer+1),"pt",size(ptbins)-1, ptbins);
      h_hitLumiMap[ichamber][ilayer] = fs->make<TH2D>(Form("hitLumiMap ch %i lay %i",ichamber, ilayer+1),"hitLumiMap",700,0,700,32,1,9);
      h_stripLumiMap[ichamber][ilayer] = fs->make<TH2D>(Form("stripLumiMap ch %i lay %i",ichamber, ilayer+1),"stripLumiMap",700,0,700,400,0,400);
      h_inRoll[ichamber][ilayer] = fs->make<TH1D>(Form("inRoll ch %i lay %i",ichamber, ilayer+1),"inRoll",8,0.5,8.5);
      h_inStrip[ichamber][ilayer] = fs->make<TH1D>(Form("inStrip ch %i lay %i",ichamber, ilayer+1),"inStrip",24,0,384);
      h_inVfat[ichamber][ilayer] = fs->make<TH2D>(Form("inVfat ch %i lay %i",ichamber, ilayer+1),"inVfat",3,0.5,3.5,8,0.5,8.5);
      h_inStripPerRoll[ichamber][ilayer] = fs->make<TH2D>(Form("inStripPerRoll ch %i lay %i",ichamber, ilayer+1),"inStripPerRoll",384,0,384,8,0.5,8.5);
      h_inPos[ichamber][ilayer] = fs->make<TH2D>(Form("inPos ch %i lay %i",ichamber, ilayer+1),"inPos",100,-70,120,100,-260,-110);
      h_inPhi[ichamber][ilayer] = fs->make<TH1D>(Form("inPhi ch %i lay %i",ichamber, ilayer+1),"inPhi",44,-1.8352,-1.1318);

      h_hitRoll[ichamber][ilayer] = fs->make<TH1D>(Form("hitRoll ch %i lay %i",ichamber, ilayer+1),"hitRoll",8,0.5,8.5);
      h_hitStrip[ichamber][ilayer] = fs->make<TH1D>(Form("hitStrip ch %i lay %i",ichamber, ilayer+1),"hitStrip",384,0,384);
      h_hitVfat[ichamber][ilayer] = fs->make<TH2D>(Form("hitVfat ch %i lay %i",ichamber, ilayer+1),"hitVfat",3,0.5,3.5,8,0.5,8.5);
      h_hitStripPerRoll[ichamber][ilayer] = fs->make<TH2D>(Form("hitStripPerRoll ch %i lay %i",ichamber, ilayer+1),"hitStripPerRoll",384,0,384,8,0.5,8.5);
      h_hitPos[ichamber][ilayer] = fs->make<TH2D>(Form("hitPos ch %i lay %i",ichamber, ilayer+1),"hitPos",100,-70,120,100,-260,-110);
      for (int ieta=0; ieta<MAXROLL;++ieta) {
        h_hitErrX[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("hitErrX ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"hitErrX",100,0,1.5);
        h_hitErrY[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("hitErrY ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"hitErrY",100,0,1.5);
        h_inErrX[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("inErrX ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"inErrX",100,0,1.5);
        h_inErrY[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("inErrY ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"inErrY",100,0,2);
        h_hitNstrip[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("hitNstrip ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"hitNstrip",10,0,10);
        h_resX_etaPart[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("resX_etaPart ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"resX",500,-3,3);
        h_resY_etaPart[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("resY_etaPart ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"resY",500,-15,15);
        h_resPhi_etaPart[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("resPhi_etaPart ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"resPhi",500,-0.03,0.03);
        h_pullX_etaPart[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("pullX_etaPart ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"pullX",500,-3,3);
        h_pullY_etaPart[ichamber][ilayer][ieta] = fs->make<TH1D>(Form("pullY_etaPart ch %i lay %i ieta %i",ichamber, ilayer+1, ieta),"pullY",500,-15,15);
      }

      h_resX[ichamber][ilayer] = fs->make<TH1D>(Form("resX ch %i lay %i",ichamber, ilayer+1),"resX",500,-3,3);
      h_resY[ichamber][ilayer] = fs->make<TH1D>(Form("resY ch %i lay %i",ichamber, ilayer+1),"resY",500,-15,15);
      h_resZ[ichamber][ilayer] = fs->make<TH1D>(Form("resZ ch %i lay %i",ichamber, ilayer+1),"resZ",500,-1,1);
      h_resPhi[ichamber][ilayer] = fs->make<TH1D>(Form("resPhi ch %i lay %i",ichamber, ilayer+1),"resPhi",500,-0.03,0.03);
      h_resXvsNstrip[ichamber][ilayer] = fs->make<TH2D>(Form("resXvsNstrip ch %i lay %i",ichamber, ilayer+1),"resXvsNstrip",250,-3,3,15,0,15);
      h_pullX[ichamber][ilayer] = fs->make<TH1D>(Form("pullX ch %i lay %i",ichamber, ilayer+1),"pullX",500,-3,3);
      h_pullY[ichamber][ilayer] = fs->make<TH1D>(Form("pullY ch %i lay %i",ichamber, ilayer+1),"pullY",500,-15,15);
      h_muon_pt_matched[ichamber][ilayer] = fs->make<TH1D>(Form("pt_matched ch %i lay %i",ichamber, ilayer+1),"pt_matched", size(ptbins)-1, ptbins);
      h_inStrip_matched[ichamber][ilayer] = fs->make<TH1D>(Form("inStrip_matched ch %i lay %i",ichamber, ilayer+1),"inStrip",24,0,384);
      h_inPhi_matched[ichamber][ilayer] = fs->make<TH1D>(Form("inPhi_matched ch %i lay %i",ichamber, ilayer+1),"inPhi_matched",44,-1.8352,-1.1318); // -1.8352 ~ -1.1318 (0.1745+0.0026) per 1 ch)
      h_inPos_matched[ichamber][ilayer] = fs->make<TH2D>(Form("inPos_matched ch %i lay %i",ichamber, ilayer+1),"inPos",100,-70,120,100,-260,-110);
      h_inPhiVsHitPhi[ichamber][ilayer] = fs->make<TH2D>(Form("inPhiVsHitPhi ch %i lay %i",ichamber, ilayer+1),"inPos",100,-1.9,-1.1,100,-1.9,-1.1);
      h_inXVsHitX[ichamber][ilayer] = fs->make<TH2D>(Form("inXVsHitX ch %i lay %i",ichamber, ilayer+1),"inPos",100,-70,120,100,-70,120);
      h_inStripVsHitStrip[ichamber][ilayer] = fs->make<TH2D>(Form("inStripVsHitStrip ch %i lay %i",ichamber, ilayer+1),"inPos",400,0,400,400,0,400);
      h_inPhiVsHitStrip[ichamber][ilayer] = fs->make<TH2D>(Form("inPhiVsHitStrip ch %i lay %i",ichamber, ilayer+1),"inPos",100,-1.9,-1.1,400,0,400);
    }
  }
}

SliceTestEfficiencyAnalysis::~SliceTestEfficiencyAnalysis(){
   //cout << "1. Total prop. points:  " << count1 << endl;
   //cout << "2. Points on boundary region:  " << count2 << endl;
   //cout << "3. Points on boundary region (fail matching): " << count3 << endl;
   //cout << "4. Matched with hit on next eta sector:  " << count4 << endl;
}

void
SliceTestEfficiencyAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  physicsEvent = false;
  nEvents++;

  b_nMuons = 0;
  b_nMuonsInGEMRegion = 0;
  b_nGEMHits = 0;
  
  b_event = iEvent.id().event();
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();

  edm::ESHandle<GEMGeometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  const GEMGeometry* GEMGeometry_ = &*hGeom;

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  
  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);
 
  iEvent.getByToken(muons_, muons);

  iEvent.getByToken(amcData_, amcData);
  iEvent.getByToken(gebStatusCol_, gebStatusCol);
  iEvent.getByToken(vfatStatusCol_, vfatStatusCol);

  if (b_run != 319347) return;

  for (auto ch : GEMGeometry_->chambers()) {
    for(auto roll : ch->etaPartitions()) {
      auto rId = roll->id();
      auto recHitsRange = gemRecHits->get(rId);
      b_nGEMHits += recHitsRange.second - recHitsRange.first;
      for (auto hit = recHitsRange.first; hit != recHitsRange.second; ++hit) {
        int vfat = hit->firstClusterStrip()/128 +1;
        h_hitLumiMap[rId.chamber()][rId.layer()-1]->Fill(b_lumi,rId.roll()+vfat/4.);
        h_stripLumiMap[rId.chamber()][rId.layer()-1]->Fill(b_lumi,hit->firstClusterStrip());
      }
      cout << rId.chamber() << "  " << rId.roll() << " : " << roll->pitch() << endl;
    }
  }
  if (b_nGEMHits == 0) return;

  int amcBx = -1;
  for (GEMAMCdataCollection::DigiRangeIterator amcsIt = amcData->begin(); amcsIt != amcData->end(); ++amcsIt){
    auto amcs = (*amcsIt).second;
    for (auto amc = amcs.first; amc != amcs.second; ++amc) {
      amcBx = amc->bx();
    }
  }

  for (auto & mu : *muons) {

    if (!mu.passed(reco::Muon::CutBasedIdTight)) continue;
    if (mu.pt() < minPt) continue;

    const reco::Track* muonTrack = 0;
    if ( mu.globalTrack().isNonnull() ) muonTrack = mu.globalTrack().get();
    else if ( mu.outerTrack().isNonnull()  ) muonTrack = mu.outerTrack().get();
    if (!muonTrack) continue;
    b_nMuons++;

    nGEMHitInMuontrack = 0;
    for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
      if ((*hit)->geographicalId().det() == 2 && (*hit)->geographicalId().subdetId() == 4) {
        nGEMHitInMuontrack++;
      }
    }
    h_nGEMHitInMuontrack->Fill(nGEMHitInMuontrack);

    bool onGEMChamber = false;
    reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
    for (auto chamber : GEMGeometry_->chambers()) {
	  if (chamber->id().chamber() == 1) continue; // ignore chammber 1
	  if (mu.eta() * chamber->id().region() < 0 ) continue;

	  TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),
	  					            chamber->surface());
	  if (!tsos.isValid()) continue;

	  GlobalPoint tsosGP = tsos.globalPosition();
      auto etaPart = findEtaPartition(chamber, tsosGP);
      if (!etaPart) continue; // no eta partition matched

      auto gemid = etaPart->id();
      auto locPos = etaPart->toLocal(tsosGP);
      auto strip = (int) etaPart->strip(locPos);
      auto vfat = ((int) strip/128)+1;
      auto shiftedlocPos = etaPart->centreOfStrip(strip);

      if (fiducialCut) { if (strip < 20 || strip > 364) continue; }
      if (amcBxCut) { if (!checkEtaPartitionGood(etaPart, amcBx)) continue; }

	  h_muon_pt[gemid.chamber()][gemid.layer()-1]->Fill(mu.pt());
	  h_inRoll[gemid.chamber()][gemid.layer()-1]->Fill(gemid.roll());
	  h_inStrip[gemid.chamber()][gemid.layer()-1]->Fill(strip);
	  h_inVfat[gemid.chamber()][gemid.layer()-1]->Fill(vfat, gemid.roll());
	  h_inStripPerRoll[gemid.chamber()][gemid.layer()-1]->Fill(strip, gemid.roll());
	  h_inPos[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.x(), tsosGP.y());
	  h_inMap->Fill(gemid.chamber()+gemid.layer()/2., gemid.roll());
	  h_inPhi[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.phi());
      h_inErrX[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(tsos.localError().positionError().xx());
      h_inErrY[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(tsos.localError().positionError().yy());

      auto matchedHit = findMatchedHit(gemRecHits.product(), gemid, shiftedlocPos);

      // searching near eta sector (boundary)
      /*
      if ( !matchedHit ){
        const TrapezoidalStripTopology* topo_(dynamic_cast<const TrapezoidalStripTopology*>(&(etaPart->topology())));
        const float stripLength(topo_->stripLength());

        const GEMEtaPartition* nextEtaPart = nullptr;
        if (gemid.roll() != 1 && (locPos.y()+2) >  stripLength/2.) nextEtaPart = chamber->etaPartition(gemid.roll()+1);
        if (gemid.roll() != 8 && (locPos.y()-2) < -stripLength/2.) nextEtaPart = chamber->etaPartition(gemid.roll()-1);
        if (nextEtaPart) matchedHit = findMatchedHit(gemRecHits.product(), nextEtaPart->id(), locPos);
      }
      */

      if ( !matchedHit ) continue;

	  auto hitLocPos = matchedHit->localPosition();
      auto hitGlobPos = etaPart->toGlobal(hitLocPos);
      auto hitStrip = matchedHit->firstClusterStrip();
      auto resX = shiftedlocPos.x() - hitLocPos.x();
      auto resY = shiftedlocPos.y() - hitLocPos.y();
      auto resZ = shiftedlocPos.z() - hitLocPos.z();
      auto resPhi = tsosGP.phi() - hitGlobPos.phi();
	  h_inPhiVsHitPhi[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.phi(), hitGlobPos.phi());
	  h_inXVsHitX[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.x(), hitGlobPos.x());
	  h_inStripVsHitStrip[gemid.chamber()][gemid.layer()-1]->Fill(strip, hitStrip);
	  h_inPhiVsHitStrip[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.phi(), hitStrip);

	  if (resX > 5.0) continue;
      onGEMChamber = true;
      LocalError && locErr = tsos.localError().positionError();
      LocalError && hitLocErr = matchedHit->localPositionError();
      auto pullX = resX / std::sqrt(hitLocErr.xx() + locErr.xx());
      auto pullY = resY / std::sqrt(hitLocErr.yy() + locErr.yy());
      h_hitErrX[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(hitLocErr.xx());
      h_hitErrY[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(hitLocErr.yy());

      //Filling histograms
      h_hitRoll[gemid.chamber()][gemid.layer()-1]->Fill(gemid.roll());
	  h_hitStrip[gemid.chamber()][gemid.layer()-1]->Fill(hitStrip);
	  h_hitVfat[gemid.chamber()][gemid.layer()-1]->Fill(((int)hitStrip/128)+1, gemid.roll());
	  h_hitStripPerRoll[gemid.chamber()][gemid.layer()-1]->Fill(hitStrip, gemid.roll());
	  h_hitPos[gemid.chamber()][gemid.layer()-1]->Fill(hitGlobPos.x(), hitGlobPos.y());
	  h_hitNstrip[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(matchedHit->clusterSize());
	  h_hitMap->Fill(gemid.chamber()+gemid.layer()/2., gemid.roll());

	  h_resX[gemid.chamber()][gemid.layer()-1]->Fill(resX);
	  h_resY[gemid.chamber()][gemid.layer()-1]->Fill(resY);
	  h_resZ[gemid.chamber()][gemid.layer()-1]->Fill(resZ);
	  h_resPhi[gemid.chamber()][gemid.layer()-1]->Fill(resPhi);
	  h_pullX[gemid.chamber()][gemid.layer()-1]->Fill(pullX);
	  h_pullY[gemid.chamber()][gemid.layer()-1]->Fill(pullY);
	  h_resX_etaPart[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(resX);
	  h_resY_etaPart[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(resY);
	  h_resPhi_etaPart[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(resPhi);
	  h_pullX_etaPart[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(pullX);
	  h_pullY_etaPart[gemid.chamber()][gemid.layer()-1][gemid.roll()-1]->Fill(pullY);

	  h_muon_pt_matched[gemid.chamber()][gemid.layer()-1]->Fill(mu.pt());
	  h_inStrip_matched[gemid.chamber()][gemid.layer()-1]->Fill(strip);
	  h_inPhi_matched[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.phi());
	  h_inPos_matched[gemid.chamber()][gemid.layer()-1]->Fill(tsosGP.x(), tsosGP.y());
    }

    if (onGEMChamber && getInvariantMass(mu)) physicsEvent = true;
  }
  //if ( physicsEvent ) cout << "physics!: " << b_event << "  " << b_nGEMHits << endl;
  t_event->Fill();
}

const GEMEtaPartition* SliceTestEfficiencyAnalysis::findEtaPartition(const GEMChamber*& chamber, GlobalPoint& tsosGP){
  for (auto etaPart : chamber->etaPartitions()) {
    const LocalPoint locPos = etaPart->toLocal(tsosGP);
    const LocalPoint locPos2D(locPos.x(), locPos.y(), 0);
    const BoundPlane& bps(etaPart->surface());
    if (!bps.bounds().inside(locPos2D)) continue;
    return etaPart;
  }
  return nullptr;
}

bool SliceTestEfficiencyAnalysis::checkEtaPartitionGood(const GEMEtaPartition* part, int amcBx)
{
  GEMDetId rId = part->id();
  auto vfats = vfatStatusCol->get(rId);
  for (auto vfat = vfats.first; vfat != vfats.second; ++vfat) {
    if (vfat->bc() != amcBx ) return false;
    if (int(vfat->quality()) != 0 ) return false;
  }
  return true;
}

const GEMRecHit* SliceTestEfficiencyAnalysis::findMatchedHit(const GEMRecHitCollection* gemRecHits, GEMDetId gemid, LocalPoint locPos)
{
  float dx = 999;
  const GEMRecHit* closestHit = nullptr;
  auto recHitsRange = gemRecHits->get(gemid);
  for (auto hit = recHitsRange.first; hit != recHitsRange.second; ++hit) {
    LocalPoint hitLocPos = hit->localPosition();
    if ( fabs(hitLocPos.x() - locPos.x()) < fabs(dx) ) {
      dx = hitLocPos.x() - locPos.x();
      closestHit = &(*hit);
    }
  }
  return closestHit;
}

//vector<float> SliceTestEfficiencyAnalysis::getInvariantMass(muons, const reco::Muon mu)
bool SliceTestEfficiencyAnalysis::getInvariantMass(const reco::Muon mu)
{
  TLorentzVector mu1tlv;
  mu1tlv.SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), 0.105658);
  
  bool massInRegion = false;
  for (auto & mu2 : *muons) { // double count?
      if (&mu == &mu2) continue;

      if (!mu2.passed(reco::Muon::CutBasedIdTight)) continue;
      if (mu2.pt() < minPt) continue;

      const reco::Track* muonTrack = 0;
      if ( mu2.globalTrack().isNonnull() ) muonTrack = mu2.globalTrack().get();
      else if ( mu2.outerTrack().isNonnull()  ) muonTrack = mu2.outerTrack().get();
      if (!muonTrack) continue;

      TLorentzVector mu2tlv;
      mu2tlv.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), 0.105658);
      auto invMass = (mu1tlv+mu2tlv).M();
      if (abs(invMass-3.096) < 0.05 || // jpsi
          abs(invMass-9.460) < 0.05 || // upsilon
          abs(invMass-91.187) < 0.5) // zboson
            massInRegion = true;
  }
  return massInRegion;
}

void SliceTestEfficiencyAnalysis::beginJob(){}
void SliceTestEfficiencyAnalysis::endJob(){}

void SliceTestEfficiencyAnalysis::beginRun(Run const& run, EventSetup const&){}
void SliceTestEfficiencyAnalysis::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(SliceTestEfficiencyAnalysis);
