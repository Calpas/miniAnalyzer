// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/EDFilter.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h" 
#include "DataFormats/TauReco/interface/PFTau.h"   
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterAssociation.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TMatrixD.h"
#include <TTree.h>
#include <TH1.h>

#include <math.h>
#include <vector>
#include <map>
#include <string>
#include <array>
#include <algorithm> 
#include <typeinfo>

using namespace std;
using namespace edm;
using namespace reco;


// user's function

// return true if the electron pass the spring15_25ns_cutbased_veto
bool spring15_25ns_cut_veto(const pat::Electron &el, const reco::Vertex &pv);

// set the selected pair lepton
template <class L1, class L2> void selpair( multimap<L1, L2> &mapl1l2);

// return the lepton isolation
double leptoniso(const pat::Tau &tau);
double leptoniso(const pat::Electron &el);
double leptoniso(const pat::Muon &mu);

// return true if a trigger match to the lepton
bool triggermatching(const edm::TriggerNames &triggerNames, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, LorentzVector lp4, vector<string> &vpathLabel);

// return true if the muon passes the medium cut
bool isMediumMuon(const pat::Muon & recoMu); 

// return true if the jet passes the loose cut
bool isLooseJet(const pat::Jet & jet); 

// return the gentau particule
const reco::GenParticle* getGenTau(const pat::Tau& patTau);

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class miniAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>{
   public:
      explicit miniAnalyzer(const edm::ParameterSet&);
      ~miniAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() 				               override;
      virtual void analyze (const edm::Event&, const edm::EventSetup&) override;
      virtual void beginRun(const edm::Run&, const edm::EventSetup&)   override;
      virtual void endRun  (const edm::Run&, const edm::EventSetup&)   override;
      virtual void endJob  () 					       override;

      // fill the preselection tree 
      bool passPostSyncCut();

      // set tree variables 
      int  setgen  (const edm::Event& iEvent,const reco::Candidate &cand);
      void setel   (const pat::Electron &);
      void settau  (const pat::Tau &);
      void settau  (const pat::Tau &, const pat::Tau &);
      void settrigweight  (const pat::Muon &);
      void setidisoweight (const pat::Muon &);
      void setdileptonveto(const vector<pat::Muon> &vmu, const reco::Vertex &pv);
      void setdileptonveto(const vector<pat::Electron> &vel, const reco::Vertex &pv);
      template<class L1, class L2> 
	void setthirdleptonveto(
	    const multimap<L1, L2>      &ml1l2,
	    const vector<pat::Muon>     &vmu,
	    const vector<pat::Electron> &vel,
	    const reco::Vertex &pv
	    ); 
      template<class L1, class L2> 
	void setjet(
	    const vector<pat::Jet> &vjet, 
	    const multimap<L1, L2> &ml1l2
	    );
      template<class L1, class L2>
	void setvar(
	    const multimap<L1, L2> &mapl1l2, 
	    const reco::Vertex     &pv,
	    const pat::MET 	   &pfmet, 
	    const pat::MET 	   &puppiMET, 
	    const TMatrixD         &covMET
	    );

      // set svfit variables
      void svfit(
	  const LorentzVector &l1, 
	  const LorentzVector &l2, 
	  const LorentzVector &met,
	  const TMatrixD      &covMET
	  );


      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<reco::GenJetCollection>        genTaulToken_;
      edm::EDGetTokenT<reco::GenJetCollection>        genTauhToken_;
      edm::EDGetTokenT<double> 		              rhoToken_; 
      edm::EDGetTokenT<GenRunInfoProduct>             runinfoToken_;
      edm::EDGetTokenT<reco::VertexCollection>        vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection>           muonToken_;
      edm::EDGetTokenT<edm::View<reco::GsfElectron>>  electronToken_ ;
      edm::EDGetTokenT<pat::TauCollection>            tauToken_;
      edm::EDGetTokenT<pat::JetCollection>            jetToken_;
      edm::EDGetTokenT<vector<PileupSummaryInfo>>     puToken_; 
      edm::EDGetTokenT<edm::ValueMap<bool>>           eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool>>           eleTightIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<float>>          mvaValuesMapToken_;
      edm::EDGetTokenT<pat::METCollection> 	      metToken_; 
      edm::EDGetTokenT<pat::METCollection> 	      metpuppiToken_; 
      edm::EDGetTokenT<edm::TriggerResults>	      triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>   triggerPrescales_;
      vector<string>                                  mutauFilterName_;
      vector<string>                                  etauFilterName_;
      vector<string>                                  tautauFilterName_;
      vector<string>                                  emuFilterName_;

      // event 
      int run, lumi, evt;  
      // pile Up
      int npv, npu; 
      double rho;
      // leg 1 (leading muon for mutau) 
      double pt_1; 	 
      double phi_1; 	 
      double eta_1; 	 
      double m_1;	 
      double q_1;	 
      double d0_1; 		  // dxy between leading track and first PV
      double dZ_1; 		  // dZ between leading track and first PV
      double mt_1; 		  // use MVAMet
      double pfmt_1; 		  // use as above but using PF MET
      double puppimt_1;           // use as above but using PUPPI MET
      double iso_1; 	
      double id_e_mva_nt_loose_1; // non-triggering electron ID MVA score  
      int    gen_match_1; 	  // type of gen particle matched to object
      double decayModeFindingOldDMs_1; 	 
      double byCombinedIsolationDeltaBetaCorrRaw3Hits_1; 	
      double byLooseCombinedIsolationDeltaBetaCorr3Hits_1;
      double byMediumCombinedIsolationDeltaBetaCorr3Hits_1;
      double byTightCombinedIsolationDeltaBetaCorr3Hits_1;
      double byVLooseIsolationMVArun2v1DBoldDMwLT_1;
      double byLooseIsolationMVArun2v1DBoldDMwLT_1;
      double byMediumIsolationMVArun2v1DBoldDMwLT_1;
      double byTightIsolationMVArun2v1DBoldDMwLT_1;
      double byVTightIsolationMVArun2v1DBoldDMwLT_1;
      double byIsolationMVA3newDMwoLTraw_1; //do not exist	 
      double byIsolationMVA3oldDMwoLTraw_1; //do not exist	 
      double byIsolationMVA3newDMwLTraw_1;  //do not exist 
      double byIsolationMVA3oldDMwLTraw_1;  //do not exist 
      double againstElectronLooseMVA6_1; 	 
      double againstElectronMediumMVA6_1; 	 
      double againstElectronTightMVA6_1; 	 
      double againstElectronVLooseMVA6_1; 	 
      double againstElectronVTightMVA6_1; 	 
      double againstMuonLoose3_1; 	 
      double againstMuonTight3_1; 	 
      double chargedIsoPtSum_1; 	 
      double neutralIsoPtSum_1; 	 
      double puCorrPtSum_1; 	 
      double trigweight_1;  
      double idisoweight_1;
      // leg 2 (trailing tau for mutau) 	 
      double pt_2; 	 
      double phi_2; 	 
      double eta_2; 	 
      double m_2; 	 
      double q_2; 	 
      double d0_2; 	
      double dZ_2; 
      double mt_2; 	
      double pfmt_2; 
      double puppimt_2;
      double iso_2; 	 
      int    gen_match_2;
      double decayModeFindingOldDMs_2; 	 
      double byCombinedIsolationDeltaBetaCorrRaw3Hits_2; 	
      double byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
      double byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
      double byTightCombinedIsolationDeltaBetaCorr3Hits_2;
      double byVLooseIsolationMVArun2v1DBoldDMwLT_2;
      double byLooseIsolationMVArun2v1DBoldDMwLT_2;
      double byMediumIsolationMVArun2v1DBoldDMwLT_2;
      double byTightIsolationMVArun2v1DBoldDMwLT_2;
      double byVTightIsolationMVArun2v1DBoldDMwLT_2;
      double byIsolationMVA3newDMwoLTraw_2; 
      double byIsolationMVA3oldDMwoLTraw_2; 
      double byIsolationMVA3newDMwLTraw_2;  
      double byIsolationMVA3oldDMwLTraw_2;  
      double againstElectronLooseMVA6_2; 	 
      double againstElectronMediumMVA6_2; 	 
      double againstElectronTightMVA6_2; 	 
      double againstElectronVLooseMVA6_2; 	 
      double againstElectronVTightMVA6_2; 	 
      double againstMuonLoose3_2; 	 
      double againstMuonTight3_2; 	 
      double chargedIsoPtSum_2; 	 
      double neutralIsoPtSum_2; 	 
      double puCorrPtSum_2; 	 
      double idisoweight_2; 
      // di-tau system 
      double dr_tt;
      double dphi_tt;
      double dr_t1j1;
      double dr_t2j1;
      double dr_t1j2;
      double dr_t2j2;

      double pt_tt;  	     // use pfmet: a la HIG-13-004 (ptl1+ptl2+MET)
      double mt_tot;     // Use mvamet
      double m_vis; 	 
      double m_sv; 	     // using MarkovChain MC integration svfit.mass()
      double mt_sv; 	     // using MarkovChain MC integration svfit.transverseMass() 
      double pt_sv;  	     // use mvamet 
      double eta_sv; 	     // use mvamet
      double phi_sv; 	     // use mvamet
      double met_sv; 	     // use mvamet: using MarkovChain MC integration svfit.fittedMET().Rho() 
      // met 	 
      double met; 	     // type 1 corrected PF MET
      double metphi; 	     // type 1 corrected PF MET phi
      double puppimet; 	     
      double puppimetphi;    
      double mvamet; 	     
      double mvametphi;      
      double pzetavis; 
      double pzetamiss;      // use mvamet
      double pfpzetamiss;    // use pfmet
      double puppipzetamiss; // use puppimet
      double mvacov00; 	     // mvamet
      double mvacov01;
      double mvacov10; 
      double mvacov11; 	
      double metcov00; 	     // pf met
      double metcov01;
      double metcov10; 
      double metcov11; 	
      // additional jets  
      int    nbtag;
      int    njets; 	     // pt>30 and abs(eta)<4.7
      int    njetspt20;      // pt>20 and abs(eta)<4.7 
      int    njetingap;      // jets passing pfJetID and pt > 30 GeV, in eta gap between the jets 
      int    njetingap20;    // jets passing pfJetID and pt > 20 GeV, in eta gap between the jets
      // leading jet sorted by pt Fill only if corrected jet pt > 20 GeV
      double mjj;
      double jdeta;
      double jdphi;
      double dijetpt;
      double dijetphi;
      double jpt_1; 	 
      double jeta_1; 	 
      double jphi_1; 	 
      double jrawf_1; 	     // factor to be applied to the jet p4 to obtain its uncorrected p4
      double jmva_1; 	     // pu Jet id score
      double jpt_2; 	 
      double jeta_2; 	 
      double jphi_2; 	 
      double jrawf_2; 	     
      double jmva_2; 	
      double hdijetphi;
      double visjeteta; // ??
      double ptvis; // ??

      // extra lepton vetos: evt is vetoed if true 	 
      int    dimuon_veto;  
      int    dielectron_veto;  
      int    extraelec_veto; 
      int    extramuon_veto; 
      double puweight; 	 
      // generator Variables 
      int    NUP;            // number of partons 
      // event Weights 	
      double weight;         // global event weight (MC xSection and num events excluded) 
      // leg 1 (leading muon for mutau) 	 
      double id_m_loose_1; 
      double id_m_medium_1;
      double id_m_tight_1; 
      double id_m_tightnovtx_1;
      double id_m_highpt_1; 
      double id_e_cut_veto_1; 	 
      double id_e_cut_loose_1; 	 
      double id_e_cut_medium_1; 	 
      double id_e_cut_tight_1; 	 
      // leg 2 (trailing tau for mutau) 	 
      double id_m_loose_2; 	 
      double id_m_medium_2; 	 
      double id_m_tight_2; 	 
      double id_m_tightnovtx_2; 	 
      double id_m_highpt_2;	 
      double id_e_cut_veto_2; 	 
      double id_e_cut_loose_2; 	 
      double id_e_cut_medium_2; 	 
      double id_e_cut_tight_2; 	 

      //map to instantiate and clear the trees 
      map<string, double*> mdouble;
      map<string, int*> mint;

      TH1I  *nevt;
      TTree *mutau, *etau, *tautau, *emu;
      TH1I  *nevt_sync;
      TTree *mutau_sync, *etau_sync, *tautau_sync, *emu_sync;

      // svfit channel flag
      string decay;

      // to retrive ele mva value
      multimap<const pat::Electron*, const double > melmva;  

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
miniAnalyzer::miniAnalyzer(const edm::ParameterSet& iConfig):
    prunedGenToken_      (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genpruned"))),
    genTaulToken_        (consumes<reco::GenJetCollection>       (iConfig.getParameter<edm::InputTag>("gentaul"))),
    genTauhToken_        (consumes<reco::GenJetCollection>       (iConfig.getParameter<edm::InputTag>("gentauh"))),
    rhoToken_	         (consumes<double>  		         (iConfig.getParameter<edm::InputTag>("rho"))),
    runinfoToken_        (consumes<GenRunInfoProduct,edm::InRun> (iConfig.getParameter<edm::InputTag>("genruninfo"))),
    vtxToken_            (consumes<reco::VertexCollection>       (iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_	         (consumes<pat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_       (consumes<edm::View<reco::GsfElectron>> (iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_            (consumes<pat::TauCollection>	         (iConfig.getParameter<edm::InputTag>("taus"))),
    jetToken_	         (consumes<pat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jets"))),
    puToken_             (consumes<vector<PileupSummaryInfo>>    (iConfig.getParameter<edm::InputTag>("pu"))),
    eleMediumIdMapToken_ (consumes<edm::ValueMap<bool>>          (iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
    eleTightIdMapToken_  (consumes<edm::ValueMap<bool>>          (iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
    mvaValuesMapToken_   (consumes<edm::ValueMap<float>>         (iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
    metToken_	         (consumes<pat::METCollection>           (iConfig.getParameter<edm::InputTag>("mets"))),
    metpuppiToken_       (consumes<pat::METCollection>           (iConfig.getParameter<edm::InputTag>("metspuppi"))),
    triggerBits_	 (consumes<edm::TriggerResults>		 (iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_      (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_    (consumes<pat::PackedTriggerPrescales>  (iConfig.getParameter<edm::InputTag>("prescales"))),
    mutauFilterName_                                             (iConfig.getParameter<std::vector<std::string> >("mutauFilterName")),
    etauFilterName_                                              (iConfig.getParameter<std::vector<std::string> >("etauFilterName")),
    tautauFilterName_                                            (iConfig.getParameter<std::vector<std::string> >("tautauFilterName")),
    emuFilterName_                                               (iConfig.getParameter<std::vector<std::string> >("emuFilterName")),
    decay{""}

{
  // now do what ever initialization is needed
  usesResource("TFileService");

  // set tree variables adress
  mint["evt"]  = &evt;
  mint["lumi"] = &lumi;
  mint["run"]  = &run;
  mint["npu"]  = &npu;
  mint["npv"]  = &npv;
  mint["gen_match_1"] = &gen_match_1;
  mint["gen_match_2"] = &gen_match_2;
  mint["nbtag"]       = &nbtag;
  mint["njets"]       = &njets;
  mint["njetspt20"]   = &njetspt20;
  mint["njetingap"]   = &njetingap;
  mint["njetingap20"] = &njetingap20 ;
  mint["dimuon_veto"]     = &dimuon_veto;
  mint["dielectron_veto"] = &dielectron_veto;
  mint["extraelec_veto"]  = &extraelec_veto;
  mint["extramuon_veto"]  = &extramuon_veto;
  mint["NUP"]  = &NUP;

  mdouble["rho"]  = &rho;
  // leg 1 (leading muon for mutau) 
  mdouble["pt_1"]   = &pt_1;
  mdouble["phi_1"]  = &phi_1;
  mdouble["eta_1"]  = &eta_1;
  mdouble["m_1"]    = &m_1;
  mdouble["q_1"]    = &q_1;
  mdouble["d0_1"]   = &d0_1;
  mdouble["dZ_1"]   = &dZ_1;
  mdouble["mt_1"]   = &mt_1;
  mdouble["pfmt_1"] = &pfmt_1;
  mdouble["puppimt_1"] = &puppimt_1;
  mdouble["iso_1"]  = &iso_1;
  mdouble["id_e_mva_nt_loose_1"]      = &id_e_mva_nt_loose_1;
  mdouble["decayModeFindingOldDMs_1"] = &decayModeFindingOldDMs_1;
  mdouble["byCombinedIsolationDeltaBetaCorrRaw3Hits_1"]    = &byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  mdouble["byLooseCombinedIsolationDeltaBetaCorr3Hits_1"]  = &byLooseCombinedIsolationDeltaBetaCorr3Hits_1;
  mdouble["byMediumCombinedIsolationDeltaBetaCorr3Hits_1"] = &byMediumCombinedIsolationDeltaBetaCorr3Hits_1;
  mdouble["byTightCombinedIsolationDeltaBetaCorr3Hits_1"]  = &byTightCombinedIsolationDeltaBetaCorr3Hits_1;
  mdouble["byVLooseIsolationMVArun2v1DBoldDMwLT_1"] = &byVLooseIsolationMVArun2v1DBoldDMwLT_1;
  mdouble["byLooseIsolationMVArun2v1DBoldDMwLT_1"]  = &byLooseIsolationMVArun2v1DBoldDMwLT_1;
  mdouble["byMediumIsolationMVArun2v1DBoldDMwLT_1"] = &byMediumIsolationMVArun2v1DBoldDMwLT_1;
  mdouble["byTightIsolationMVArun2v1DBoldDMwLT_1"]  = &byTightIsolationMVArun2v1DBoldDMwLT_1;
  mdouble["byVTightIsolationMVArun2v1DBoldDMwLT_1"] = &byVTightIsolationMVArun2v1DBoldDMwLT_1;
  mdouble["byIsolationMVA3newDMwoLTraw_1"] = &byIsolationMVA3newDMwoLTraw_1;
  mdouble["byIsolationMVA3oldDMwoLTraw_1"] = &byIsolationMVA3oldDMwoLTraw_1;
  mdouble["byIsolationMVA3newDMwLTraw_1"]  = &byIsolationMVA3newDMwLTraw_1;
  mdouble["byIsolationMVA3oldDMwLTraw_1"]  = &byIsolationMVA3oldDMwLTraw_1;
  mdouble["againstElectronLooseMVA6_1"]  = &againstElectronLooseMVA6_1;
  mdouble["againstElectronMediumMVA6_1"] = &againstElectronMediumMVA6_1;
  mdouble["againstElectronTightMVA6_1"]  = &againstElectronTightMVA6_1;
  mdouble["againstElectronVLooseMVA6_1"] = &againstElectronVLooseMVA6_1;
  mdouble["againstElectronVTightMVA6_1"] = &againstElectronVTightMVA6_1;
  mdouble["againstMuonLoose3_1"] = &againstMuonLoose3_1;
  mdouble["againstMuonTight3_1"] = &againstMuonTight3_1;
  mdouble["chargedIsoPtSum_1"]   = &chargedIsoPtSum_1;
  mdouble["neutralIsoPtSum_1"]   = &neutralIsoPtSum_1;
  mdouble["puCorrPtSum_1"] = &puCorrPtSum_1;
  mdouble["trigweight_1"]  = &trigweight_1;
  mdouble["idisoweight_1"] = &idisoweight_1;
  // leg 2 (trailing tau for mutau) 	 
  mdouble["pt_2"]  = &pt_2;
  mdouble["phi_2"] = &phi_2;
  mdouble["eta_2"] = &eta_2;
  mdouble["m_2"]   = &m_2;
  mdouble["q_2"]   = &q_2;
  mdouble["d0_2"]  = &d0_2;
  mdouble["dZ_2"]  = &dZ_2;
  mdouble["mt_2"]  = &mt_2;
  mdouble["pfmt_2"]    = &pfmt_2;
  mdouble["puppimt_2"] = &puppimt_2;
  mdouble["iso_2"] = &iso_2;
  mdouble["decayModeFindingOldDMs_2"] = &decayModeFindingOldDMs_2;
  mdouble["byCombinedIsolationDeltaBetaCorrRaw3Hits_2"]    = &byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  mdouble["byLooseCombinedIsolationDeltaBetaCorr3Hits_2"]  = &byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
  mdouble["byMediumCombinedIsolationDeltaBetaCorr3Hits_2"] = &byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
  mdouble["byTightCombinedIsolationDeltaBetaCorr3Hits_2"]  = &byTightCombinedIsolationDeltaBetaCorr3Hits_2;
  mdouble["byVLooseIsolationMVArun2v1DBoldDMwLT_2"] = &byVLooseIsolationMVArun2v1DBoldDMwLT_2;
  mdouble["byLooseIsolationMVArun2v1DBoldDMwLT_2"]  = &byLooseIsolationMVArun2v1DBoldDMwLT_2;
  mdouble["byMediumIsolationMVArun2v1DBoldDMwLT_2"] = &byMediumIsolationMVArun2v1DBoldDMwLT_2;
  mdouble["byTightIsolationMVArun2v1DBoldDMwLT_2"]  = &byTightIsolationMVArun2v1DBoldDMwLT_2;
  mdouble["byVTightIsolationMVArun2v1DBoldDMwLT_2"] = &byVTightIsolationMVArun2v1DBoldDMwLT_2;
  mdouble["byIsolationMVA3newDMwoLTraw_2"] = &byIsolationMVA3newDMwoLTraw_2;
  mdouble["byIsolationMVA3oldDMwoLTraw_2"] = &byIsolationMVA3oldDMwoLTraw_2;
  mdouble["byIsolationMVA3newDMwLTraw_2"]  = &byIsolationMVA3newDMwLTraw_2;
  mdouble["byIsolationMVA3oldDMwLTraw_2"]  = &byIsolationMVA3oldDMwLTraw_2;
  mdouble["againstElectronLooseMVA6_2"]  = &againstElectronLooseMVA6_2;
  mdouble["againstElectronMediumMVA6_2"] = &againstElectronMediumMVA6_2;
  mdouble["againstElectronTightMVA6_2"]  = &againstElectronTightMVA6_2;
  mdouble["againstElectronVLooseMVA6_2"] = &againstElectronVLooseMVA6_2;
  mdouble["againstElectronVTightMVA6_2"] = &againstElectronVTightMVA6_2;
  mdouble["againstMuonLoose3_2"] = &againstMuonLoose3_2;
  mdouble["againstMuonTight3_2"] = &againstMuonTight3_2;
  mdouble["chargedIsoPtSum_2"]   = &chargedIsoPtSum_2;
  mdouble["neutralIsoPtSum_2"]   = &neutralIsoPtSum_2;
  mdouble["puCorrPtSum_2"] = &puCorrPtSum_2;
  mdouble["idisoweight_2"] = &idisoweight_2;
  // di-tau system 
  mdouble["dr_tt"]   = &dr_tt;
  mdouble["dphi_tt"] = &dphi_tt;
  mdouble["dr_t1j1"] = &dr_t1j1;
  mdouble["dr_t2j1"] = &dr_t2j1;
  mdouble["dr_t1j2"] = &dr_t1j2;
  mdouble["dr_t2j2"] = &dr_t2j2;
  mdouble["pt_tt"]   = &pt_tt;
  mdouble["mt_tot"]  = &mt_tot;
  mdouble["m_vis"]   = &m_vis;
  mdouble["m_sv"]    = &m_sv;
  mdouble["mt_sv"]   = &mt_sv;
  mdouble["pt_sv"]   = &pt_sv;
  mdouble["eta_sv"]  = &eta_sv;
  mdouble["phi_sv"]  = &phi_sv;
  mdouble["met_sv"]  = &met_sv;
  // met 	 
  mdouble["met"]     = &met;
  mdouble["metphi"]  = &metphi;
  mdouble["puppimet"]    = &puppimet;
  mdouble["puppimetphi"] = &puppimetphi;
  mdouble["mvamet"]      = &mvamet;
  mdouble["mvametphi"]   = &mvametphi;
  mdouble["pzetavis"]    = &pzetavis;
  mdouble["pzetamiss"]   = &pzetamiss;
  mdouble["pfpzetamiss"] = &pfpzetamiss;
  mdouble["mvacov00"] = &mvacov00;
  mdouble["mvacov01"] = &mvacov01;
  mdouble["mvacov10"] = &mvacov10;
  mdouble["mvacov11"] = &mvacov11;
  mdouble["metcov00"] = &metcov00;
  mdouble["metcov01"] = &metcov01;
  mdouble["metcov10"] = &metcov10;
  mdouble["metcov11"] = &metcov11;

  mdouble["mjj"]      = &mjj;
  mdouble["jdeta"]    = &jdeta;
  mdouble["jdphi"]    = &jdphi;
  mdouble["dijetpt"]  = &dijetpt;
  mdouble["dijetphi"] = &dijetphi;
  mdouble["jpt_1"]    = &jpt_1;
  mdouble["jeta_1"]   = &jeta_1;
  mdouble["jphi_1"]   = &jphi_1;
  mdouble["jrawf_1"]  = &jrawf_1;
  mdouble["jmva_1"]   = &jmva_1;
  mdouble["jpt_2"]    = &jpt_2;
  mdouble["jeta_2"]   = &jeta_2;
  mdouble["jphi_2"]   = &jphi_2;
  mdouble["jrawf_2"]  = &jrawf_2;
  mdouble["jmva_2"]   = &jmva_2;
  mdouble["puweight"] = &puweight;
  mdouble["weight"]   = &weight;
  mdouble["id_m_loose_1"]      = &id_m_loose_1;
  mdouble["id_m_medium_1"]     = &id_m_medium_1;
  mdouble["id_m_tight_1"]      = &id_m_tight_1;
  mdouble["id_m_tightnovtx_1"] = &id_m_tightnovtx_1;
  mdouble["id_m_highpt_1"]     = &id_m_highpt_1;
  mdouble["id_e_cut_veto_1"]   = &id_e_cut_veto_1;
  mdouble["id_e_cut_loose_1"]  = &id_e_cut_loose_1;
  mdouble["id_e_cut_medium_1"] = &id_e_cut_medium_1;
  mdouble["id_e_cut_tight_1"]  = &id_e_cut_tight_1;
  mdouble["id_m_loose_2"]      = &id_m_loose_2;
  mdouble["id_m_medium_2"]     = &id_m_medium_2;
  mdouble["id_m_tight_2"]      = &id_m_tight_2;
  mdouble["id_m_tightnovtx_2"] = &id_m_tightnovtx_2;
  mdouble["id_m_highpt_2"]     = &id_m_highpt_2;
  mdouble["id_e_cut_veto_2"]   = &id_e_cut_veto_2;
  mdouble["id_e_cut_loose_2"]  = &id_e_cut_loose_2;
  mdouble["id_e_cut_medium_2"] = &id_e_cut_medium_2;
  mdouble["id_e_cut_tight_2"]  = &id_e_cut_tight_2;

  //initialise tree variables
  for(auto it=mdouble.begin(); it!=mdouble.end(); ++it) *(it->second)=-10000;
  for(auto it=mint.begin();    it!=mint.end();    ++it) *(it->second)=-10000;
}


miniAnalyzer::~miniAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------


void
miniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // initial evt's number
  nevt->Fill(0);
  //nevt_sync->Fill(0);

  // event ID 
  run = iEvent.id().run(); lumi = iEvent.id().luminosityBlock(); evt = iEvent.id().event();

  // vertices
  edm::Handle<reco::VertexCollection> vertices;   iEvent.getByToken(vtxToken_, vertices);
  npv = vertices->size();
  const reco::Vertex &pv = vertices->front();

  // pileup
  if (!iEvent.isRealData()) {
    Handle<std::vector< PileupSummaryInfo>> puInfo; iEvent.getByToken(puToken_, puInfo);
    npu=puInfo->begin()->getTrueNumInteractions();
  }
  edm::Handle<double> rhoHandle; iEvent.getByToken(rhoToken_, rhoHandle); rho = *rhoHandle;

  //Handle < LHEEventProduct > lheEPHandle;
  //lheEPHandle.getByLabel (iEvent, "externalLHEProducer");
  //NUP=lheEPHandle->hepeup().NUP; //not sure this is correct!!
  //cout<<"NUP: "<<NUP<<endl;

  //pfmet
  edm::Handle<pat::METCollection> mets; iEvent.getByToken(metToken_, mets);
  const pat::MET &pfmet = mets->front();
  met=pfmet.pt(); 	  
  metphi=pfmet.phi();
  TMatrixD covMET(2, 2); 
  covMET[0][0] = metcov00 = pfmet.getSignificanceMatrix()(0, 0); 
  covMET[0][1] = metcov01 = pfmet.getSignificanceMatrix()(0, 1); 
  covMET[1][0] = metcov10 = pfmet.getSignificanceMatrix()(1, 0);
  covMET[1][1] = metcov11 = pfmet.getSignificanceMatrix()(1, 1); 	

  //puppimet
  edm::Handle<pat::METCollection> metspuppi; iEvent.getByToken(metpuppiToken_, metspuppi);
  const pat::MET &puppiMET = metspuppi->front();
  puppimet=puppiMET.pt(); 	  
  puppimetphi=puppiMET.phi();

  //mvamet??
  //mvacov00 = mvamet.getSignificanceMatrix()(0, 0); 
  //mvacov01 = mvamet.getSignificanceMatrix()(0, 1); 
  //mvacov10 = mvamet.getSignificanceMatrix()(1, 0);
  //mvacov11 = mvamet.getSignificanceMatrix()(1, 1); 	

  // handle object
  edm::Handle<pat::MuonCollection>  muons; 		 iEvent.getByToken(muonToken_, muons); 
  edm::Handle<edm::View<reco::GsfElectron> > electrons;  iEvent.getByToken(electronToken_, electrons);
  edm::Handle<edm::ValueMap<bool>>  medium_id_decisions; iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  edm::Handle<edm::ValueMap<bool>>  tight_id_decisions;  iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
  edm::Handle<edm::ValueMap<float>> mvaValues;           iEvent.getByToken(mvaValuesMapToken_, mvaValues);
  edm::Handle<pat::TauCollection>   taus;   		 iEvent.getByToken(tauToken_, taus);   
  edm::Handle<pat::JetCollection>   jets;   		 iEvent.getByToken(jetToken_, jets);   


  ////////////////////////
  //		        // 
  // Objet preselection //
  //		        //	
  ////////////////////////
  // pre-selected object
  vector<pat::Muon>     vmu;
  vector<pat::Electron> vel;
  vector<pat::Tau>      vtau;
  vector<pat::Jet>      vjet;


  //////////
  //      //
  // Muon //
  //      //
  //////////
  //cout<<"\nentering patmu, size: "<<muons->size()<<endl;
  for (const pat::Muon &mu : *muons) {
    if(mu.pt()<=10)continue; 
    if(fabs(mu.eta())>=2.4) continue;
    if(fabs(mu.bestTrack()->dxy(pv.position()))>=0.045) continue;
    if(fabs(mu.bestTrack()->dz (pv.position()))>=0.2) continue;
    if(!isMediumMuon(mu)) continue;
    vmu.push_back(mu);
  }
  //for(auto i:vmu) cout<<"selmu pt/eta/phi/iso: "<<i.pt()<<" : "<<i.eta()<<" : "<<i.phi()<<" : "<< leptoniso(i) <<endl;


  //////////////
  //          //
  // Electron //
  //          //
  //////////////
  //cout<<"\nentering patele, size: "<<electrons->size()<<endl;
  // (ECAL barrel: |eta| < 1.479; endcaps: 1.48 < |eta| < 3.0)
  melmva.clear(); // retrieve the electron mva via it reference
  vel.reserve(20); // to keep the reference strore in the vector and avoid reallocation!!

  //cout<<"\nentering patel, size: "<<electrons->size()<<endl;
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);
    if(el->pt()<=13)continue; 
    if(fabs(el->eta())>=2.5) continue;
    if(fabs(el->bestTrack()->dxy(pv.position()))>=0.045) continue;
    if(fabs(el->bestTrack()->dz (pv.position()))>=0.2) continue;
    double elmva {(*mvaValues)[el]};
    if      ( fabs(el->superCluster()->eta())<0.8   && 			                 	elmva < 0.967083 ) continue; //barrel 1
    else if ( fabs(el->superCluster()->eta())>0.8   && fabs(el->superCluster()->eta())<1.479 && elmva < 0.929117 ) continue; //barrel 2
    else if ( fabs(el->superCluster()->eta())>1.479 && fabs(el->superCluster()->eta())<3     && elmva < 0.726311 ) continue; //endcap
    if(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)>1)  continue;
    const edm::Ptr<pat::Electron> elPatPtr(el);
    if(!elPatPtr->passConversionVeto()) continue;
    vel.push_back(*elPatPtr);
    melmva.insert(std::pair<const pat::Electron*, const double >( &vel.back(), elmva)); // keep the same el ref of the vector!!
  }
  //for(auto &i : vel) cout<<"selel pt/eta/phi/iso: "<<i.pt()<<" : "<<i.eta()<<" : "<<i.phi()<<" : "<< leptoniso(i)<<endl;


  /////////
  //     //
  // Tau //
  //     //
  /////////
  //cout<<"\nentering pattau, size: "<<taus->size()<<endl;
  for (const pat::Tau &tau : *taus) {
    if(tau.pt()<=20) continue;
    if(fabs(tau.eta())>=2.3) continue;
    if(tau.tauID("decayModeFinding")<=0.5) continue;
    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
    if(fabs(packedLeadTauCand->dz())>=0.2) continue; //The PackedCandidate::dz() method is wrt. the first PV by default 
    if(abs(tau.charge())!=1) continue;
    vtau.push_back(tau);
  }
  //for(auto i:vtau) cout<<"seltau pt/eta/phi/iso: "<<i.pt()<<" : "<<i.eta()<<" : "<<i.phi()<<" : "<< -leptoniso(i) <<endl;


  /////////
  //     //
  // Jet //
  //     //
  /////////
  //cout<<"\nentering patjet, size: "<<jets->size()<<endl;
  for (const pat::Jet &jet : *jets) {
    if(jet.pt()<=20)         continue;
    if(fabs(jet.eta())>=4.7) continue;
    if(!isLooseJet(jet))     continue;
    vjet.push_back(jet);
  }

  //cout<<"vmu/vel/vtau/vjet size: "<<vmu.size()<<"/"<<vel.size()<<"/"<<vtau.size()<<"/"<<vjet.size()<<endl;


  /////////////
  //         //
  // Trigger //
  //         //
  ///////////// 
  edm::Handle<edm::TriggerResults> triggerBits;		              iEvent.getByToken(triggerBits_, triggerBits);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects; iEvent.getByToken(triggerObjects_, triggerObjects);
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales; 	      iEvent.getByToken(triggerPrescales_, triggerPrescales);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  // lepton trigger flag
  bool mutrig, etrig, tautrig, tau1trig, tau2trig;


  /////////
  //     //
  // EMu //
  //     //
  /////////
  multimap<const pat::Muon*,     const pat::Electron*> mmuel;
  multimap<const pat::Electron*, const pat::Muon*>     melmu;
  if(vel.size()>=1 && vmu.size()>=1){
    for( const auto &el:vel ){
      etrig = triggermatching(triggerNames, triggerObjects, el.p4(), emuFilterName_);
      for( const auto &mu:vmu ){
	if(deltaR( el.p4(), mu.p4() )<=0.3 ) continue;
	mutrig = triggermatching(triggerNames, triggerObjects, mu.p4(), emuFilterName_);
	if(!etrig && !mutrig) continue;
	decay="em";  
	// fill tree only if a lepton has pt >24
	if(max(el.pt(), mu.pt())<=24) continue;
	mmuel.insert(std::pair<const pat::Muon *, const pat::Electron *>(&mu, &el)); 
      }// vmu
    }// vel 
    selpair(mmuel);
    if(mmuel.size()!=0){
      cout<<"\nem.....";
      // swap such that el=leg1 and mu=leg2 t make it sounds EMU!! 
      melmu.insert(std::pair<const pat::Electron *, const pat::Muon *>(mmuel.begin()->second, mmuel.begin()->first));
      setvar(melmu, pv, pfmet, puppiMET, covMET); // commom lepton variable. Tree cleaning is done here!!
      setel(*melmu.begin()->first); // el variable
      setjet(vjet, melmu); // jet variable
      setthirdleptonveto(melmu, vmu, vel, pv);
      gen_match_1 = setgen(iEvent, *melmu.begin()->first);
      gen_match_2 = setgen(iEvent, *melmu.begin()->second);
      //emu_sync ->Fill();
      //cout<<"\nselemu: "<<evt<<" : "<<pt_1<<" : "<<eta_1<<" : "<<phi_1<<" : "<<pt_2<<" : "<<eta_2<<" : "<<phi_2<<endl<<endl;
      if(passPostSyncCut()) emu ->Fill();
    }
  }// EMu


  ///////////
  //	   // 
  // MuTau //
  //	   //	
  ///////////
  multimap<const pat::Muon*, const pat::Tau*> mmutau;
  if(vmu.size()>=1 && vtau.size()>=1){
    for( const auto &mu:vmu ){
      if(mu.pt()<=23)continue; // new sync
      if(fabs(mu.eta())>=2.4) continue; // new sync
      mutrig = triggermatching(triggerNames, triggerObjects, mu.p4(), mutauFilterName_); // for now always return true
      for( const auto &tau:vtau ){
	if(deltaR( mu.p4(), tau.p4() )<=0.5 ) continue;
	tautrig = triggermatching(triggerNames, triggerObjects, tau.p4(), mutauFilterName_);
	if(!mutrig && !tautrig) continue; // for now always be true
	decay="mt";  
	// store preselected pairs 
	mmutau.insert(std::pair<const pat::Muon *, const pat::Tau *>(&mu, &tau)); 
      }// vtau
    }// vmu 
    // get the selected pair vs pt/iso 
    selpair(mmutau);

    if(mmutau.size()!=0){
      cout<<"\nmt.....";
      setvar(mmutau, pv, pfmet, puppiMET, covMET); 
      settau(*mmutau.begin()->second); 
      setjet(vjet, mmutau); 
      setdileptonveto(vmu, pv);
      setthirdleptonveto(mmutau, vmu, vel, pv);	
      gen_match_1 = setgen(iEvent, *mmutau.begin()->first);
      gen_match_2 = setgen(iEvent, *mmutau.begin()->second);
      //mutau_sync ->Fill();
      //cout<<"\nselmutau: "<<evt<<" : "<<pt_1<<" : "<<eta_1<<" : "<<phi_1<<" : "<<pt_2<<" : "<<eta_2<<" : "<<phi_2<<endl<<endl;
      if(passPostSyncCut()) mutau ->Fill();
    }
  }// MuTau


  //////////
  //      //
  // ETau //
  //      //
  //////////
  multimap<const pat::Electron*, const pat::Tau*> meltau;
  if(vel.size()>=1 && vtau.size()>=1){
    for( const auto &el:vel ){
      if(el.pt()<=26)continue; // new sync !!
      if(fabs(el.eta())>=2.1) continue;
      etrig = triggermatching(triggerNames, triggerObjects, el.p4(), etauFilterName_);
      for( const auto &tau:vtau ){
	if(deltaR( el.p4(), tau.p4() )<=0.5 ) continue;
	tautrig = triggermatching(triggerNames, triggerObjects, tau.p4(), etauFilterName_);
	if(!etrig && !tautrig) continue;
	decay="et";  
	meltau.insert(std::pair<const pat::Electron *, const pat::Tau *>(&el, &tau)); 
      }
    }
    selpair(meltau);
    if(meltau.size()!=0){
      cout<<"\net.....";
      setvar(meltau, pv, pfmet, puppiMET, covMET); // commom lepton variable. Tree cleaning is done here!!
      setel(*meltau.begin()->first); // el variable
      settau(*meltau.begin()->second); // tau variable
      setjet(vjet, meltau); // jet variable
      setdileptonveto(vel, pv);
      setthirdleptonveto(meltau, vmu, vel, pv);	
      gen_match_1 = setgen(iEvent, *meltau.begin()->first);
      gen_match_2 = setgen(iEvent, *meltau.begin()->second);
      //etau_sync ->Fill();
      //cout<<"\nseletau: "<<evt<<" : "<<pt_1<<" : "<<eta_1<<" : "<<phi_1<<" : "<<pt_2<<" : "<<eta_2<<" : "<<phi_2<<endl<<endl;
      if(passPostSyncCut()) etau ->Fill();
    }
  }// ETau


  ////////////
  //        //
  // TauTau //
  //        //
  ////////////
  multimap<const pat::Tau*, const pat::Tau*> mtautau;
  multimap<const pat::Tau*, const pat::Tau*> tmpmtautau;
  if( vtau.size()>=2){
    for( const auto &tau1:vtau ){
      if(tau1.pt()<=40) continue;
      if(fabs(tau1.eta())>=2.1) continue;
      //tau1trig = triggermatching(triggerNames, triggerObjects, tau1.p4(), tauFilterName_);
      for( const auto &tau2:vtau ){
	if(tau2.pt()<=40) continue;
	if(fabs(tau2.eta())>=2.1) continue;
	if(deltaR( tau1.p4(), tau2.p4() )<=0.5 ) continue;
	//tau2trig = triggermatching(triggerNames, triggerObjects, tau1.p4(), tauFilterName_);
	//if(!tau1trig && !tau2trig) continue;
	decay  = "tt";  
	mtautau.insert(std::pair<const pat::Tau *, const pat::Tau *>(&tau1, &tau2)); 
      }
    }
    selpair(mtautau);
    if(mtautau.size()!=0){
      cout<<"\nTauTau.....\n";

      // leading tau has highest pt
      if(mtautau.begin()->first->pt()< mtautau.begin()->second->pt()){ // if leg2 has highest pt
	tmpmtautau.insert(std::pair<const pat::Tau *, const pat::Tau *>(mtautau.begin()->second, mtautau.begin()->first)); // tmp...
	mtautau.swap(tmpmtautau); // sawp leg1 and leg2
      }

      // set and fill tree variables
      setvar(mtautau, pv, pfmet, puppiMET, covMET); 
      settau(*mtautau.begin()->first, *mtautau.begin()->second);
      setjet(vjet, mtautau); // jet variable
      setthirdleptonveto(mtautau, vmu, vel, pv);	
      gen_match_1 = setgen(iEvent, *mtautau.begin()->first);
      gen_match_2 = setgen(iEvent, *mtautau.begin()->second);
      //tautau_sync ->Fill();
      //cout<<"\nseltautau: "<<evt<<" : "<<pt_1<<" : "<<eta_1<<" : "<<phi_1<<" : "<<pt_2<<" : "<<eta_2<<" : "<<phi_2<<endl<<endl;
      if(passPostSyncCut()) tautau ->Fill();
    }
  }// TauTau

}//analysis


// ------------ user codes -----------


//////////////////////////////////////////
//			  	        //	
// set the gen match flag for el and mu //
//			  	        //	
//////////////////////////////////////////
int  miniAnalyzer::setgen(const edm::Event& iEvent, const reco::Candidate &cand){

  //cout<<"setgen()\n";
  Handle<edm::View<reco::GenParticle> > pruned; iEvent.getByToken(prunedGenToken_, pruned);  
  Handle<reco::GenJetCollection> gentaul;       iEvent.getByToken(genTaulToken_,   gentaul);
  Handle<reco::GenJetCollection> gentauh; 	iEvent.getByToken(genTauhToken_,   gentauh);

  //cout<<"prun size: "<< pruned ->size() <<endl;
  //cout<<"taul size: "<< gentaul->size()<<endl;
  //cout<<"tauh size: "<< gentauh->size()<<endl;

  int gid{6}; // return 6 if failed all other
  double tmpdr{1000}; // test deltaR
  double dr;

  // loop over gen part and select the closest
  for(const auto &gen : *pruned){
    dr = deltaR(gen.p4(), cand.p4());
    if(dr       >= 0.2) continue;
    if(gen.pt() <=   8) continue;
    if(!gen.statusFlags().isPrompt()) continue;

    if     (dr<tmpdr && abs(gen.pdgId())==11) { tmpdr=dr; gid=1; }
    else if(dr<tmpdr && abs(gen.pdgId())==13) { tmpdr=dr; gid=2; }
  }
  //cout<<"gid_prun: "<<gid<<endl;

  // loop over lepton tau
  for(const auto &gen : *gentaul){
    dr = deltaR(gen.p4(), cand.p4());
    if(dr >= 0.2) continue;
    if(dr >  tmpdr) continue;
    if(gen.pt() <= 8) continue;
    tmpdr=dr; 

    // it is a jet, so we need to look at it's constituants to see if el or mu
    std::vector< const GenParticle * > vgp{gen.getGenConstituents()};
    //cout<<"gentaul const size(): "<<vgp.size()<<endl;
    for(auto &gencons : vgp) {
      if     (abs(gencons->pdgId())==11) gid=3;
      else if(abs(gencons->pdgId())==13) gid=4; 
    }
  }
  //cout<<"gid_l: "<<gid<<endl;

  // loop over had tau
  for(const auto &gen : *gentauh){
    dr = deltaR(gen.p4(), cand.p4());
    if(dr >= 0.2) continue;
    if(dr >  tmpdr) continue;
    if(gen.pt() <= 15) continue;
    gid=5;
  }
  //cout<<"gid_h: "<<gid<<endl;

  return gid;
}


/////////////////////////////////
//			       //
// fill tree for jet variables //
//			       //
/////////////////////////////////
template <class L1, class L2> 
void miniAnalyzer::setjet(const vector<pat::Jet> &vjet, const multimap<L1, L2> &ml1l2){
  if(vjet.size()==0) return;
  decltype(*ml1l2.begin()->first)  l1{*ml1l2.begin()->first};
  decltype(*ml1l2.begin()->second) l2{*ml1l2.begin()->second};
  // additional jets vs the selected pair	 
  int countjet{0}; nbtag=0; njets=0; njetspt20=0;
  vector<LorentzVector> vj1, vj2;

  //cout<<"setjet(): vjetsize "<<vjet.size()<<endl;
  for( const auto &jet:vjet ){
    //cout<<"entering patjet..."<<endl;
    if(deltaR( l1.p4(), jet.p4() )<=0.5 ) continue;
    if(deltaR( l2.p4(), jet.p4() )<=0.5 ) continue;
    if(jet.pt()<=20) continue;
    if(fabs(jet.eta())>=4.7) continue;
    njetspt20++; 
    if(jet.pt()>30) njets++; 	 
    if(fabs(jet.eta())<2.4 && jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.8) nbtag++; 	 
    if(countjet==0) vj1.push_back(jet.p4()); 
    if(countjet==1) vj2.push_back(jet.p4()); 
    countjet++;
  }// vjet

  if(vj1.size()!=0){
    jpt_1  = vj1[0].pt();
    jeta_1 = vj1[0].eta(); 	 
    jphi_1 = vj1[0].phi();

    dr_t1j1  = deltaR(l1.p4(), vj1[0]);
    dr_t2j1  = deltaR(l2.p4(), vj1[0]);
  }
  if(vj2.size()!=0){
    jpt_2  = vj2[0].pt(); 	 
    jeta_2 = vj2[0].eta(); 	 
    jphi_2 = vj2[0].phi(); 

    dr_t1j2  = deltaR(l1.p4(), vj2[0]);
    dr_t2j2  = deltaR(l2.p4(), vj2[0]);

    mjj   = (vj1[0]+vj2[0]).M();
    jdeta = vj1[0].Eta()-vj2[0].Eta();
    jdphi = vj1[0].Phi()-vj2[0].Phi();

    dijetpt  = jpt_1+jpt_2;
    dijetphi = jphi_1+jphi_2;
    
  }
} // setjet


/////////////////////////////////
//			       //
// fill tree for ele variables //
//			       //
/////////////////////////////////
void miniAnalyzer::setel(const pat::Electron &el){ id_e_mva_nt_loose_1 = melmva.find(&el)->second; };


//////////////////////////////////////
//			            //
// fill tree for tau, tau variables //
//			            //
//////////////////////////////////////
void miniAnalyzer::settau(const pat::Tau &tau1, const pat::Tau &tau2){
  //cout<<"settau(tau, tau): \n";
  decayModeFindingOldDMs_1 = tau1.tauID("decayModeFinding"); 	 
  byCombinedIsolationDeltaBetaCorrRaw3Hits_1    = tau1.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"); 	
  byLooseCombinedIsolationDeltaBetaCorr3Hits_1  = tau1.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
  byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = tau1.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
  byTightCombinedIsolationDeltaBetaCorr3Hits_1  = tau1.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
  byVLooseIsolationMVArun2v1DBoldDMwLT_1 = tau1.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
  byLooseIsolationMVArun2v1DBoldDMwLT_1  = tau1.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
  byMediumIsolationMVArun2v1DBoldDMwLT_1 = tau1.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
  byTightIsolationMVArun2v1DBoldDMwLT_1  = tau1.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
  byVTightIsolationMVArun2v1DBoldDMwLT_1 = tau1.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
  againstElectronLooseMVA6_1  = tau1.tauID("againstElectronLooseMVA6"); 	 
  againstElectronMediumMVA6_1 = tau1.tauID("againstElectronMediumMVA6"); 	 
  againstElectronTightMVA6_1  = tau1.tauID("againstElectronTightMVA6"); 	 
  againstElectronVLooseMVA6_1 = tau1.tauID("againstElectronVLooseMVA6"); 	 
  againstElectronVTightMVA6_1 = tau1.tauID("againstElectronVTightMVA6"); 	 
  againstMuonLoose3_1 = tau1.tauID("againstMuonLoose3"); 	 
  againstMuonTight3_1 = tau1.tauID("againstMuonTight3"); 	 
  chargedIsoPtSum_1   = tau1.tauID("chargedIsoPtSum");
  neutralIsoPtSum_1   = tau1.tauID("neutralIsoPtSum");
  puCorrPtSum_1       = tau1.tauID("puCorrPtSum");
  iso_1 = chargedIsoPtSum_1 + max(0., neutralIsoPtSum_1 - 0.2* puCorrPtSum_1);	
  pat::PackedCandidate const* packedLeadTauCand1 = dynamic_cast<pat::PackedCandidate const*>(tau1.leadChargedHadrCand().get());
  d0_1   = packedLeadTauCand1->dxy();
  dZ_1   = packedLeadTauCand1->dz();

  decayModeFindingOldDMs_2 = tau2.tauID("decayModeFinding"); 	 
  byCombinedIsolationDeltaBetaCorrRaw3Hits_2    = tau2.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"); 	
  byLooseCombinedIsolationDeltaBetaCorr3Hits_2  = tau2.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
  byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = tau2.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
  byTightCombinedIsolationDeltaBetaCorr3Hits_2  = tau2.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
  byVLooseIsolationMVArun2v1DBoldDMwLT_2 = tau2.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
  byLooseIsolationMVArun2v1DBoldDMwLT_2  = tau2.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
  byMediumIsolationMVArun2v1DBoldDMwLT_2 = tau2.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
  byTightIsolationMVArun2v1DBoldDMwLT_2  = tau2.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
  byVTightIsolationMVArun2v1DBoldDMwLT_2 = tau2.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
  againstElectronLooseMVA6_2  = tau2.tauID("againstElectronLooseMVA6"); 	 
  againstElectronMediumMVA6_2 = tau2.tauID("againstElectronMediumMVA6"); 	 
  againstElectronTightMVA6_2  = tau2.tauID("againstElectronTightMVA6"); 	 
  againstElectronVLooseMVA6_2 = tau2.tauID("againstElectronVLooseMVA6"); 	 
  againstElectronVTightMVA6_2 = tau2.tauID("againstElectronVTightMVA6"); 	 
  againstMuonLoose3_2 = tau2.tauID("againstMuonLoose3"); 	 
  againstMuonTight3_2 = tau2.tauID("againstMuonTight3"); 	 
  chargedIsoPtSum_2   = tau2.tauID("chargedIsoPtSum");
  neutralIsoPtSum_2   = tau2.tauID("neutralIsoPtSum");
  puCorrPtSum_2       = tau2.tauID("puCorrPtSum");
  iso_2 = chargedIsoPtSum_2 + max(0., neutralIsoPtSum_2 - 0.2* puCorrPtSum_2);	
  pat::PackedCandidate const* packedLeadTauCand2 = dynamic_cast<pat::PackedCandidate const*>(tau2.leadChargedHadrCand().get());
  d0_2   = packedLeadTauCand2->dxy();
  dZ_2   = packedLeadTauCand2->dz();
}// settau


/////////////////////////////////
//			       //
// fill tree for tau variables //
//			       //
/////////////////////////////////
void miniAnalyzer::settau(const pat::Tau &tau){
  //cout<<"settau(tau): \n";
  decayModeFindingOldDMs_2 = tau.tauID("decayModeFinding"); 	 
  byCombinedIsolationDeltaBetaCorrRaw3Hits_2    = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"); 	
  byLooseCombinedIsolationDeltaBetaCorr3Hits_2  = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
  byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
  byTightCombinedIsolationDeltaBetaCorr3Hits_2  = tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
  byVLooseIsolationMVArun2v1DBoldDMwLT_2 = tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
  byLooseIsolationMVArun2v1DBoldDMwLT_2  = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
  byMediumIsolationMVArun2v1DBoldDMwLT_2 = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
  byTightIsolationMVArun2v1DBoldDMwLT_2  = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
  byVTightIsolationMVArun2v1DBoldDMwLT_2 = tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
  againstElectronLooseMVA6_2  = tau.tauID("againstElectronLooseMVA6"); 	 
  againstElectronMediumMVA6_2 = tau.tauID("againstElectronMediumMVA6"); 	 
  againstElectronTightMVA6_2  = tau.tauID("againstElectronTightMVA6"); 	 
  againstElectronVLooseMVA6_2 = tau.tauID("againstElectronVLooseMVA6"); 	 
  againstElectronVTightMVA6_2 = tau.tauID("againstElectronVTightMVA6"); 	 
  againstMuonLoose3_2 = tau.tauID("againstMuonLoose3"); 	 
  againstMuonTight3_2 = tau.tauID("againstMuonTight3"); 	 
  chargedIsoPtSum_2   = tau.tauID("chargedIsoPtSum");
  neutralIsoPtSum_2   = tau.tauID("neutralIsoPtSum");
  puCorrPtSum_2       = tau.tauID("puCorrPtSum");
  iso_2 = chargedIsoPtSum_2 + max(0., neutralIsoPtSum_2 - 0.2* puCorrPtSum_2);	
  pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
  d0_2   = packedLeadTauCand->dxy();
  dZ_2   = packedLeadTauCand->dz();
}// settau


/////////////////////////////////////
//			           //
// fill tree for general variables //
//			           //
/////////////////////////////////////
template<class L1, class L2> 
  void 
miniAnalyzer::setvar(
    const multimap<L1, L2> &ml1l2, 
    const reco::Vertex     &pv,
    const pat::MET 	   &pfmet, 
    const pat::MET 	   &puppiMET, 
    const TMatrixD         &covMET)
{
  // variables to be skip during variable tree cleaning
  vector<string> vskip{"weight", "evt", "lumi", "run", "rho", "npu", "npv", "NUP",
    "met", "metphi", "metcov00", "metcov01", "metcov10", "metcov11",
    "mvamet", "mvametphi", "mvacov00", "mvacov01", "mvacov10", "mvacov11",
    "puppimet", "puppimetphi"
  };
  // set var to -10000
  for(auto it=mdouble.begin(); it!=mdouble.end(); ++it){ if(find(vskip.begin(), vskip.end(), it->first) != vskip.end()) continue; *(it->second)=-10000; }
  for(auto it=mint.begin();    it!=mint.end();    ++it){ if(find(vskip.begin(), vskip.end(), it->first) != vskip.end()) continue; *(it->second)=-10000; }

  // retrieve selected pair
  decltype(*ml1l2.begin()->first)  l1{*ml1l2.begin()->first};
  decltype(*ml1l2.begin()->second) l2{*ml1l2.begin()->second};

  dr_tt   = deltaR(l1.p4(), l2.p4());
  dphi_tt = deltaPhi(l1.phi(), l2.phi());

  // leg 1 (leading tau for tt, electon for et, em muon for mt) 
  pt_1   = l1.pt(); 	 
  phi_1  = l1.phi(); 	 
  eta_1  = l1.eta(); 	 
  m_1    = l1.mass();
  q_1    = l1.charge();
  //mt_1 = ??
  pfmt_1    = sqrt(2*l1.pt()*met*(1-cos(deltaPhi(l1.phi(), metphi)))); 
  puppimt_1 = sqrt(2*l1.pt()*puppimet*(1-cos(deltaPhi(l1.phi(), puppimetphi))));
  if(decay!="tt"){d0_1 = l1.bestTrack()->dxy(pv.position()); dZ_1 = l1.bestTrack()->dz(pv.position()); iso_1 = leptoniso(l1);} 

  //Leg 2 (trailing tau for tt, tau for et, mt, muon for em)  	 
  pt_2   = l2.pt(); 	 
  phi_2  = l2.phi(); 	 
  eta_2  = l2.eta(); 	 
  m_2    = l2.mass(); 	 
  q_2    = l2.charge();
  //mt_2 = ??
  pfmt_2 = sqrt(2*l2.pt()*pfmet.pt()*(1-cos(deltaPhi(l2.phi(), pfmet.phi())))); 
  puppimt_2 = sqrt(2*l2.pt()*puppimet*(1-cos(deltaPhi(l2.phi(), puppimetphi))));
  if(decay=="em"){d0_2 = l2.bestTrack()->dxy(pv.position()); dZ_2 = l2.bestTrack()->dz(pv.position()); iso_2  = leptoniso(l2);}

  //di-tau system 
  pt_tt  = l1.pt()+l2.pt()+pfmet.pt();  	     

  mt_tot = sqrt(2*l1 .pt()*pfmet.pt()*(1-cos(deltaPhi(l1.phi(), pfmet.phi()))) + 
      2*l2.pt()*pfmet.pt()*(1-cos(deltaPhi(l2.phi(), pfmet.phi()))) +
      2*l1.pt()*l2.pt()   *(1-cos(deltaPhi(l1.phi(), l2.phi()))));  

  m_vis  = (l1.p4()+l2.p4()).mass();

  svfit(l1.p4(), l2.p4(), pfmet.p4(), covMET);

  pzetavis    = (1/sqrt(pow(cos(l1.phi())+cos(l2.phi()), 2) + pow(sin(l1.phi())+sin(l2.phi()), 2))) * 
    ((l1.px()+l2.px()) * (cos(l1.phi())+cos(l2.phi())) + (sin(l1.phi())+cos(l2.phi())) * (l1.py()+l2.py())); 

  //pzetamiss = (1/sqrt(pow(cos(l1.phi())+cos(l2.phi()), 2) + pow(sin(l1.phi())+sin(l2.phi()), 2))) * 
  //  (mvamet.px() * (cos(l1.phi())+cos(l2.phi())) + mvamet.py() * (sin(l1.phi())+sin(l2.phi()))); 

  pfpzetamiss = (1/sqrt(pow(cos(l1.phi())+cos(l2.phi()), 2) + pow(sin(l1.phi())+sin(l2.phi()), 2))) * 
    (pfmet.px() * (cos(l1.phi())+cos(l2.phi())) + pfmet.py() * (sin(l1.phi())+sin(l2.phi()))); 

  puppipzetamiss = (1/sqrt(pow(cos(l1.phi())+cos(l2.phi()), 2) + pow(sin(l1.phi())+sin(l2.phi()), 2))) * 
    (puppiMET.px() * (cos(l1.phi())+cos(l2.phi())) + puppiMET.py() * (sin(l1.phi())+sin(l2.phi()))); 

}// setvar



///////////////////////////////////////////////////////
//					             //	
// return true if the evt pass the preselection cuts //
//					             //	
///////////////////////////////////////////////////////
bool  miniAnalyzer::passPostSyncCut()
{
  //cout<<"passPostSyncCut(): \n";

  if(nbtag!=0) return false;

  if (decay=="mt")
    return 
      iso_1<0.15 && 
      againstElectronVLooseMVA6_2>0.5 &&
      againstMuonTight3_2>0.5 &&
      byLooseCombinedIsolationDeltaBetaCorr3Hits_2>0.5 && 
      dimuon_veto<1 && 
      extraelec_veto<1 && 
      extramuon_veto<1;

  else if (decay=="et")
    return 
      iso_1<0.1 && 
      againstElectronTightMVA6_2>0.5 && 
      againstMuonLoose3_2>0.5 && 
      byTightCombinedIsolationDeltaBetaCorr3Hits_2>0.5 && 
      dielectron_veto<1 && 
      extraelec_veto<1 && 
      extramuon_veto<1;

  else if(decay=="tt") 
    return 
      againstElectronVLooseMVA6_1>0.5 && 
      againstMuonLoose3_1>0.5 &&
      byTightCombinedIsolationDeltaBetaCorr3Hits_1>0.5 && 
      againstElectronVLooseMVA6_2>0.5 && 
      againstMuonLoose3_2>0.5 && 
      byTightCombinedIsolationDeltaBetaCorr3Hits_2>0.5 && 
      extraelec_veto<1 && 
      extramuon_veto<1;

  else if(decay=="em") 
    return 
      iso_1<0.15 && 
      iso_2<0.2 &&  
      extraelec_veto<1 && 
      extramuon_veto<1;

  else { cout<<"passPostSyncCut(): no decay return 0!!"; return 0;}
}


//////////////////////////////////
//				//	
// set the selected pair lepton //
// 				//
//////////////////////////////////
template <class L1, class L2> void selpair( multimap<L1, L2> &mapl1l2){
  double iso, miniso{1000}, pt, maxpt{-1000};

  //cout<<"\nselpair(): <<endl; 
  if(mapl1l2.size()<=1) return;

  //find min iso l1
  for(auto it = mapl1l2.begin(); it != mapl1l2.end(); it++){
    iso = leptoniso(*it->first);
    if(iso<miniso) miniso=iso;
  }
  //erase key in map with iso> min iso 
  for(auto it = mapl1l2.begin(); it != mapl1l2.end(); ){
    iso = leptoniso(*it->first);
    if(iso!=miniso) it = mapl1l2.erase(it);
    else ++it;
  }

  if(mapl1l2.size()>=2){

    //find max pt l1 
    for(auto it = mapl1l2.begin(); it != mapl1l2.end(); it++){
      pt = it->first->pt();
      if(pt>maxpt) maxpt=pt;
    }
    //erase key in map with pt< max pt  
    for(auto it = mapl1l2.begin(); it != mapl1l2.end(); ){
      pt = it->first->pt();
      if(pt!=maxpt) it = mapl1l2.erase(it);
      else ++it;
    }
  }

  if(mapl1l2.size()>=2){
    //find min iso l2
    miniso=1000;
    for(auto it = mapl1l2.begin(); it != mapl1l2.end(); it++){
      iso = leptoniso(*it->second);
      if(iso<miniso) miniso=iso;
    }
    //erase key in map with iso> max iso 
    for(auto it = mapl1l2.begin(); it != mapl1l2.end(); ){
      iso = leptoniso(*(it->second));
      if(iso!=miniso) it = mapl1l2.erase(it);
      else ++it;
    }
  }

  if(mapl1l2.size()>=2){
    maxpt=-1000;
    //find max pt l2 
    for(auto it = mapl1l2.begin(); it != mapl1l2.end(); it++){
      pt = it->second->pt();
      if(pt>maxpt) maxpt=pt;
    }
    //erase key in map with pt< maxpt  
    for(auto it = mapl1l2.begin(); it != mapl1l2.end(); ){
      pt = it->second->pt();
      if(pt!=maxpt) it = mapl1l2.erase(it);
      else ++it;
    }
  }
} // selpair


////////////////////////////////////////////////////
//						  //	
// overloading function, return  lepton isolation //
// 						  //
////////////////////////////////////////////////////
double leptoniso(const pat::Tau &tau){
  return -tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw"); //tauiso better for hight MVA value!!
}
double leptoniso(const pat::Electron &el){
  return (el.pfIsolationVariables().sumChargedHadronPt +
      max(0.0, el.pfIsolationVariables().sumNeutralHadronEt +
	el.pfIsolationVariables().sumPhotonEt - 
	0.5 * el.pfIsolationVariables().sumPUPt)) / el.pt();
}
double leptoniso(const pat::Muon &mu){
  return (mu.pfIsolationR04().sumChargedHadronPt + 
      std::max(0., (mu.pfIsolationR04().sumNeutralHadronEt + 
	  mu.pfIsolationR04().sumPhotonEt - 
	  0.5*mu.pfIsolationR04().sumPUPt ))) / mu.pt(); 
}


//////////////////////////////////////////////////
//						//	
// return true if a trigger match to the lepton // 
// 				          	//
//////////////////////////////////////////////////
bool triggermatching(const edm::TriggerNames &triggerNames, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, LorentzVector lp4, vector<string> &vpathLabel) {

  for(auto &pathLabel : vpathLabel){
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
      //cout<<"triggermatching(): unpackPathNames... \n";
      obj.unpackPathNames(triggerNames); 
      // Print trigger path
      vector<string> vpath{obj.pathNames()};
      //for(auto &i : vpath) cout<<"triggermatching(): pathNames "<<i<<endl;
      if (obj.hasPathName(pathLabel)) {
	//cout<<"triggermatching(): pathmatch :"<<pathLabel<<endl;
	//cout<<"obj.pdgid: "<<obj.pdgId()<<endl;
	double dR = deltaR(lp4, obj.p4());
	if (dR < 0.5) return true;
      }
    }
  }
  return true; // should be false to apply the trigger matching
}


////////////////////////////////
//			      //	
// set extra mu/ele veto flag //
//			      //	
////////////////////////////////
template<class L1, class L2> 
void miniAnalyzer::setthirdleptonveto(const multimap<L1, L2> &ml1l2, const vector<pat::Muon> &vmu, const vector<pat::Electron> &vel, const reco::Vertex &pv){

  //decltype(*ml1l2.begin()->first)  l1{*ml1l2.begin()->first};
  //decltype(*ml1l2.begin()->second) l2{*ml1l2.begin()->second};

  int nm{0}, ne{0};
  for(auto &mu : vmu){
    if(mu.pt()<=10)continue;
    if(fabs(mu.eta())>=2.4) continue;
    if(mu.bestTrack()->dxy(pv.position()) >=0.045 ) continue; 		 
    if(mu.bestTrack()->dz (pv.position()) >=0.2) continue; 
    if(!isMediumMuon(mu)) continue;
    if(leptoniso(mu)>=0.3) continue;
    nm++;
  }

  for(auto &el : vel){
    if(el.pt()<=10)continue;
    if(fabs(el.eta())>=2.5) continue;
    if(el.bestTrack()->dxy(pv.position()) >=0.045 ) continue; 		 
    if(el.bestTrack()->dz (pv.position()) >=0.2) continue; 

    // map btw vel and map(vel, mva)
    decltype (*melmva.find(&el)->first) mel {*melmva.find(&el)->first};
    double mmva {melmva.find(&el)->second};
    if      ( fabs(mel.superCluster()->eta())<0.8   && 			                	mmva <0.913286 ) continue; //barrel 1
    else if ( fabs(mel.superCluster()->eta())>0.8   && fabs(mel.superCluster()->eta())<1.479 && mmva <0.805013 ) continue; //barrel 2
    else if ( fabs(mel.superCluster()->eta())>1.479 && fabs(mel.superCluster()->eta())<3     && mmva <0.358969 ) continue; //endcap

    if(!el.passConversionVeto()) continue;
    if(el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) >1) continue;
    if(leptoniso(el)>=0.3) continue;
    ne++;
  }

  if(decay=="mt"){
    if(nm>=2) extramuon_veto=1;
    if(ne!=0) extraelec_veto=1;
  }
  if(decay=="et"){ 
    if(ne>=2) extraelec_veto=1;
    if(nm!=0) extramuon_veto=1;
  }
  if(decay=="em"){
    if(nm>=2) extramuon_veto=1;
    if(ne>=2) extraelec_veto=1;
  }
  if(decay=="tt"){
    if(nm!=0) extramuon_veto=1;
    if(ne!=0) extraelec_veto=1;
  }
}

//////////////////////////
//			//	
// set dimuon veto flag //
//			//	
//////////////////////////
void miniAnalyzer::setdileptonveto(const vector< pat::Muon> &vmu, const reco::Vertex &pv)
{
  //cout<<"dimuveto()"<<endl;
  for(const auto &mu1:vmu){
    if(mu1.pt() <= 15)  continue; 
    if(fabs(mu1.eta()) >= 2.4) continue;
    if(!mu1.isGlobalMuon())continue;    
    if(!mu1.isTrackerMuon())continue;  
    if(!mu1.isPFMuon())  continue;      
    if(fabs(mu1.bestTrack()->dz(pv.position()))  >= 0.2) continue;  
    if(fabs(mu1.bestTrack()->dxy(pv.position())) >= 0.045) continue;
    if(leptoniso(mu1) >= 0.3) continue;
    for(const auto &mu2:vmu){
      if(deltaR(mu1.p4(), mu2.p4())<=0.15) continue;
      if(mu2.charge()*mu1.charge()>=0) continue;
      if(mu2.pt() <= 15)  continue; 
      if(fabs(mu2.eta()) >= 2.4) continue;
      if(!mu2.isGlobalMuon())continue;    
      if(!mu2.isTrackerMuon())continue;  
      if(!mu2.isPFMuon())  continue;      
      if(fabs(mu2.bestTrack()->dz(pv.position()))  >= 0.2) continue;  
      if(fabs(mu2.bestTrack()->dxy(pv.position())) >= 0.045) continue;
      if(leptoniso(mu2) >= 0.3) continue;
      dimuon_veto=1; // if found a couple, set the flag to true
    }
  }
}


//////////////////////////
//			//	
// set dielec veto flag //
//			//	
//////////////////////////
void miniAnalyzer::setdileptonveto(const vector<pat::Electron> &vel, const reco::Vertex &pv)
{
  //cout<<"dielveto()"<<endl;
  for(const auto &el1:vel){
    if(el1.pt() <= 15) continue;  
    if(fabs(el1.eta()) >= 2.5 ) continue;
    if(!spring15_25ns_cut_veto(el1, pv)) continue;
    if(fabs(el1.bestTrack()->dz(pv.position()))  >= 0.2) continue;
    if(fabs(el1.bestTrack()->dxy(pv.position())) >= 0.045) continue;
    if(leptoniso(el1) >= 0.3) continue;
    for(const auto &el2:vel){
      if(deltaR(el2.p4(), el1.p4())<=0.15) continue;
      if(el2.charge()*el1.charge()>0) continue;
      if(el2.pt() <= 15) continue;  
      if(fabs(el2.eta()) >= 2.5 ) continue;
      if(!spring15_25ns_cut_veto(el2, pv)) continue;
      if(fabs(el2.bestTrack()->dz(pv.position()))  >= 0.2) continue;
      if(fabs(el2.bestTrack()->dxy(pv.position())) >= 0.045) continue;
      if(leptoniso(el2) >= 0.3) continue;
      dielectron_veto=1; 
    } 
  } 
}


//////////////////////////////////////////////////////////////////////
//								    //	
// return true if the electron pass the spring15_25ns_cutbased_veto //
//								    //	
//////////////////////////////////////////////////////////////////////
bool spring15_25ns_cut_veto(const pat::Electron &el, const reco::Vertex &pv){
  // POG Spring15 25ns cut-based "Veto" ID 
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_8_0
  bool spring15Veto;
  if(fabs(el.superCluster()->eta())<=1.479){
    spring15Veto = (
	el.full5x5_sigmaIetaIeta()                <  0.0114 && 
	fabs(el.deltaEtaSuperClusterTrackAtVtx()) <  0.0152 &&
	fabs(el.deltaPhiSuperClusterTrackAtVtx()) <  0.216  &&
	el.hadronicOverEm()                       <  0.181  && 
	//relIsoWithEA  			    <  0.161  && \\ remove for consistency
	fabs(1/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy()) <  0.207  && // ooEmooP
	//fabs(el.bestTrack()->dxy(pv.position()))  <  0.0564 && \\ remove for consistency
	//fabs(el.bestTrack()->dz(pv.position()))   <  0.472  && \\ remove for consistency
	el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <=2 &&
	el.passConversionVeto());
  }
  else{
    spring15Veto = (
	el.full5x5_sigmaIetaIeta()                <  0.0352 && 
	fabs(el.deltaEtaSuperClusterTrackAtVtx()) <  0.0113 &&
	fabs(el.deltaPhiSuperClusterTrackAtVtx()) <  0.237  &&
	el.hadronicOverEm()                       <  0.116  && 
	//relIsoWithEA  			    <  0.193  && \\ remove for consistency
	fabs(1/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy()) <  0.174 && // ooEmooP
	//fabs(el.bestTrack()->dxy(pv.position()))  <  0.222  && \\ remove for consistency
	//fabs(el.bestTrack()->dz(pv.position()))   <  0.921  && \\ remove for consistency
	el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <=3 &&
	el.passConversionVeto());
  }
  return spring15Veto;
}


///////////////////////////////////////////////////
//						 //
// return true if the muon passes the medium cut //
//						 //
///////////////////////////////////////////////////
bool isMediumMuon(const pat::Muon & recoMu) 
{
  bool goodGlob = 
    recoMu.isGlobalMuon()                           && 
    recoMu.globalTrack()->normalizedChi2() < 3      && 
    recoMu.combinedQuality().chi2LocalPosition < 12 && 
    recoMu.combinedQuality().trkKink < 20; 
  bool isMedium = 
    muon::isLooseMuon(recoMu) && 
    recoMu.innerTrack()->validFraction() > 0.49     && 
    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}


//////////////////////////////////////////////////
//	   				        //
// return true if the jet passes the medium cut //
//					        //
//////////////////////////////////////////////////
bool isLooseJet(const pat::Jet & jet){
  double NHF	            = jet.neutralHadronEnergyFraction();
  double NEMF		    = jet.neutralEmEnergyFraction();
  double CHF	            = jet.chargedHadronEnergyFraction();
  double CEMF  	            = jet.chargedEmEnergyFraction();
  double NumConst           = jet.chargedMultiplicity()+jet.neutralMultiplicity();
  double NumNeutralParticle = jet.neutralMultiplicity();
  double CHM 		    = jet.chargedMultiplicity(); 

  bool looseJet = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(jet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(jet.eta())>2.4) && fabs(jet.eta())<=2.7) ||
    (NEMF<0.90 && NumNeutralParticle>2 && fabs(jet.eta())>2.7 && fabs(jet.eta())<=3.0) ||
    (NEMF<0.90 && NumNeutralParticle>10 && fabs(jet.eta())>3.0);

  return looseJet;
} 


/////////////////////////////////
//			       //
// return the gentau particule //
//			       //
/////////////////////////////////
const reco::GenParticle* getGenTau(const pat::Tau& patTau)
{
  vector<reco::GenParticleRef> genParticles = patTau.genParticleRefs();
  for ( auto &gen : genParticles){
    if ( gen.isAvailable() ) {
      const reco::GenParticleRef& genParticle = (gen);
      if ( abs(genParticle->pdgId())==15 ) return genParticle.get();
    }
  }
  return 0;
}


/////////////////////////////
//			   //
// set the svfit variables //
//			   //
/////////////////////////////
void miniAnalyzer::svfit(
    const LorentzVector &l1, 
    const LorentzVector &l2, 
    const LorentzVector &met,
    const TMatrixD      &covMET)
{
  svFitStandalone::kDecayType l1Type, l2Type;
  if ( decay=="mt" ) {
    l1Type = svFitStandalone::kTauToMuDecay;
    l2Type = svFitStandalone::kTauToHadDecay;
  } 
  else if ( decay=="et" ) {
    l1Type = svFitStandalone::kTauToElecDecay;
    l2Type = svFitStandalone::kTauToHadDecay;
  } 
  else if ( decay=="tt" ) {
    l1Type = svFitStandalone::kTauToHadDecay;
    l2Type = svFitStandalone::kTauToHadDecay;
  } 
  else if ( decay=="mm" ) {
    l1Type = svFitStandalone::kTauToMuDecay;
    l2Type = svFitStandalone::kTauToMuDecay;
  } 
  else if ( decay=="ee" ) {
    l1Type = svFitStandalone::kTauToElecDecay;
    l2Type = svFitStandalone::kTauToElecDecay;
  } 
  else if ( decay=="em" ) {
    l1Type = svFitStandalone::kTauToElecDecay;
    l2Type = svFitStandalone::kTauToMuDecay;
  } 
  else {
    cerr << "svfit(): invalid channel!!" << endl;
    assert(0);
  }

  // define lepton four vectors
  vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l1Type, l1.Pt(), l1.Eta(), l1.Phi(), l1.M())); 
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(l2Type, l2.Pt(), l2.Eta(), l2.Phi(), l2.M()));
  SVfitStandaloneAlgorithm algo(measuredTauLeptons, met.Px(), met.Py(), covMET, 0);
  algo.addLogM(false);  

  edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);  
  TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
  algo.shiftVisPt(true, inputFile_visPtResolution);
  algo.integrateMarkovChain();

  m_sv   = algo.getMass(); // Full SVFit mass - return value is in units of GeV
  mt_sv  = algo.transverseMass(); // Transverse SVFit mass
  pt_sv  = algo.pt();
  eta_sv = algo.eta();
  phi_sv = algo.phi();
  met_sv = algo.fittedMET().Rho();

  //if ( algo.isValidSolution() ) cout << "svfit(): m/mt= " << m_sv << "/"<< mt_sv << endl;
  //else cout << "svfit(): status of NLL is not valid [" << algo.isValidSolution() << "]" << endl;

  delete inputFile_visPtResolution;
}// svfit()


// ------------ method called once each job just before starting event loop  ------------
  void 
miniAnalyzer::beginJob()
{

  edm::Service<TFileService> fs;
  TFileDirectory subdirsync = fs->mkdir("sync");
  nevt_sync   = subdirsync.make<TH1I>("nevt", "Initial number of event", 1, 0, 1);
  mutau_sync  = new TTree("mt", "mt");
  etau_sync   = new TTree("et", "et");
  tautau_sync = new TTree("tt", "tt");
  emu_sync    = new TTree("em", "em");

  TFileDirectory subdirh2taus = fs->mkdir("h2taus");
  nevt   = subdirh2taus.make<TH1I>("nevt", "Initial number of event", 1, 0, 1);
  mutau  = new TTree("mt", "mt");
  etau   = new TTree("et", "et");
  tautau = new TTree("tt", "tt");
  emu    = new TTree("em", "em");

  for(auto it=mint.begin(); it!=mint.end(); ++it) {
    mutau ->Branch((it->first).c_str(), it->second);
    etau  ->Branch((it->first).c_str(), it->second);
    tautau->Branch((it->first).c_str(), it->second);
    emu   ->Branch((it->first).c_str(), it->second);
  }
  for(auto it=mdouble.begin(); it!=mdouble.end(); ++it) {
    mutau ->Branch((it->first).c_str(), it->second);
    etau  ->Branch((it->first).c_str(), it->second);
    tautau->Branch((it->first).c_str(), it->second);
    emu   ->Branch((it->first).c_str(), it->second);
  }

  for(auto it=mint.begin(); it!=mint.end(); ++it) {
    mutau_sync ->Branch((it->first).c_str(), it->second);
    etau_sync  ->Branch((it->first).c_str(), it->second);
    tautau_sync->Branch((it->first).c_str(), it->second);
    emu_sync   ->Branch((it->first).c_str(), it->second);
  }
  for(auto it=mdouble.begin(); it!=mdouble.end(); ++it) {
    mutau_sync ->Branch((it->first).c_str(), it->second);
    etau_sync  ->Branch((it->first).c_str(), it->second);
    tautau_sync->Branch((it->first).c_str(), it->second);
    emu_sync   ->Branch((it->first).c_str(), it->second);
  }
}


// ------------ method called after each end of run  ------------

void miniAnalyzer::beginRun(const edm::Run &run, const edm::EventSetup &es)
{
}


void miniAnalyzer::endRun(const edm::Run &run, const edm::EventSetup &es)
{
  // event cs
  edm::Handle<GenRunInfoProduct> runinfo; run.getByToken(runinfoToken_, runinfo);
  weight = runinfo->crossSection();
  //cout<<"cross section: "<<weight<<endl;
  //cout<<"externalXSecLO (): "<<runinfo->externalXSecLO().value()<<endl;
  //cout<<"externalXSecNLO (): "<<runinfo->externalXSecNLO().value()<<endl;
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
miniAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
miniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniAnalyzer);



/*
int  miniAnalyzer::setgen(const edm::Event& iEvent,const reco::Candidate &cand){

  cout<<"setgen()\n";
  Handle<edm::View<reco::GenParticle> > pruned; iEvent.getByToken(prunedGenToken_, pruned);  
  Handle<reco::GenJetCollection> gentaul;       iEvent.getByToken(genTaulToken_,   gentaul);
  Handle<reco::GenJetCollection> gentauh; 	iEvent.getByToken(genTauhToken_,   gentauh);

  cout<<"prun size: "<< pruned ->size() <<endl;
  cout<<"taul size: "<< gentaul->size()<<endl;
  cout<<"tauh size: "<< gentauh->size()<<endl;

  int gid{6}; // return 6 if failed all other
  double tmpdr{1000}; // test deltaR

  // loop over gen part and select the closest
  for(const auto &gen : *pruned){
    double dr{deltaR(gen.p4(), cand.p4())};
    if(dr       >= 0.2) continue;
    if(gen.pt() <=   8) continue;
    if(!gen.statusFlags().isPrompt()) continue;

    if     (dr<tmpdr && abs(gen.pdgId())==11) { tmpdr=dr; gid=1; }
    else if(dr<tmpdr && abs(gen.pdgId())==13) { tmpdr=dr; gid=2; }
  }
  //cout<<"gid1: "<<gid<<endl;

  // loop over lepton tau
  for(const auto &gen : *gentaul){

    // it is a get, so we need to look at it's constituants
    std::vector< const GenParticle * > vgp{gen.getGenConstituents()};
    cout<<"gentaul size(): "<<vgp.size()<<endl;

    for(auto &gencons : vgp) {
      double dr{deltaR(gencons->p4(), cand.p4())};
      if(dr            >= 0.2) continue;
      if(gencons->pt() <=   8) continue;
      if(!gencons->statusFlags().isDirectPromptTauDecayProduct()) continue;

      if     (dr<tmpdr && abs(gencons->pdgId())==11) { tmpdr=dr; gid=3; }
      else if(dr<tmpdr && abs(gencons->pdgId())==13) { tmpdr=dr; gid=4; }
    }
  }

  // loop over had tau
  for(const auto &gen : *gentauh){

    double dr{deltaR(gen.p4(), cand.p4())};
    if(dr >= 0.2) continue;
    if(dr>tmpdr) continue;
    tmpdr=dr;

    std::vector< const GenParticle * > vgp{gen.getGenConstituents()};
    //cout<<"gentauhconst size(): "<<vgp.size()<<endl;

    double sumpt{0};
    for(auto &gencons : vgp) {
      if(abs(gencons->pdgId()==11) || abs(gencons->pdgId()==13)) continue;
      cout<<"id: "<<gencons->pdgId()<<endl;
      cout<<"pt: "<<gencons->pt()<<endl;
      cout<<"flag promt: "<<gencons->statusFlags().isPrompt()<<endl;
      cout<<"flag tau: "<<gencons->statusFlags().isDirectPromptTauDecayProduct()<<endl;
      //if(!gencons->statusFlags().isPrompt()) continue;
      if(!gencons->statusFlags().isDirectPromptTauDecayProduct()) continue;
      cout<<"pass flag!!\n";
      sumpt+=gencons->pt();
    }
    cout<<"sumpt: "<<sumpt<<endl;
    if (sumpt >15 && dr<tmpdr) { tmpdr=dr; gid=5; }
  }
  cout<<"gid: "<<gid<<endl;

  return gid;
}
*/
