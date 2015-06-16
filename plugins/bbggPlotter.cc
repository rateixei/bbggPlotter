// -*- C++ -*-
//
// Package:    flashgg/bbggPlotter
// Class:      bbggPlotter
// 
/**\class bbggPlotter bbggPlotter.cc flashgg/bbggPlotter/plugins/bbggPlotter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rafael Teixeira De Lima
//         Created:  Tue, 16 Jun 2015 17:11:20 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class bbggPlotter : public edm::EDAnalyzer {
   public:
      explicit bbggPlotter(const edm::ParameterSet&);
      ~bbggPlotter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      //Parameter tokens
      EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
      EDGetTokenT<View<Photon> > photonToken_;
      EDGetTokenT<View<Jet> > thejetToken_;
      EDGetTokenT<View<reco::Vertex> > vertexToken_;

      //Thresholds
      vector<double> ph_pt;
      vector<double> ph_eta;
      vector<double> ph_hoe;
      vector<double> ph_sieie;
      vector<double> ph_r9;
      vector<double> ph_chIso;
      vector<double> ph_nhIso;
      vector<double> ph_phIso;
      vector<double> ph_elVeto;
      vector<bool> ph_doID;

      vector<double> diph_pt;
      vector<double> diph_eta;
      vector<double> diph_mass;

      vector<double> jt_pt;
      vector<double> jt_eta;
      vector<double> jt_bDis;
      vector<bool> jt_doPU;

      vector<double> dijt_pt;
      vector<double> dijt_eta;
      vector<double> dijt_mass; 

      vector<double> cand_mass;

      //OutFile & Hists
      TFile* outFile;
      std::map<std::string, TH1F*> hists;

};

bbggPlotter::bbggPlotter(const edm::ParameterSet& iConfig)
diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getUntrackedParameter<InputTag> ( "DiPhotonTag", InputTag( "flashggDiPhotons" ) ) ) ),
photonToken_( consumes<View<flashgg::Photon> >( iConfig.getUntrackedParameter<InputTag>( "PhotonTag", InputTag( "flashggPhotons" ) ) ) ),
thejetToken_( consumes<View<flashgg::Jet> >( iConfig.getUntrackedParameter<InputTag>( "JetTag", InputTag( "flashggJets" ) ) ) ),
vertexToken_( consumes<View<reco::Vertex> >( iConfig.getUntrackedParameter<InputTag> ( "VertexTag", InputTag( "offlinePrimaryVertices" ) ) ) )

{
   //now do what ever initialization is needed
//Default values for thresholds
      vector<double> def_ph_pt;
      vector<double> def_ph_eta;
      vector<double> def_ph_hoe;
      vector<double> def_ph_sieie;
      vector<double> def_ph_r9;
      vector<double> def_ph_chIso;
      vector<double> def_ph_nhIso;
      vector<double> def_ph_phIso;
      vector<double> def_ph_elVeto;
      vector<bool> def_ph_doID;

      vector<double> def_diph_pt;
      vector<double> def_diph_eta;
      vector<double> def_diph_mass;

      vector<double> def_jt_pt;
      vector<double> def_jt_eta;
      vector<double> def_jt_bDis;
      vector<bool> def_jt_doPU;

      vector<double> def_dijt_pt;
      vector<double> def_dijt_eta;
      vector<double> def_dijt_mass;

      vector<double> def_cand_mass;

      def_ph_pt.push_back(10.);         def_ph_pt.push_back(10.);
      def_ph_eta.push_back(0.);         def_ph_eta.push_back(0.);
      def_ph_hoe.push_back(-1.);        def_ph_hoe.push_back(-1.);
      def_ph_sieie.push_back(-1.);      def_ph_sieie.push_back(-1.);
      def_ph_r9.push_back(-1.);         def_ph_r9.push_back(-1.);
      def_ph_chIso.push_back(-1.);      def_ph_chIso.push_back(-1.);
      def_ph_nhIso.push_back(-1.);      def_ph_nhIso.push_back(-1.);
      def_ph_phIso.push_back(-1.);      def_ph_phIso.push_back(-1.);
      def_ph_elVeto.push_back(-1.);     def_ph_elVeto.push_back(-1.);
      def_ph_doID.push_back(false);     def_ph_doID.push_back(false);

      def_diph_pt.push_back(10.);       def_diph_pt.push_back(10.);
      def_diph_eta.push_back(0.);       def_diph_eta.push_back(0.);
      def_diph_mass.push_back(0.);      def_diph_mass.push_back(1000.);
			
      def_jt_pt.push_back(10.);         def_jt_pt.push_back(10.);
      def_jt_eta.push_back(0.);         def_jt_eta.push_back(0.);
      def_jt_bDis.push_back(0.);        def_jt_bDis.push_back(0.);
      def_jt_doPU.push_back(false);     def_jt_doPU.push_back(false);

      def_dijt_pt.push_back(10.);       def_dijt_pt.push_back(10.);
      def_dijt_eta.push_back(0.);       def_dijt_eta.push_back(0.);
      def_dijt_mass.push_back(0.);      def_dijt_mass.push_back(1000.);

      def_cand_mass.push_back(0.);	def_cand_mass.push_back(2000.);

//Get thresholds from config file
      ph_pt     = iConfig.getUntrackedParameter<vector<double > >("PhotonPt", def_ph_pt);
      ph_eta    = iConfig.getUntrackedParameter<vector<double > >("PhotonEta", def_ph_eta);
      ph_hoe    = iConfig.getUntrackedParameter<vector<double > >("PhotonHoverE", def_ph_hoe);
      ph_sieie  = iConfig.getUntrackedParameter<vector<double > >("PhotonSieie", def_ph_sieie);
      ph_r9     = iConfig.getUntrackedParameter<vector<double > >("PhotonR9", def_ph_r9);
      ph_chIso  = iConfig.getUntrackedParameter<vector<double > >("PhotonChargedIso", def_ph_chIso);
      ph_nhIso  = iConfig.getUntrackedParameter<vector<double > >("PhotonNeutralIso", def_ph_nhIso);
      ph_phIso  = iConfig.getUntrackedParameter<vector<double > >("PhotonPhotonIso", def_ph_phIso);
      ph_elVeto = iConfig.getUntrackedParameter<vector<double > >("PhotonElectronVeto", def_ph_elVeto);
      ph_doID   = iConfig.getUntrackedParameter<vector<double > >("PhotonDoID", def_ph_doID),;

      diph_pt   = iConfig.getUntrackedParameter<vector<double > >("DiPhotonPt", def_diph_pt);
      diph_eta  = iConfig.getUntrackedParameter<vector<double > >("DiPhotonEta", def_diph_eta);
      diph_mass = iConfig.getUntrackedParameter<vector<double > >("DiPhotonMassWindow", def_diph_mass),;

      jt_pt     = iConfig.getUntrackedParameter<vector<double > >("JetPt", def_jt_pt);
      jt_eta    = iConfig.getUntrackedParameter<vector<double > >("JetEta", def_jt_eta);
      jt_bDis   = iConfig.getUntrackedParameter<vector<double > >("JetBDiscriminant", def_jt_bDis);
      jt_doPU   = iConfig.getUntrackedParameter<vector<double > >("JetDoPUID", def_jt_doPU),;

      dijt_pt   = iConfig.getUntrackedParameter<vector<double > >("DiJetPt", def_dijt_pt);
      dijt_eta  = iConfig.getUntrackedParameter<vector<double > >("DiJetEta", def_dijt_eta);
      dijt_mass = iConfig.getUntrackedParameter<vector<double > >("DiJetMassWindow", def_dijt_mass);

      cand_mass = iConfig.getUntrackedParameter<vector<double > >("CandidateMassWindow", def_cand_mass);

}


bbggPlotter::~bbggPlotter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
bbggPlotter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
bbggPlotter::beginJob()
{
	outFile = new TFile("out.root", "RECREATE");

	hists["dipho_pt"] 	= new TH1F("dipho_pt", "DiPhoton p_{T}; p_{T}(#gamma#gamma) (GeV); Events", 100, 0, 300);
	hists["dipho_eta"] 	= new TH1F("dipho_eta", "DiPhoton #eta; #eta(#gamma#gamma); Events", 100, -5, 5);
	hists["dipho_mass"] 	= new TH1F("dipho_mass", "DiPhoton Mass; M(#gamma#gamma); Events", 100, diph_mass[0], diph_mass[1]);

	hists["dijet_pt"] 	= new TH1F("dijet_pt", "DiJet p_{T}; p_{T}(jj) (GeV); Events", 100, 0, 300);
	hists["dijet_eta"] 	= new TH1F("dijet_eta", "DiJet #eta; #eta(jj); Events", 100, -5, 5);
	hists["dijet_mass"] 	= new TH1F("dijet_mass", "DiJet Mass; M(jj); Events", 100, dijt_mass[0], dijt_mass[1]);

	hists["cand_pt"] 	= new TH1F("cand_pt", "DiHiggs Candidate (jj#gamma#gamma) p_{T}; p_{T}(jj#gamma#gamma) (GeV); Events", 100, 0, 500);
	hists["cand_eta"] 	= new TH1F("cand_eta", "DiHiggs Candidate (jj#gamma#gamma) #eta; #eta(jj#gamma#gamma) (GeV); Events", 100, -5, 5);
	hists["cand_mass"] 	= new TH1F("cand_mass", "DiHiggs Candidate (jj#gamma#gamma) Mass; M(jj#gamma#gamma) (GeV); Events", 100, cand_mass[0], cand_mass[1]);

	hists["pho1_pt"] 	= new TH1F("pho1_pt", "Leading Photon p_{T}; p_{T}(leading #gamma) (GeV); Events", 100, 10, 150);
	hists["pho1_eta"] 	= new TH1F("pho1_eta", "Leading Photon #eta; #eta(leading #gamma); Events", 100, -5., 5.);
	hists["pho1_hoe"] 	= new TH1F("pho1_hoe", "Leading Photon H/E; H/E(leading #gamma); Events", 100, -0.01, 0.1);
	hists["pho1_sieie"] 	= new TH1F("pho1_sieie", "Leading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(leading #gamma); Events", 100, 0.0, 0.04);
	hists["pho1_r9"] 	= new TH1F("pho1_r9", "Leading Photon R9; R9(leading #gamma); Events", 100, 0, 1.1);
	hists["pho1_chiso"] 	= new TH1F("pho1_chiso", "Leading Photon Corrected Charged Isolation; Corrected Charged Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_nhiso"] 	= new TH1F("pho1_nhiso", "Leading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_phiso"] 	= new TH1F("pho1_phiso", "Leading Photon Corrected Photon Isolation; Corrected Photon Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_elveto"] 	= new TH1F("pho1_elveto", "Leading Photon Electron Veto; Electron Veto(leading #gamma); events", 8, -1, 3);

	hists["pho2_pt"] 	= new TH1F("pho2_pt", "SubLeading Photon p_{T}; p_{T}(subLeading #gamma) (GeV); Events", 100, 10, 150);
	hists["pho2_eta"] 	= new TH1F("pho2_eta", "SubLeading Photon #eta; #eta(subLeading #gamma); Events", 100, -5., 5.);
	hists["pho2_hoe"] 	= new TH1F("pho2_hoe", "SubLeading Photon H/E; H/E(subLeading #gamma); Events", 100, -0.01, 0.1);
	hists["pho2_sieie"] 	= new TH1F("pho2_sieie", "SubLeading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(subLeading #gamma); Events", 100, 0.0, 0.04);
	hists["pho2_r9"] 	= new TH1F("pho2_r9", "SubLeading Photon R9; R9(subLeading #gamma); Events", 100, 0, 1.1);
	hists["pho2_chiso"] 	= new TH1F("pho2_chiso", "SubLeading Photon Corrected Charged Isolation; Corrected Charged Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_nhiso"] 	= new TH1F("pho2_nhiso", "SubLeading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_phiso"] 	= new TH1F("pho2_phiso", "SubLeading Photon Corrected Photon Isolation; Corrected Photon Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_elveto"] 	= new TH1F("pho2_elveto", "SubLeading Photon Electron Veto; Electron Veto(subLeading #gamma); events", 8, -1, 3);


	hists["jet1_pt"] 	= new TH1F("jet1_pt", "Leading Jet p_{T}; p_{T}(leading jet) (GeV); Events", 100, 10, 150);
	hists["jet1_eta"] 	= new TH1F("jet1_eta", "Leading Jet #eta; #eta(leading jet); Events", 100, -5., 5.);
	hists["jet1_bDis"] 	= new TH1F("jet1_bDis", "Leading Jet b-Discriminant; b-Discriminant(leading jet); Events", 100, -0.01, 1.01);

	hists["jet2_pt"] 	= new TH1F("jet2_pt", "SubLeading Jet p_{T}; p_{T}(subLeading jet) (GeV); Events", 100, 10, 150);
	hists["jet2_eta"] 	= new TH1F("jet2_eta", "SubLeading Jet #eta; #eta(subLeading jet); Events", 100, -5., 5.);
	hists["jet2_bDis"] 	= new TH1F("jet2_bDis", "SubLeading Jet b-Discriminant; b-Discriminant(subLeading jet); Events", 100, -0.01, 1.01);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbggPlotter::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bbggPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbggPlotter);
