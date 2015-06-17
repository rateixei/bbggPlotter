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
#include <vector.h>
#include <map.h>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

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

      // ---------- extra methods... ---------------------
      double getCHisoToCutValue(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho, float rho);
      double getNHisoToCutValue(flashgg::Photon* pho, float rho);
      double getPHisoToCutValue(flashgg::Photon* pho, float rho);
      double getEA( float eta, int whichEA);


      // ----------member data ---------------------------
      //Parameter tokens
      EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
      EDGetTokenT<View<Jet> > thejetToken_;
      edm::InputTag rhoFixedGrid_
      std::string bTagType;

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
      int n_bJets;

      vector<double> dijt_pt;
      vector<double> dijt_eta;
      vector<double> dijt_mass; 

      vector<double> cand_pt;
      vector<double> cand_eta;
      vector<double> cand_mass;

      //OutFile & Hists
      TFile* outFile;
      std::map<std::string, TH1F*> hists;

};

bbggPlotter::bbggPlotter(const edm::ParameterSet& iConfig)
diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getUntrackedParameter<InputTag> ( "DiPhotonTag", InputTag( "flashggDiPhotons" ) ) ) ),
thejetToken_( consumes<View<flashgg::Jet> >( iConfig.getUntrackedParameter<InputTag>( "JetTag", InputTag( "flashggJets" ) ) ) )
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
      int def_n_bJets;

      vector<double> def_dijt_pt;
      vector<double> def_dijt_eta;
      vector<double> def_dijt_mass;

      vector<double> def_cand_pt;
      vector<double> def_cand_eta;
      vector<double> def_cand_mass;

      std::string def_bTagType;

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

      def_n_bJets = 0;

      def_dijt_pt.push_back(10.);       def_dijt_pt.push_back(10.);
      def_dijt_eta.push_back(0.);       def_dijt_eta.push_back(0.);
      def_dijt_mass.push_back(0.);      def_dijt_mass.push_back(1000.);

      def_cand_pt.push_back(0.);
      def_cand_eta.push_back(0.);
      def_cand_mass.push_back(0.);		def_cand_mass.push_back(2000.);


      def_bTagType = "pfCombinedInclusiveSecondaryVertexV2BJetTags";

//Get thresholds from config file
      ph_pt     = iConfig.getUntrackedParameter<vector<double > >("PhotonPtOverDiPhotonMass", def_ph_pt);
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

      jt_pt     = iConfig.getUntrackedParameter<vector<double > >("JetPtOverDiJetMass", def_jt_pt);
      jt_eta    = iConfig.getUntrackedParameter<vector<double > >("JetEta", def_jt_eta);
      jt_bDis   = iConfig.getUntrackedParameter<vector<double > >("JetBDiscriminant", def_jt_bDis);
      jt_doPU   = iConfig.getUntrackedParameter<vector<double > >("JetDoPUID", def_jt_doPU);

      n_bJets = iConfig.getUntrackedParameter<int>("n_bJets", def_n_bJets);

      dijt_pt   = iConfig.getUntrackedParameter<vector<double > >("DiJetPt", def_dijt_pt);
      dijt_eta  = iConfig.getUntrackedParameter<vector<double > >("DiJetEta", def_dijt_eta);
      dijt_mass = iConfig.getUntrackedParameter<vector<double > >("DiJetMassWindow", def_dijt_mass);

      cand_pt 	= iConfig.getUntrackedParameter<vector<double > >("CandidatePt", def_cand_mass);
      cand_eta 	= iConfig.getUntrackedParameter<vector<double > >("CandidateEta", def_cand_mass);
      cand_mass = iConfig.getUntrackedParameter<vector<double > >("CandidateMassWindow", def_cand_mass);

      rhoFixedGrid_  = iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" );

      bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", bTagType );

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
	 
   double dipho_pt = -1, dipho_eta = -1, dipho_phi = -1; dipho_mass = -1;
   
   double dijet_pt = -1, dijet_eta = -1, dipho_phi = -1; dijet_mass = -1;
   
   double cand4_pt = -1, cand4_eta = -1, cand4_phi = -1; cand4_mass = -1;
   
   double pho1_pt = -1, pho1_eta = -1, pho1_phi = -1;
   double pho1_hoe = -1, pho1_sieie = -1, pho1_r9 = -1, pho1_chiso = -1, pho1_nhiso = -1, pho1_phiso = -1, pho1_elveto = -1;
   
   double pho2_pt = -1, pho2_eta = -1, pho2_phi = -1;
   double pho2_hoe = -1, pho2_sieie = -1, pho2_r9 = -1, pho2_chiso = -1, pho2_nhiso = -1, pho2_phiso = -1, pho2_elveto = -1;
   
   double jet1_pt = -1, jet1_eta = -1, jet1_phi = -1; jet1_bDis = -1;
   double jet2_pt = -1, jet2_eta = -1, jet2_phi = -1; jet2_bDis = -1;

   Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
   evt.getByToken( diPhotonToken_, diPhotons );
   Handle<View<flashgg::Jet> > theJets;
   evt.getByToken( thejetToken_, theJets );
   Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
   evt.getByLabel( rhoFixedGrid_, rhoHandle );
   const double rhoFixedGrd = *( rhoHandle.product() );

   bool isValidDiPhotonCandidate = false;

   edm::Ptr<reco::Vertex> CandVtx;
   edm::Ptr<flashgg::DiPhotonCandidate> diphoCand;

   //Begin DiPhoton Loop/Selection -----------------------------------------------------------
   for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ )
   {
	 edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
	 
	 dipho_pt = dipho->pt(); dipho_eta = dipho->eta(); dipho_phi = dipho->phi(); dipho_mass = dipho->mass();
		 
	 if(dipho_mass < diph_mass[0] || dipho_mass > diph_mass[1]) continue;
	 if(fabs(dipho_eta) > diph_eta[0] ) continue;
	 if(dipho_pt < diph_pt[0] ) continue;
		 
	 pho1_pt = dipho->leadingPhoton()->pt();			pho2_pt = dipho->subLeadingPhoton()->pt();
	 pho1_eta = dipho->leadingPhoton()->superCluster()->eta();	pho2_eta = dipho->subLeadingPhoton()->superCluster()->eta();
	 pho1_phi = dipho->leadingPhoton()->superCluster()->phi();	pho2_phi = dipho->subLeadingPhoton()->superCluster()->phi();
	 pho1_hoe = dipho->leadingPhoton()->hadronicOverEm(); 		pho2_hoe = dipho->subLeadingPhoton()->hadronicOverEm();
	 pho1_sieie = dipho->leadingPhoton()->full5x5_sigmaIetaIeta();	pho2_sieie = dipho->subLeadingPhoton()->full5x5_sigmaIetaIeta();
	 pho1_r9 = dipho->leadingPhoton()->r9();			pho2_r9 = dipho->subLeadingPhoton()->r9();
	 pho1_elveto = dipho->leadingPhoton()->passElectronVeto();	pho2_elveto = dipho->subLeadingPhoton()->passElectronVeto();

	 pho1_chiso = bbgPlotter::getCHisoToCutValue( dipho, 0, rhoFixedGrd);
	 pho2_chiso = bbgPlotter::getCHisoToCutValue( dipho, 1, rhoFixedGrd);
	 pho1_nhiso = bbgPlotter::getNHisoToCutValue( dipho->leadingPhoton(), rhoFixedGrd );
	 pho2_nhiso = bbgPlotter::getNHisoToCutValue( dipho->subLeadingPhoton(), rhoFixedGrd );
	 pho1_phiso = bbgPlotter::getPHisoToCutValue( dipho->leadingPhoton(), rhoFixedGrd );
	 pho2_phiso = bbgPlotter::getPHisoToCutValue( dipho->subLeadingPhoton(), rhoFixedGrd );

		 
	 if( pho1_pt < dipho_mass*ph_pt[0] ) continue;
	 if( fabs(pho1_eta) > ph_eta[1] ) continue;
	 if( pho2_pt < dipho->mass()*ph_pt[1] ) continue;
	 if( fabs(pho2_eta) > ph_eta[1] ) continue;
		 
	 bool pho1_id = true, pho2_id = true;
	 if( ph_doID[0] )
	 {
	   int pho1Index = 0;
	   if( fabs(pho1_eta) > ph_eta[0] ) pho1Index = 1;
			 
	   if( pho1_hoe > ph_hoe[pho1Index] ) 	pho1_id = false;
	   if( pho1_sieie > ph_sieie[pho1Index] ) pho1_id = false;
	   if( pho1_chiso > ph_chIso[pho1Index] ) pho1_id = false;
	   if( pho1_nhiso > ph_nhIso[pho1Index] ) pho1_id = false;
	   if( pho1_phiso > ph_phIso[pho1Index] ) pho1_id = false;
	   if( pho1_elveto != ph_elVeto[0] ) 	pho1_id = false;
	}
	if( ph_doID[1] )
	{
	   int pho2Index = 2;
	   if( fabs(pho2_eta) > ph_eta[0] ) pho2Index = 3;
		 
	   if( pho2_hoe > ph_hoe[pho1Index] ) 	pho2_id = false;
	   if( pho2_sieie > ph_sieie[pho1Index] ) pho2_id = false;
	   if( pho2_chiso > ph_chIso[pho1Index] ) pho2_id = false;
	   if( pho2_nhiso > ph_nhIso[pho1Index] ) pho2_id = false;
	   if( pho2_phiso > ph_phIso[pho1Index] ) pho2_id = false;
	   if( pho2_elveto != ph_elVeto[1] ) 	pho2_id = false;
	}

	if(pho1_id == true && pho2_id == true){
	  isValidDiPhotonCandidate = true;
  	  CandVtx = dipho->vtx();
	  diphoCand = dipho;
	  break;
	}

   }
   if( isValidDiPhotonCandidate == false ) return;
   //End DiPhoton Loop/Selection -----------------------------------------------------------
   
   //Begin Jets Loop/Selection -----------------------------------------------------------
   vector<edm::Ptr<flashgg::Jet>> bJet;
   vector<edm::Ptr<flashgg::Jet>> lJet;
   int nJet1 = 0, nJet2 = 0;
   for( unsigned int jetIndex = 0; jetIndex < theJets->size(); jetIndex++ )
   {
   	edm::Ptr<flashgg::Jet> jet = theJets->ptrAt( jetIndex );
   	bool isJet1 = true, isJet2 = true;
	
   	if(jet->pt() < jt_pt[0]) isJet1 = false;
   	if(fabs(jet->eta()) > jt_eta[0] ) isJet1 = false;
   	if( jt_doPU[0] && jet->passesPuJetId(CandVtx) == 0 ) isJet1 = false;
	
   	if(jet->pt() < jt_pt[1]) isJet2 = false;
   	if(fabs(jet->eta()) > jt_eta[1] ) isJet2 = false;
   	if( jt_doPU[1] && jet->passesPuJetId(CandVtx) == 0 ) isJet2 = false;
	
   	if(isJet1) nJet1++;
   	if(isJet1 == false && isJet2) nJet2++;
	
   	if( jet->bDiscriminator(bTagType) > jt_bDis[0] ) bJet.push_back(jet);
   	if( jet->bDiscriminator(bTagType) < jt_bDis[0]  && jet->bDiscriminator(bTagType) > jt_bDis[1] ) lJet.push_back(jet);
   }

   int totJets = nJet1 + nJet2;
   if( nJet1 < 1 || totJets < 2 ) return;
   if( nJet.size() < n_bJets ) return;

   edm::Ptr<flashgg::Jet> jet1, jet2;
   TLorentzVector DiJet(0,0,0,0);
   double dijetPt_ref = 0;
   bool hasDiJet = false;

   if(bJet.size() > 1)
   {
   	for(unsigned int jt1 = 0; jt1 < bJet.size(); jt1++){
   		for(unsigned int jt2 = jt1+1; jt2 < bJet.size(); jt2++){
   			TLorentzVector dijet = bJet[jt1]->p4() + bJet[jt2]->p4();
   			if(dijet.pt() > dijetPt_ref && dijet.pt() > dijt_pt[0] && fabs(dijet.Eta()) < dijt_eta[0] ){
				hasDiJet = true;
   				dijetPt_ref = dijet.pt();
   				DiJet = dijet;
   				if( bJet[jt1]->pt() > bJet[jt2]->pt() ) {
   					jet1 = bJet[jt1];
   					jet2 = bJet[jt2];
   				} else {
   					jet2 = bJet[jt1];
   					jet1 = bJet[jt2];
   				}

   			}
		}
   	}
   } else
   {
   	for(unsigned int jt1 = 0; jt1 < bJet.size(); jt1++){
   		for(unsigned int jt2 = 0; jt2 < lJet.size(); jt2++){
   			TLorentzVector dijet = bJet[jt1]->p4() + lJet[jt2]->p4();
   			if(dijet.pt() > dijetPt_ref && dijet.pt() > dijt_pt[0] && fabs(dijet.Eta()) < dijt_eta[0] ){
				hasDiJet = true;
   				dijetPt_ref = dijet.pt();
   				DiJet = dijet;
   				if( bJet[jt1]->pt() > lJet[jt2]->pt() ) {
   					jet1 = bJet[jt1];
   					jet2 = lJet[jt2];
   				} else {
   					jet2 = bJet[jt1];
   					jet1 = lJet[jt2];
   				}
   			}
		}
   	}
   }
   
   if( hasDiJet == false ) return;
   
   dijet_pt = DiJet.Pt();
   dijet_eta = DiJet.Eta();
   dijet_phi = DiJet.Phi();
   dijet_mass = DiJet.M();
	
   jet1_pt = jet1->pt();
   jet1_eta = jet1->eta();
   jet1_phi = jet1->phi();
   jet1_bDis = jet1->bDiscriminator(bTagType);
   jet1_PUid = jet1->passesPuJetId(CandVtx);
	
   jet2_pt = jet2->pt();
   jet2_eta = jet2->eta();
   jet2_phi = jet2->phi();
   jet2_bDis = jet2->bDiscriminator(bTagType);
   jet2_PUid = jet2->passesPuJetId(CandVtx);
   //End Jets Loop/Selection -----------------------------------------------------------
   
   //Candidate assignment --------------------------------------------------------------
   TLorentzVector HHCandidate = DiJet + diphoCand->p4();
   cand4_pt = HHCandidate.Pt();
   cand4_eta = HHCandidate.Eta();
   cand4_phi = HHCandidate.Phi();
   cand4_mass = HHCandidate.M();
   if(cand4_pt < cand_pt[0] ) return;
   if(fabs(cand4_eta) > cand_eta[0] ) return;
   if(cand4_mass < cand_mass[0] || cand4_mass > cand_mass[1] ) return;
   
   //END Candidate assignment ---------------------------------------------------------- 
	 
   hists["dipho_pt"]->Fill(dipho_pt);
   hists["dipho_eta"]->Fill(dipho_eta);
   hists["dipho_phi"]->Fill(dipho_phi);
   hists["dipho_mass"]->Fill(dipho_mass);
	
   hists["dijet_pt"]->Fill(dijet_pt);
   hists["dijet_eta"]->Fill(dijet_eta);
   hists["dijet_phi"]->Fill(dijet_phi);
   hists["dijet_mass"]->Fill(dijet_mass);
	
   hists["cand4_pt"]->Fill(cand4_pt);
   hists["cand4_eta"]->Fill(cand4_eta);
   hists["cand4_phi"]->Fill(cand4_phi);
   hists["cand4_mass"]->Fill(cand4_mass);
	
   hists["pho1_pt"]->Fill(pho1_pt);
   hists["pho1_eta"]->Fill(pho1_eta);
   hists["pho1_phi"]->Fill(pho1_phi);
	
   hists["pho1_hoe"]->Fill(pho1_hoe);
   hists["pho1_sieie"]->Fill(pho1_sieie);
   hists["pho1_r9"]->Fill(pho1_r9);
   hists["pho1_chiso"]->Fill(pho1_chiso);
   hists["pho1_nhiso"]->Fill(pho1_nhiso);
   hists["pho1_phiso"]->Fill(pho1_phiso);
   hists["pho1_elveto"]->Fill(pho1_elveto);
	
   hists["pho2_pt"]->Fill(pho2_pt);
   hists["pho2_eta"]->Fill(pho2_eta);
   hists["pho2_phi"]->Fill(pho2_phi);
	
   hists["pho2_hoe"]->Fill(pho2_hoe);
   hists["pho2_sieie"]->Fill(pho2_sieie);
   hists["pho2_r9"]->Fill(pho2_r9);
   hists["pho2_chiso"]->Fill(pho2_chiso);
   hists["pho2_nhiso"]->Fill(pho2_nhiso);
   hists["pho2_phiso"]->Fill(pho2_phiso);
   hists["pho2_elveto"]->Fill(pho2_elveto);
	
   hists["jet1_pt"]->Fill(jet1_pt);
   hists["jet1_eta"]->Fill(jet1_eta);
   hists["jet1_phi"]->Fill(jet1_phi);
   hists["jet1_bDis"]->Fill(jet1_bDis);
   hists["jet1_PUid"]->Fill(jet1_PUid);
	
   hists["jet2_pt"]->Fill(jet2_pt);
   hists["jet2_eta"]->Fill(jet2_eta);
   hists["jet2_phi"]->Fill(jet2_phi);
   hists["jet2_bDis"]->Fill(jet2_bDis);
   hists["jet2_PUid"]->Fill(jet2_PUid);

}

// ------------ Extra methods... -------------
//
double bbgPlotter::getCHisoToCutValue(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho, float rho)
{
	double PFIso, eta;
	if(whichPho == 0) {
		PFIso = dipho->leadingView().pfChIso03WrtChosenVtx();
		eta = dipho->leadingPhoton()->superCluster()->eta();
	}
	if(whichPho == 1) {
		PFIso = dipho->subLeadingView().pfChIso03WrtChosenVtx();
		eta = dipho->subLeadingPhoton()->superCluster()->eta();
	}
	
	double EA = bbgPlotter::getEA(eta, 0);
	double finalValue = max(PFIso - rho*EA, 0.);
	return finalValue;
}

double bbgPlotter::getNHisoToCutValue(flashgg::Photon* pho, float rho)
{
	double PFIso = pho->pfNeutIso03();
	double eta = pho->superCluster()->eta();
	double EA = bbgPlotter::getEA(eta, 1);
	double extraFactor = 0;
	if(fabs(eta) < 1.479) extraFactor = exp(0.0028*pho->pt()+0.5408);
	if(fabs(eta) > 1.479) extraFactor = 0.01725*pho->pt();
	double finalValue = max(PFIso - rho*EA, 0.) - extraFactor;
	return finalValue;
}

double bbgPlotter::getPHisoToCutValue(flashgg::Photon* pho, float rho)
{
	double PFIso = pho->pfPhoIso03();
	double eta = pho->superCluster()->eta();
	double EA = bbgPlotter::getEA(eta, 2);
	double extraFactor = 0;
	if(fabs(eta) < 1.479) extraFactor = 0.0014*pho->pt();
	if(fabs(eta) > 1.479) extraFactor = 0.0091*pho->pt();
	double finalValue = max(PFIso - rho*EA, 0.) - extraFactor;
	return finalValue;
}

double bbgPlotter::getEA( float eta, int whichEA){
        if(whichEA < 0 || whichEA > 2){
                cout << "WRONG EA TYPE" << endl;
                return -1;
        }

        float EA[7][3];

        EA[0][0] = 0.0234; EA[0][1] = 0.0053; EA[0][2] = 0.078;
        EA[1][0] = 0.0189; EA[1][1] = 0.0103; EA[1][2] = 0.0629;
        EA[2][0] = 0.0171; EA[2][1] = 0.0057; EA[2][2] = 0.0264;
        EA[3][0] = 0.0129; EA[3][1] = 0.0070; EA[3][2] = 0.0462;
        EA[4][0] = 0.0110; EA[4][1] = 0.0152; EA[4][2] = 0.0740;
        EA[5][0] = 0.0074; EA[5][1] = 0.0232; EA[5][2] = 0.0924;
        EA[6][0] = 0.0035; EA[6][1] = 0.1709; EA[6][2] = 0.1484;

        float feta = fabs(eta);

        if(feta > 0.000 && feta < 1.000 ) return EA[0][whichEA];
        if(feta > 1.000 && feta < 1.479 ) return EA[1][whichEA];
        if(feta > 1.479 && feta < 2.000 ) return EA[2][whichEA];
        if(feta > 2.000 && feta < 2.200 ) return EA[3][whichEA];
        if(feta > 2.200 && feta < 2.300 ) return EA[4][whichEA];
        if(feta > 2.300 && feta < 2.400 ) return EA[5][whichEA];
        if(feta > 2.400 && feta < 10.00 ) return EA[6][whichEA];

        return -1;
}


// ------------ method called once each job just before starting event loop  ------------
void 
bbggPlotter::beginJob()
{
	outFile = new TFile("out.root", "RECREATE");

	hists["dipho_pt"] 	= new TH1F("dipho_pt", "DiPhoton p_{T}; p_{T}(#gamma#gamma) (GeV); Events", 100, 0, 300);
	hists["dipho_eta"] 	= new TH1F("dipho_eta", "DiPhoton #eta; #eta(#gamma#gamma); Events", 100, -5, 5);
	hists["dipho_phi"]	= new TH1F("dipho_phi", "DiPhoton #phi; #phi(#gamma#gamma); Events", 100, -6.5, 6.5);
	hists["dipho_mass"] 	= new TH1F("dipho_mass", "DiPhoton Mass; M(#gamma#gamma); Events", 100, diph_mass[0], diph_mass[1]);

	hists["dijet_pt"] 	= new TH1F("dijet_pt", "DiJet p_{T}; p_{T}(jj) (GeV); Events", 100, 0, 300);
	hists["dijet_eta"] 	= new TH1F("dijet_eta", "DiJet #eta; #eta(jj); Events", 100, -5, 5);
	hists["dijet_phi"]	= new TH1F("dijet_phi", "DiJet #phi; #phi(jj); Events", 100, -6.5, 6.5);
	hists["dijet_mass"] 	= new TH1F("dijet_mass", "DiJet Mass; M(jj); Events", 100, dijt_mass[0], dijt_mass[1]);

	hists["cand_pt"] 	= new TH1F("cand_pt", "DiHiggs Candidate (jj#gamma#gamma) p_{T}; p_{T}(jj#gamma#gamma) (GeV); Events", 100, 0, 500);
	hists["cand_eta"] 	= new TH1F("cand_eta", "DiHiggs Candidate (jj#gamma#gamma) #eta; #eta(jj#gamma#gamma); Events", 100, -5, 5);
	hists["cand_phi"]	= new TH1F("cand_phi", "DiHiggs Candidate (jj#gamma#gamma) #phi; #phi(jj#gamma#gamma); Events", 100, -6.5, 6.5); 
	hists["cand_mass"] 	= new TH1F("cand_mass", "DiHiggs Candidate (jj#gamma#gamma) Mass; M(jj#gamma#gamma) (GeV); Events", 100, cand_mass[0], cand_mass[1]);

	hists["pho1_pt"] 	= new TH1F("pho1_pt", "Leading Photon p_{T}; p_{T}(leading #gamma) (GeV); Events", 100, 10, 150);
	hists["pho1_eta"] 	= new TH1F("pho1_eta", "Leading Photon #eta; #eta(leading #gamma); Events", 100, -5., 5.);
	hists["pho1_phi"]	= new TH1F("pho1_phi", "Leading Photon #phi; #phi(leading #gamma); Events", 100, -6.5, 6.5);

	hists["pho1_hoe"] 	= new TH1F("pho1_hoe", "Leading Photon H/E; H/E(leading #gamma); Events", 100, -0.01, 0.1);
	hists["pho1_sieie"] 	= new TH1F("pho1_sieie", "Leading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(leading #gamma); Events", 100, 0.0, 0.04);
	hists["pho1_r9"] 	= new TH1F("pho1_r9", "Leading Photon R9; R9(leading #gamma); Events", 100, 0, 1.1);
	hists["pho1_chiso"] 	= new TH1F("pho1_chiso", "Leading Photon Corrected Charged Isolation; Corrected Charged Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_nhiso"] 	= new TH1F("pho1_nhiso", "Leading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_phiso"] 	= new TH1F("pho1_phiso", "Leading Photon Corrected Photon Isolation; Corrected Photon Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_elveto"] 	= new TH1F("pho1_elveto", "Leading Photon Electron Veto; Electron Veto(leading #gamma); events", 8, -1, 3);

	hists["pho2_pt"] 	= new TH1F("pho2_pt", "SubLeading Photon p_{T}; p_{T}(subLeading #gamma) (GeV); Events", 100, 10, 150);
	hists["pho2_eta"] 	= new TH1F("pho2_eta", "SubLeading Photon #eta; #eta(subLeading #gamma); Events", 100, -5., 5.);
	hists["pho2_phi"]	= new TH1F("pho2_phi", "SubLeading Photon #phi; #phi(subleading #gamma); Events", 100, -6.5, 6.5);

	hists["pho2_hoe"] 	= new TH1F("pho2_hoe", "SubLeading Photon H/E; H/E(subLeading #gamma); Events", 100, -0.01, 0.1);
	hists["pho2_sieie"] 	= new TH1F("pho2_sieie", "SubLeading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(subLeading #gamma); Events", 100, 0.0, 0.04);
	hists["pho2_r9"] 	= new TH1F("pho2_r9", "SubLeading Photon R9; R9(subLeading #gamma); Events", 100, 0, 1.1);
	hists["pho2_chiso"] 	= new TH1F("pho2_chiso", "SubLeading Photon Corrected Charged Isolation; Corrected Charged Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_nhiso"] 	= new TH1F("pho2_nhiso", "SubLeading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_phiso"] 	= new TH1F("pho2_phiso", "SubLeading Photon Corrected Photon Isolation; Corrected Photon Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_elveto"] 	= new TH1F("pho2_elveto", "SubLeading Photon Electron Veto; Electron Veto(subLeading #gamma); events", 8, -1, 3);


	hists["jet1_pt"] 	= new TH1F("jet1_pt", "Leading Jet p_{T}; p_{T}(leading jet) (GeV); Events", 100, 10, 150);
	hists["jet1_eta"] 	= new TH1F("jet1_eta", "Leading Jet #eta; #eta(leading jet); Events", 100, -5., 5.);
	hists["jet1_phi"]	= new TH1F("jet1_phi", "Leading Jet #phi; #phi(leading jet); Events", 100, -6.5, 6.5);
	hists["jet1_bDis"] 	= new TH1F("jet1_bDis", "Leading Jet b-Discriminant; b-Discriminant(leading jet); Events", 100, -0.01, 1.01);
	hists["jet1_PUid"] 	= new TH1F("jet1_PUid", "Leading Jet PU ID; PU ID(leading jet); Events", 8, -1, 3);

	hists["jet2_pt"] 	= new TH1F("jet2_pt", "SubLeading Jet p_{T}; p_{T}(subLeading jet) (GeV); Events", 100, 10, 150);
	hists["jet2_eta"] 	= new TH1F("jet2_eta", "SubLeading Jet #eta; #eta(subLeading jet); Events", 100, -5., 5.);
	hists["jet2_phi"]	= new TH1F("jet2_phi", "SubLeading Jet #phi; #phi(subleading jet); Events", 100, -6.5, 6.5);
	hists["jet2_bDis"] 	= new TH1F("jet2_bDis", "SubLeading Jet b-Discriminant; b-Discriminant(subLeading jet); Events", 100, -0.01, 1.01);
	hists["jet2_PUid"] 	= new TH1F("jet2_PUid", "SubLeading Jet PU ID; PU ID(subleading jet); Events", 8, -1, 3);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbggPlotter::endJob() 
{
	outFile->cd();
	for(std::map<std::string, TH1F*>::iterator it = hists.begin(); it != hists.end(); ++it)
	{
		cout << "Saving histogram... " << it->first << endl;
		it->second->Write();
	}
	outFile->Close();
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
