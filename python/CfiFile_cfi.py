import FWCore.ParameterSet.Config as cms

#### Comments refer to the structure
#### of the vector below the comment line

bbggAnalyzer = cms.EDAnalyzer('bbggPlotter',
	DiPhotonTag=cms.untracked.InputTag('flashggDiPhotons'),
	JetTag=cms.untracked.InputTag('flashggJets'),
	rhoFixedGridCollection=cms.untracked.InputTag('fixedGridRhoAll'),
	#0: Pho1, 1: Pho2
	PhotonPtOverDiPhotonMass=cms.untracked.vdouble(0., 0.),
	#0: First upper boundary (EB), 1: second upper boundary (EE) (only (1) is used to cut on both)
	PhotonEta=cms.untracked.vdouble(1.479, 2.400),
	#0: Pho1 EB, 1: Pho1 EE, 2: Pho2 EB, 3: Pho2 EE
	PhotonHoverE=cms.untracked.vdouble(0.012, 0.023, 0.012, 0.023 ),
	#0: Pho1 EB, 1: Pho1 EE, 2: Pho2 EB, 3: Pho2 EE
	PhotonSieie=cms.untracked.vdouble(0.0100, 0.0267, 0.0100, 0.0267 ),
	#0: Pho1 EB, 1: Pho1 EE, 2: Pho2 EB, 3: Pho2 EE
	PhotonR9=cms.untracked.vdouble(0., 0.),
	#0: Pho1 EB, 1: Pho1 EE, 2: Pho2 EB, 3: Pho2 EE
	PhotonChargedIso=cms.untracked.vdouble(1.79, 1.09, 1.79, 1.09 ),
	#0: Pho1 EB, 1: Pho1 EE, 2: Pho2 EB, 3: Pho2 EE
	PhotonNeutralIso=cms.untracked.vdouble(0.16, 4.31, 0.16, 4.31),
	#0: Pho1 EB, 1: Pho1 EE, 2: Pho2 EB, 3: Pho2 EE
	PhotonPhotonIso=cms.untracked.vdouble(1.90, 1.90, 1.90, 1.90),
	#0: Pho1, 1: Pho2
	PhotonElectronVeto=cms.untracked.vdouble(0., 0.),
	#0: Pho1, 1: Pho2
	PhotonDoID=cms.untracked.vdouble(0, 0,),
	#0: lower boundary for dipho pt
	DiPhotonPt=cms.untracked.vdouble(0.),
	#0: upper boundary
	DiPhotonEta=cms.untracked.vdouble(2.4),
	#0: DiPhoton mass window lower boundary, 1: upper boundary
	DiPhotonMassWindow=cms.untracked.vdouble(80., 165.),
	#0: jet1, 1: jet2
	JetPtOverDiJetMass=cms.untracked.vdouble(10., 10.),
	#0: jet1, 1: jet2
	JetEta=cms.untracked.vdouble(2.5, 2.5),
	#0: highest b-tag requirement, 1: lowest b-tag requirement
	JetBDiscriminant=cms.untracked.vdouble(0.5, 0.5),
	#0: jet1, 1: jet2
	JetDoPUID=cms.untracked.vdouble(1, 1),
	#Number of required jets passing requirements in JetBDiscriminant
	n_bJets=cms.untracked.int32(0),
	#0: lower boundary for dijet pt
	DiJetPt=cms.untracked.vdouble(10.),
	#0: upper boundary for dijet pt
	DiJetEta=cms.untracked.vdouble(20.),
	#0: DiJet mass window lower boundary, 1: upper boundary
	DiJetMassWindow=cms.untracked.vdouble(40., 200.),
	#0: 4-candidate mass window lower boundary, 1: upper boundary
	CandidateMassWindow=cms.untracked.vdouble(0, 1000.),
	#string for btag algorithm
	bTagType=cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags')
)
