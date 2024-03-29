//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1EventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "B1Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction *runAction, G4bool LightOutFlag, G4bool OnlyEnterSiPM)
	: G4UserEventAction(),
	  fRunAction(runAction),
	  fEdep(0.),
	  fEdepSiPM(0.),
	  fEdkin(0.),
	  fNumHitsDet(0),
	  fPrePterNo(0),
	  fPreProbeNo(0),
	  fPostAbsNo(0),
	  fEdepEle(0.),
	  fEdepPos(0.),
	  fEdepFot(0.),
	  fEdepSiPMpos(0.),
	  fEdepSiPMfot(0.),
	  fEnteringParticle(0),
	  fSourceExitPassCounter(0.),
	  fPterPassCounter(0.),
	  fPostAbsPassCounter(0.),
	  fPreProbePassCounter(0.),
	  fNSourceExit(0.),
	  fSourceExitStoreTrackID(0),
	  fPterStoreTrackID(0),
	  fPostAbsStoreTrackID(0),
	  fPreProbeStoreTrackID(0),
	  fEnterPterFlag(0),
	  fLightOutFlag(LightOutFlag),
	  fOnlyEnterSiPM(OnlyEnterSiPM)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event *)
{
	// G4cout << "SANREMO " << fOnlyEnterSiPM << G4endl;
	fEdep = 0.;
	fEdepSiPM = 0.;
	fEdkin = 0.;
	(fRunAction->GetPrePterEn()).clear();
	(fRunAction->GetPrePterPart()).clear();
	(fRunAction->GetPrePterX()).clear();
	(fRunAction->GetPrePterY()).clear();
	(fRunAction->GetPrePterZ()).clear();
	(fRunAction->GetPreProbeEn()).clear();
	(fRunAction->GetPreProbePart()).clear();
	(fRunAction->GetPostAbsEn()).clear();
	(fRunAction->GetPostAbsPart()).clear();
	(fRunAction->GetPterEn()).clear();
	(fRunAction->GetPterEnPrim()).clear();
	(fRunAction->GetPterPartPrim()).clear();
	(fRunAction->GetPterTime()).clear();
	(fRunAction->GetPterX()).clear();
	(fRunAction->GetPterY()).clear();
	(fRunAction->GetPterZ()).clear();
	(fRunAction->GetPterPart()).clear();

	(fRunAction->GetAnnihX()).clear();
	(fRunAction->GetAnnihY()).clear();
	(fRunAction->GetAnnihZ()).clear();

	(fRunAction->GetAnnihT()).clear();

	(fRunAction->GetPostAbsPart()).clear();

	(fRunAction->GetSourceCosX()).clear();
	(fRunAction->GetSourceCosY()).clear();
	(fRunAction->GetSourceCosZ()).clear();

	(fRunAction->GetSourceEn()).clear();
	(fRunAction->GetSourcePart()).clear();
	(fRunAction->GetSourceIsotope()).clear();

	(fRunAction->SetMotherIsotope(-10));
	(fRunAction->SetMotherEnergy(-10));
	(fRunAction->SetMotherPID(999));
	(fRunAction->SetMotherTime(0));

	(fRunAction->GetExitEn()).clear();
	(fRunAction->GetExitX()).clear();
	(fRunAction->GetExitY()).clear();
	(fRunAction->GetExitZ()).clear();
	(fRunAction->GetExitCosX()).clear();
	(fRunAction->GetExitCosY()).clear();
	(fRunAction->GetExitCosZ()).clear();
	(fRunAction->GetExitPart()).clear();
	(fRunAction->GetExitParentID()).clear();
	(fRunAction->GetExitOrigX()).clear();
	(fRunAction->GetExitOrigY()).clear();
	(fRunAction->GetExitOrigZ()).clear();

	(fRunAction->GetExitProcess()).clear();

	(fRunAction->GetEAbsComp()).clear();
	(fRunAction->GetEAbsSiPMComp()).clear();
	(fRunAction->GetSiPMPart()).clear();
	(fRunAction->GetSiPMEne()).clear();
	(fRunAction->GetSiPMVX()).clear();
	(fRunAction->GetSiPMVY()).clear();
	(fRunAction->GetSiPMVZ()).clear();

	fNumHitsDet = 0;
	fPrePterNo = 0;
	fPreProbeNo = 0;
	fPostAbsNo = 0;

	fEdepEle = 0.;
	fEdepPos = 0;
	fEdepFot = 0.;
	fEdepSiPMpos = 0.,
	fEdepSiPMfot = 0.,
	fEnteringParticle = 0;
	fNSourceExit = 0;
	fSourceExitPassCounter = 0;
	fPterPassCounter = 0;
	fPostAbsPassCounter = 0;
	fPreProbePassCounter = 0;
	fSourceExitStoreTrackID = 0;
	fPterStoreTrackID = 0;
	fPostAbsStoreTrackID = 0;
	fPreProbeStoreTrackID = 0;
	fNPMT = 0;
	fEnterPterFlag = 0;
	fSourceReg = "";

	//	fSourceX=0;
	//	fSourceY=0;
	//	fSourceZ=0;
	/*
		 fSourceEne=0;
	 */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event *evento)
{

	fRunAction->AddEdep(fEdep);
	//	fRunAction->AddEdepSiPM(fEdepSiPM);
	fRunAction->AddEdkin(fEdkin);

	(fRunAction->GetEAbsComp()).push_back(fEdepEle / keV);
	(fRunAction->GetEAbsComp()).push_back(fEdepPos / keV);
	(fRunAction->GetEAbsComp()).push_back(fEdepFot / keV);
	//	(fRunAction->GetEAbsSiPMComp()).push_back(fEdepSiPMpos/keV);
	//	(fRunAction->GetEAbsSiPMComp()).push_back(fEdepSiPMfot/keV);

	G4int NevTot = fRunAction->GetEventNumber();

	if ((10 * evento->GetEventID()) % NevTot == 0)
	{
		G4cout << "Progress status: " << (evento->GetEventID() / (G4double)NevTot) * 100 << " %, Nev= " << evento->GetEventID() << ", NTotEv= " << NevTot << G4endl;
	}
	// get analysis manager
	auto analysisManager = G4RootAnalysisManager::Instance();

	// fill ntuple
	analysisManager->FillNtupleDColumn(0, 0, fEdep / keV);
	// analysisManager->FillNtupleDColumn(0, 34, fEdepSiPM/keV);
	analysisManager->FillNtupleDColumn(0, 2, fNumHitsDet); // number of hits into the detector
	analysisManager->FillNtupleDColumn(0, 11, fPrePterNo);
	analysisManager->FillNtupleDColumn(0, 17, fPreProbeNo);
	analysisManager->FillNtupleDColumn(0, 20, fPostAbsNo);

	analysisManager->FillNtupleDColumn(0, 25, fSourceX / mm);
	analysisManager->FillNtupleDColumn(0, 26, fSourceY / mm);
	analysisManager->FillNtupleDColumn(0, 27, fSourceZ / mm);

	analysisManager->FillNtupleIColumn(0, 34, fNPMT);
	analysisManager->FillNtupleIColumn(0, 35, fEnterPterFlag);

	// if(!fLightOutFlag||fEdep>0) analysisManager->AddNtupleRow(0);
	if (
		(!fLightOutFlag && !fOnlyEnterSiPM) // NON ho richiesto la flag light
		|| fEdep > 0						// L'energia depositata è >0
											// || (!fLightOutFlag && !fOnlyEnterSiPM) // NON ho chiesto di salvare solo quello che entra nel SiPm
	)
	{
		// G4cout << "SANREMO ENTRO" << G4endl;
		analysisManager->AddNtupleRow(0);
	}
	else
	{
		// G4cout << "SANREMO NON ENTRO" << G4endl;
	}

	if (evento->GetEventID() <= 1e5)
	{ // to write to proper ntuple all the source particles info
		analysisManager->FillNtupleDColumn(1, 0, fSourceX / mm);
		analysisManager->FillNtupleDColumn(1, 1, fSourceY / mm);
		analysisManager->FillNtupleDColumn(1, 2, fSourceZ / mm);
		analysisManager->FillNtupleSColumn(1, 9, fSourceReg);
		analysisManager->FillNtupleDColumn(1, 23, fNSourceExit);

		analysisManager->AddNtupleRow(1);
	}

	// Keep interesting events for offline visualization
	////	if (fPterPassCounter==0 && fPreProbePassCounter>0) {
	//			if (fPrePterNo>0 && fPreProbeNo==0) {
	////		if (fPterPassCounter>fPreProbePassCounter) {
	//		G4cout<<"Evento con tracce in pter ma 0 in Probe: "<<evento->GetEventID()<<" "<<fPrePterNo<<" "<< fPreProbeNo  <<G4endl;
	//		G4Event* evt = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
	//		evt->KeepTheEvent();
	//
	//	}
	//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
