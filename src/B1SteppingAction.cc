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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1RunAction.hh"
#include "B1DetectorConstruction.hh"


#include "G4Step.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

#include "B1Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction, B1RunAction* runAction, G4double AbsHoleDiam, G4int GaSet)
: G4UserSteppingAction(),
fEventAction(eventAction),
fScoringVolume(0),
runStepAction(runAction),
fAbsHoleDiam(AbsHoleDiam),
fGaSet(GaSet)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
	
	G4VPhysicalVolume* ThisVol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4VPhysicalVolume* NextVol = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	
	G4int debug=0;
	
	//In case of GaSet 2/3 look for annihilation points (since it's probably Gallium)
	if ((fGaSet==2 ||fGaSet==3) && step->GetTrack()->GetDynamicParticle() ->GetPDGcode() == -11 && step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="annihil") {
		G4Event* evt = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
		evt->KeepTheEvent();
		(runStepAction->GetRunAnnihX()).push_back(step->GetPostStepPoint()->GetPosition().x()/mm);
		(runStepAction->GetRunAnnihY()).push_back(step->GetPostStepPoint()->GetPosition().y()/mm);
		(runStepAction->GetRunAnnihZ()).push_back(step->GetPostStepPoint()->GetPosition().z()/mm);
		(runStepAction->GetRunAnnihT()).push_back(step->GetPostStepPoint()->GetLocalTime()/ns);
	}
	
	// ########################################
	// ###################### Optical Photons ENTERING SiPm
	if(step->GetTrack()->GetDynamicParticle() ->GetPDGcode()== 0 && NextVol && ThisVol->GetName()=="Pter" && NextVol->GetName()=="SiPm") {
//		G4cout<<"FOTONE OTTICO ENTRA IN SiPm"<<G4endl;
		fEventAction->AddNPMT(1);
	}
	// ######################
	// ########################################
	
	
	// ########################################
	// ###################### Interactions in SiPm
	if (ThisVol->GetName()=="SiPm") {
		fEventAction->AddEdepSiPM(step->GetTotalEnergyDeposit());
//		G4cout<<"INTERAZIONE NEL SIPM: DepEne= "<<step->GetTotalEnergyDeposit()/keV<<" Part= "<<step->GetTrack()->GetDynamicParticle() ->GetPDGcode()<<G4endl;
	}
	// ######################
	// ########################################
	
	
	// ########################################
	// ###################### ENTERING Pter (from wherever)
	
	//	if((NextVol && ThisVol->GetName()=="FrontShield" && NextVol->GetName()=="Pter")|| (NextVol && ThisVol->GetName()=="World" && NextVol->GetName()=="Pter")) { //what enters Pter (either from FrontShield or world)
	if((NextVol && ThisVol->GetName()!="Pter" && NextVol->GetName()=="Pter")) { //what enters Pter (form every different volume)
		
		if (debug) G4cout<<"\nCIAODEBUG\n Particella entrata in PTER - fEventAction->GetEnteringParticle() ERA = "<<fEventAction->GetEnteringParticle();
		fEventAction->SetEnteringParticle(step->GetTrack()->GetDynamicParticle()->GetPDGcode());
		if (debug) G4cout<<" SETTO fEventAction->GetEnteringParticle()= "<<fEventAction->GetEnteringParticle()<<G4endl<<G4endl;
		
		// Check if the current track had already enetered somehow the PTER (to avoid double counting), otherwise increase the counter
		if (fEventAction->GetStoreTrackIDPter()==step->GetTrack()->GetTrackID()) { //if I already saw this track exiting the source...
			fEventAction->AddPassCounterPter(1);  //increase the counter
			//			G4cout<<"PterDEBUG CONTROLLA "<<fEventAction->GetStoreTrackIDPter()<<", PassCounter= "<<fEventAction->GetPassCounterPter()<<G4endl;
		}else {
			fEventAction->SetStoreTrackIDPter(step->GetTrack()->GetTrackID());
			//			G4cout<<"PterDEBUG PRIMO PASSAGGIO!! "<<fEventAction->GetStoreTrackIDPter()<<", PassCounter= "<<fEventAction->GetPassCounterPter()<<G4endl;
			//            if (fEventAction->GetPassCounter()!=0) G4cout<<"MERDAAAAA Primo passaggio di"<<fEventAction->GetStoreTrackID()<<" ma con PassCounter= "<<fEventAction->GetPassCounter()<<G4endl;
		}
		
		// Salvo le info solo della prima volta che una particella entra in pter
		if (fEventAction->GetPassCounterPter()==0) {
			G4double eKinPre = step->GetPostStepPoint()->GetKineticEnergy();
			//Fill vector
			(runStepAction->GetRunEnPre()).push_back(eKinPre/keV);
			fEventAction->AddNoPre(1); //update the counter of particles entering Pter in the event
			(runStepAction->GetRunPart()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode()); //add PID of particle enetering Pter
		}
	}
	
	// ###################### END ENTERING Pter
	// ########################################
	
	
	//Modified on 2017-11-17 by collamaf: now the condition works for both cases: with or without Cu collimator.
	//If there is not collimator save what goes from source to dummy. If there is a collimator save what goes from world (the hole) into dummy
	
	
	// ########################################
	// ###################### EXITING SOURCE
	if( NextVol && ThisVol->GetName()!="DummyExitSorg" && NextVol->GetName()=="DummyExitSorg" && step->GetPreStepPoint()->GetMomentumDirection().z()>0)
	{
			 
			 //collamaf: to avoid double counting same track going back and forth, check if I already counted it
			 if (fEventAction->GetSourceExitStoreTrackID()==step->GetTrack()->GetTrackID()) { //if I already saw this track exiting the source...
				 fEventAction->AddSourceExitPassCounter(1);  //increase the counter
			 }else {
				 fEventAction->SetSourceExitStoreTrackID(step->GetTrack()->GetTrackID());
			 }
			 
			 // Salvo le info solo della prima volta che una particella esce dalla sorgente
			 if (fEventAction->GetSourceExitPassCounter()==0) {
				 fEventAction->AddNSourceExit(1);
			 (runStepAction->GetRunEnExit()).push_back(step->GetPostStepPoint()->GetKineticEnergy()/keV);
			 (runStepAction->GetRunXExit()).push_back(step->GetPostStepPoint()->GetPosition().x()/mm);
			 (runStepAction->GetRunYExit()).push_back(step->GetPostStepPoint()->GetPosition().y()/mm);
			 (runStepAction->GetRunZExit()).push_back(step->GetPostStepPoint()->GetPosition().z()/mm);
			 (runStepAction->GetRunCosXExit()).push_back(step->GetPreStepPoint()->GetMomentumDirection().x());
			 (runStepAction->GetRunCosYExit()).push_back(step->GetPreStepPoint()->GetMomentumDirection().y());
			 (runStepAction->GetRunCosZExit()).push_back(step->GetPreStepPoint()->GetMomentumDirection().z());
			 (runStepAction->GetRunPartExit()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
			 (runStepAction->GetRunParentIDExit()).push_back(step->GetTrack()->GetParentID());
			 
			 if (step->GetTrack()->GetCreatorProcess()) {
				 (runStepAction->GetRunExitProcess().push_back((step->GetTrack()->GetCreatorProcess()->GetProcessType())));
			 } else {
				 (runStepAction->GetRunExitProcess().push_back(-17));
			 }
			 (runStepAction->GetRunPartPostAbs()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		 }
		 
		 /*
			We have to use PreStepPoint to save the exit cosines, otherwise we already have particles flipped..
			*/
	 }
	
	
	if (NextVol && ((fAbsHoleDiam>=0 && fGaSet == 2 &&  (ThisVol->GetName()=="SourceExtGa" && NextVol->GetName()=="Absorber") ) )) { //what actually exits the source
		
		//collamaf: to avoid double counting same track going back and forth, check if I already counted it
		if (fEventAction->GetStoreTrackIDDummy2()==step->GetTrack()->GetTrackID()) { //if I already saw this track exiting the source...
			fEventAction->AddPassCounterDummy2(1);  //increase the counter
		}else {
			fEventAction->SetStoreTrackIDDummy2(step->GetTrack()->GetTrackID());
		}
		
		// Salvo le info solo della prima volta che una particella esce dalla sorgente
		if (fEventAction->GetPassCounterDummy2()==0) {
			(runStepAction->GetRunPreAbsEn()).push_back(step->GetPostStepPoint()->GetKineticEnergy()/keV);
			(runStepAction->GetRunPartPreAbs()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
			(runStepAction->GetRunPartExit()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		}
	}
	
	
	
	
	
	if (!fScoringVolume) {
		const B1DetectorConstruction* detectorConstruction
		= static_cast<const B1DetectorConstruction*>
		(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		fScoringVolume = detectorConstruction->GetScoringVolume();
	}
	
	if (0 && runStepAction->GetMotherIsotope() != 0 && runStepAction->GetMotherIsotope() !=1) G4cout<<"PterDEBUG PROVA STEPPING  MotherIsotope Val= "<< runStepAction->GetMotherIsotope()
		<<G4endl;
	
	// get volume of the current step
	G4LogicalVolume* volume
	= step->GetPreStepPoint()->GetTouchableHandle()
	->GetVolume()->GetLogicalVolume();
	
	
	// ########################################
	// ###################### INSIDE Pter - Per each hit into sensitive detector
	// check if we are in scoring volume
	if (volume == fScoringVolume && step->GetTrack()->GetDynamicParticle() ->GetPDGcode()!=0 ) {
		fEventAction->SetEnterPterFlag();
		//pixel information collection
		G4int CopyNB=step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
		fEventAction->AddNo(1);
//		G4cout<<"Buccato drendo"<<G4endl;
		G4ThreeVector pixCenter;
		G4TouchableHandle touchHandle =step->GetPreStepPoint()->GetTouchableHandle();
		G4ThreeVector vec_origin(0.,0.,0.);
		G4ThreeVector globalPos = touchHandle->GetHistory()-> GetTopTransform().Inverse().TransformPoint(vec_origin);
		pixCenter = globalPos;
		
		if (CopyNB>0) {
			//fill vectors
			(runStepAction->GetRunPixNo()).push_back(CopyNB);
			//			(runStepAction->GetRunPixEneDep()).push_back(step->GetTotalEnergyDeposit()/keV);
			(runStepAction->GetRunPixXpos()).push_back(pixCenter.getX()/mm);
			(runStepAction->GetRunPixYpos()).push_back(pixCenter.getY()/mm);
		}
		
		// collect energy deposited in this step
		G4StepPoint* postPoint = step->GetPostStepPoint();
		G4double edepStep = step->GetTotalEnergyDeposit();
		G4ThreeVector post=postPoint->GetPosition();
		
		//Fill vector
		(runStepAction->GetRunEnPter()).push_back(step->GetTotalEnergyDeposit()/keV);
		(runStepAction->GetRunEnPterPrim()).push_back(runStepAction->GetMotherEnergy());
		(runStepAction->GetRunPartPterPrim()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		//		(runStepAction->GetRunEnPterTime()).push_back(step->GetTrack()->GetLocalTime()/ns);
		(runStepAction->GetRunEnPterTime()).push_back(step->GetTrack()->GetGlobalTime()/ns-runStepAction->GetMotherTime());
		//		G4cout<<"PterDEBUG  MotherTime= "<< runStepAction->GetMotherTime()<<" PostDiff= "<<  step->GetTrack()->GetGlobalTime()/ns-runStepAction->GetMotherTime() <<G4endl;
		(runStepAction->GetRunXPter()).push_back(step->GetPreStepPoint()->GetPosition().x()/mm);
		(runStepAction->GetRunYPter()).push_back(step->GetPreStepPoint()->GetPosition().y()/mm);
		(runStepAction->GetRunZPter()).push_back(step->GetPreStepPoint()->GetPosition().z()/mm);
		(runStepAction->GetRunPartPter()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		
		if (debug)  G4cout<<"CIAODEBUG Ho un rilascio di energia ("<< step->GetTotalEnergyDeposit()/keV<<" [KeV]) dovuto ad una particella entrata nel CMOS di tipo: "<<fEventAction->GetEnteringParticle()<<G4endl;
		
		if (fEventAction->GetEnteringParticle()==11) {  //if son of electron
			fEventAction->AddEdepEle(step->GetTotalEnergyDeposit());
		}
		else if (fEventAction->GetEnteringParticle()==-11) {  //if son of positron
			fEventAction->AddEdepPos(step->GetTotalEnergyDeposit());
		} else if (fEventAction->GetEnteringParticle()==22) {  //if son of photon
			fEventAction->AddEdepFot(step->GetTotalEnergyDeposit());
			if (debug&&step->GetTotalEnergyDeposit()>0) G4cout<<"CONTROLLA"<<G4endl;
		}
		
		fEventAction->AddEdep(edepStep);
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
 	if( NextVol && ( (fAbsHoleDiam<0 &&  ( (ThisVol->GetName()=="SourceSR" && NextVol->GetName()=="Dummy") || (ThisVol->GetName()=="SourceExtY" && NextVol->GetName()=="Dummy") || (ThisVol->GetName()=="SourceExtGa" && NextVol->GetName()=="Dummy"))) || ( (fAbsHoleDiam>=0 &&   (ThisVol->GetName()=="World" && NextVol->GetName()=="Dummy") ) )) ) { //what actually exits the source
 */

/*
 if (fGaSet==2 && step->GetTrack()->GetDynamicParticle() ->GetPDGcode() == -11 && step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="annihil" && (NextVol->GetName()=="ProbeContainer" )) {
 G4Event* evt = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
 evt->KeepTheEvent();
 (runStepAction->GetRunAnnihX()).push_back(step->GetPostStepPoint()->GetPosition().x()/mm);
 (runStepAction->GetRunAnnihY()).push_back(step->GetPostStepPoint()->GetPosition().y()/mm);
 (runStepAction->GetRunAnnihZ()).push_back(step->GetPostStepPoint()->GetPosition().z()/mm);
 }
 */


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Nuova Condizione

/*
 if( NextVol && ( ((fAbsHoleDiam<0 || fAbsHoleDiam>=0 ) &&  ( (ThisVol->GetName()=="SourceSR" && (NextVol->GetName()=="Dummy" || NextVol->GetName()=="CuCollimator" || NextVol->GetName()=="World")) || (ThisVol->GetName()=="SourceExtY" && (NextVol->GetName()=="Dummy"|| NextVol->GetName()=="ABSaround" || NextVol->GetName()=="ABSbehind" || NextVol->GetName()=="CuCollimator")) || (ThisVol->GetName()=="SourceExtGa" && (NextVol->GetName()=="GaContainer" || NextVol->GetName()=="Dummy3" || NextVol->GetName()=="Dummy2" || NextVol->GetName()=="Absorber"))) ) ) ) { //what actually exits the source
 */

//Vecchia condizione

/*
 if( NextVol && ( (fAbsHoleDiam<0 &&  ( (ThisVol->GetName()=="SourceSR" && NextVol->GetName()=="Dummy") || (ThisVol->GetName()=="SourceExtY" && NextVol->GetName()=="Dummy") || (ThisVol->GetName()=="SourceExtGa" && NextVol->GetName()=="Dummy")  || (ThisVol->GetName()=="SourceExtGa" && NextVol->GetName()=="Dummy2") )) || ( (fAbsHoleDiam>=0 && fGaSet == 2 &&  (ThisVol->GetName()=="Absorber" && NextVol->GetName()=="Dummy2") ) ) || ( (fAbsHoleDiam>=0 && fGaSet == 1 && (ThisVol->GetName()=="CuCollimator" && NextVol->GetName()=="Dummy") ) )   ) )
 */

//Modifica alla vecchia condizione per dummy3

/*
 if( NextVol && ( (fAbsHoleDiam<0 &&  ( (ThisVol->GetName()=="SourceSR" && NextVol->GetName()=="Dummy") || (ThisVol->GetName()=="SourceExtY" && NextVol->GetName()=="Dummy") || (ThisVol->GetName()=="SourceExtGa" && NextVol->GetName()=="Dummy")  || (ThisVol->GetName()=="SourceExtGa" && NextVol->GetName()=="Dummy2") )) || ( (fAbsHoleDiam<0 && fGaSet == 3 &&  (ThisVol->GetName()=="Dummy3" && NextVol->GetName()=="Dummy2") ) ) || ( (fAbsHoleDiam>=0 && fGaSet == 2 &&  (ThisVol->GetName()=="Absorber" && NextVol->GetName()=="Dummy2") ) ) || ( (fAbsHoleDiam>=0 && fGaSet == 1 && (ThisVol->GetName()=="CuCollimator" && NextVol->GetName()=="Dummy") ) ) || ( (fAbsHoleDiam>=0 && fGaSet == 3 &&  (ThisVol->GetName()=="Absorber" && NextVol->GetName()=="Dummy2") ) )  ) )
 */

