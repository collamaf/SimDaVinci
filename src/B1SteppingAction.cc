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

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction, B1RunAction* RunningAction, G4double AbsHoleDiam, G4int GaSet, G4int SourceChoice)
: G4UserSteppingAction(),
fEventAction(eventAction),
fScoringVolume(0),
fRunningAction(RunningAction),
fAbsHoleDiam(AbsHoleDiam),
fGaSet(GaSet),
fSourceChoice(SourceChoice)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
	if (fGaSet) {
		//to mute the "unused" warning, could be useful in the future!
	}
	G4VPhysicalVolume* ThisVol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4VPhysicalVolume* NextVol = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	
	G4int debug=0;
	
#pragma mark Annihilation
	//In case of Beta+ source look for annihilation points
	if (fSourceChoice==4 && step->GetTrack()->GetDynamicParticle()->GetPDGcode() == -11 && step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="annihil") {
		G4Event* evt = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
		evt->KeepTheEvent();
		(fRunningAction->GetAnnihX()).push_back(step->GetPostStepPoint()->GetPosition().x()/mm);
		(fRunningAction->GetAnnihY()).push_back(step->GetPostStepPoint()->GetPosition().y()/mm);
		(fRunningAction->GetAnnihZ()).push_back(step->GetPostStepPoint()->GetPosition().z()/mm);
		(fRunningAction->GetAnnihT()).push_back(step->GetPostStepPoint()->GetLocalTime()/ns);
	}
	
#pragma mark Optical Photons
	// ################################################################################
	// ###################### Optical Photons ENTERING SiPm
	if(step->GetTrack()->GetDynamicParticle() ->GetPDGcode()== 0 && NextVol && ThisVol->GetName()=="Pter" && NextVol->GetName()=="SiPm") {
		//		G4cout<<"FOTONE OTTICO ENTRA IN SiPm"<<G4endl;
		fEventAction->AddNPMT(1);
	}
	// ######################
	// ################################################################################
	
	//TODO: test
#pragma mark SiPm
	// ################################################################################
	// ###################### Interactions in SiPm
	if (ThisVol->GetName()=="SiPm") {
		fEventAction->AddEdepSiPM(step->GetTotalEnergyDeposit());
		//		G4cout<<"INTERAZIONE NEL SIPM: DepEne= "<<step->GetTotalEnergyDeposit()/keV<<" Part= "<<step->GetTrack()->GetDynamicParticle() ->GetPDGcode()<<G4endl;
	}
	// ######################
	// ################################################################################
	
#pragma mark Entering Pter
	// ################################################################################
	// ###################### ENTERING Pter (from wherever)
	if((NextVol && ThisVol->GetName()!="Pter" && NextVol->GetName()=="Pter")) { //what enters Pter (form every different volume)
		
		if (debug) G4cout<<"\nCIAODEBUG\n Particella entrata in PTER - fEventAction->GetEnteringParticle() ERA = "<<fEventAction->GetEnteringParticle();
		fEventAction->SetEnteringParticle(step->GetTrack()->GetDynamicParticle()->GetPDGcode()); //TODO: sistemare per evitare errori dovuti al fatto che è l'ultima traccia creata la prima ad essere tracciata..
		if (debug) G4cout<<" SETTO fEventAction->GetEnteringParticle()= "<<fEventAction->GetEnteringParticle()<<G4endl<<G4endl;
		
		// Check if the current track had already enetered somehow the PTER (to avoid double counting), otherwise increase the counter
		if (fEventAction->GetPterStoreTrackID()==step->GetTrack()->GetTrackID()) { //if I already saw this track exiting the source...
			fEventAction->AddPterPassCounter(1);  //increase the counter
			//			G4cout<<"PterDEBUG CONTROLLA "<<fEventAction->GetPterStoreTrackID()<<", PassCounter= "<<fEventAction->GetPterPassCounter()<<G4endl;
		}else {
			fEventAction->SetPterStoreTrackID(step->GetTrack()->GetTrackID());
			//			G4cout<<"PterDEBUG PRIMO PASSAGGIO!! "<<fEventAction->GetPterStoreTrackID()<<", PassCounter= "<<fEventAction->GetPterPassCounter()<<G4endl;
			//            if (fEventAction->GetPassCounter()!=0) G4cout<<"MERDAAAAA Primo passaggio di"<<fEventAction->GetStoreTrackID()<<" ma con PassCounter= "<<fEventAction->GetPassCounter()<<G4endl;
		}
		
		// Salvo le info solo della prima volta che una particella entra in pter
		if (fEventAction->GetPterPassCounter()==0) {
			(fRunningAction->GetPrePterEn()).push_back(step->GetPostStepPoint()->GetKineticEnergy()/keV);
			fEventAction->AddPrePterNo(1); //update the counter of particles entering Pter in the event
			(fRunningAction->GetPrePterPart()).push_back(step->GetTrack()->GetDynamicParticle()->GetPDGcode()); //add PID of particle enetering Pter
		}
	}
	// ###################### END ENTERING Pter
	// ################################################################################
	
#pragma mark Exiting Source
	// ################################################################################
	// ###################### EXITING SOURCE
	if( NextVol && ThisVol->GetName()!="DummyExitSorg" && NextVol->GetName()=="DummyExitSorg" && step->GetPreStepPoint()->GetMomentumDirection().z()>0) //New (2019.01.16) logic: if I go from somewhere into DummyExitSorg with a negative Z direction..
	{
		//collamaf: to avoid double counting same track going back and forth, check if I already counted it
		if (fEventAction->GetSourceExitStoreTrackID()==step->GetTrack()->GetTrackID()) { //if I already saw this track exiting the source...
			fEventAction->AddSourceExitPassCounter(1);  //increase the counter: number of times that this track exits the source (will not be written in scoring, but used to check if it is the first crossing
		}else {
			fEventAction->SetSourceExitStoreTrackID(step->GetTrack()->GetTrackID());
		}
		
		// Salvo le info solo della prima volta che una particella esce dalla sorgente
		if (fEventAction->GetSourceExitPassCounter()==0) {
			fEventAction->AddNSourceExit(1); //contatore di quante tracce escono dalla sorgente
			(fRunningAction->GetExitEn()).push_back(step->GetPostStepPoint()->GetKineticEnergy()/keV);
			(fRunningAction->GetExitX()).push_back(step->GetPostStepPoint()->GetPosition().x()/mm);
			(fRunningAction->GetExitY()).push_back(step->GetPostStepPoint()->GetPosition().y()/mm);
			(fRunningAction->GetExitZ()).push_back(step->GetPostStepPoint()->GetPosition().z()/mm);
			(fRunningAction->GetExitCosX()).push_back(step->GetPreStepPoint()->GetMomentumDirection().x());
			(fRunningAction->GetExitCosY()).push_back(step->GetPreStepPoint()->GetMomentumDirection().y());
			(fRunningAction->GetExitCosZ()).push_back(step->GetPreStepPoint()->GetMomentumDirection().z());
			(fRunningAction->GetExitPart()).push_back(step->GetTrack()->GetDynamicParticle()->GetPDGcode());
			(fRunningAction->GetExitParentID()).push_back(step->GetTrack()->GetParentID());
			
			if (step->GetTrack()->GetCreatorProcess()) {
				(fRunningAction->GetExitProcess().push_back((step->GetTrack()->GetCreatorProcess()->GetProcessType())));
			} else {
				(fRunningAction->GetExitProcess().push_back(-17));
			}
		}
	}
	// ###################### END EXITING SOURCE
	// ################################################################################
	
#pragma mark Exiting Absorber
	// ################################################################################
	// ###################### EXITING ABSORBER (if any)
	if(fAbsHoleDiam>=0 && NextVol && ThisVol->GetName()=="Absorber" && NextVol->GetName()=="DummyExitAbs") {
		
		//to avoid double counting same track going back and forth, check if I already counted it
		if (fEventAction->GetPostAbsStoreTrackID()==step->GetTrack()->GetTrackID()) { //if I already saw this track exiting the absorber...
			fEventAction->AddPostAbsPassCounter(1);  //increase the counter
		}else {
			fEventAction->SetPostAbsStoreTrackID(step->GetTrack()->GetTrackID());
		}
		
		// Salvo le info solo della prima volta che una particella esce dall'assorbitore
		if (fEventAction->GetPostAbsPassCounter()==0) {
			fEventAction->AddPostAbsNo(1); //update the counter of particles entering Pter in the event
			(fRunningAction->GetPostAbsEn()).push_back(step->GetPostStepPoint()->GetKineticEnergy()/keV);
			(fRunningAction->GetPostAbsPart()).push_back(step->GetTrack()->GetDynamicParticle()->GetPDGcode());
		}
	}
	// ###################### END EXITING ABSORBER (if any)
	// ################################################################################
	
	
#pragma mark Entering Probe
	// ################################################################################
	// ###################### ENTERING PROBE
	if(NextVol && ThisVol->GetName()=="DummyEnterProbe" && NextVol->GetName()=="FrontShield") {
		
		//to avoid double counting same track going back and forth, check if I already counted it
		if (fEventAction->GetPreProbeStoreTrackID()==step->GetTrack()->GetTrackID()) { //if I already saw this track exiting the absorber...
			fEventAction->AddPreProbePassCounter(1);  //increase the counter
		}else {
			fEventAction->SetPreProbeStoreTrackID(step->GetTrack()->GetTrackID());
		}
		
		// Salvo le info solo della prima volta che una particella esce dall'assorbitore
		if (fEventAction->GetPreProbePassCounter()==0) {
			fEventAction->AddPreProbeNo(1); //update the counter of particles entering Pter in the event
			(fRunningAction->GetPreProbeEn()).push_back(step->GetPostStepPoint()->GetKineticEnergy()/keV);
			(fRunningAction->GetPreProbePart()).push_back(step->GetTrack()->GetDynamicParticle()->GetPDGcode());
		}
	}
	// ###################### END ENTERING PROBE
	// ################################################################################
	
#pragma mark Inside Pter
	//Retrieve scoring volume
	if (!fScoringVolume) {
		const B1DetectorConstruction* detectorConstruction
		= static_cast<const B1DetectorConstruction*>
		(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		fScoringVolume = detectorConstruction->GetScoringVolume();
	}
	
	// ########################################
	// ###################### INSIDE Pter - Per each hit into sensitive detector
	// check if we are in scoring volume
	if (ThisVol->GetLogicalVolume() == fScoringVolume && step->GetTrack()->GetDynamicParticle() ->GetPDGcode()!=0 ) {
		fEventAction->SetEnterPterFlag(); //Something entered pter in this event
		fEventAction->AddNumHitsDet(1); //Update the counter of number of interactions in detector
		
		// collect energy deposited in this step
		G4double edepStep = step->GetTotalEnergyDeposit();
		
		//Fill vector
		(fRunningAction->GetPterEn()).push_back(step->GetTotalEnergyDeposit()/keV);
		(fRunningAction->GetPterEnPrim()).push_back(fRunningAction->GetMotherEnergy());
		(fRunningAction->GetPterPartPrim()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		//		(fRunningAction->GetPterTime()).push_back(step->GetTrack()->GetLocalTime()/ns);
		(fRunningAction->GetPterTime()).push_back(step->GetTrack()->GetGlobalTime()/ns-fRunningAction->GetMotherTime());
		//		G4cout<<"PterDEBUG  MotherTime= "<< fRunningAction->GetMotherTime()<<" PostDiff= "<<  step->GetTrack()->GetGlobalTime()/ns-fRunningAction->GetMotherTime() <<G4endl;
		(fRunningAction->GetPterX()).push_back(step->GetPreStepPoint()->GetPosition().x()/mm);
		(fRunningAction->GetPterY()).push_back(step->GetPreStepPoint()->GetPosition().y()/mm);
		(fRunningAction->GetPterZ()).push_back(step->GetPreStepPoint()->GetPosition().z()/mm);
		(fRunningAction->GetPterPart()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		
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
