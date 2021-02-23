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
// $Id: B1ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B1ActionInitialization.cc
/// \brief Implementation of the B1ActionInitialization class

#include "B1ActionInitialization.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1RunAction.hh"
#include "B1EventAction.hh"
#include "B1SteppingAction.hh"
#include "B1StackingAction.hh"
#include "B1DetectorConstruction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ActionInitialization::B1ActionInitialization(B1DetectorConstruction* det, G4double x0, G4double ZValue, G4double AbsHoleDiam, G4double TBR, G4int SourceSelect, G4int SensorChoice, G4double SourceDiameter, G4double SourceThickness, G4String FileName, G4int GaSet, G4double CaseDepth, G4String ExtSourceFile, G4bool LightOutFlag)
  : G4VUserActionInitialization(), fDetector(det),fX0Scan(x0), fZValue(ZValue), fAbsHoleDiam(AbsHoleDiam), fTBR(TBR),  	fSourceSelect(SourceSelect), fSensorChoice(SensorChoice), fSourceDiameter(SourceDiameter) ,fSourceThickness(SourceThickness), fFileName(FileName), fGaSet(GaSet), fCaseDepth(CaseDepth), fExtSourceFile(ExtSourceFile), fLightOutFlag(LightOutFlag)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ActionInitialization::~B1ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ActionInitialization::BuildForMaster() const
{
  B1RunAction* runAction = new B1RunAction(fX0Scan, fZValue, fAbsHoleDiam, fTBR, fSourceSelect, fSensorChoice, fFileName);
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ActionInitialization::Build() const
{
//  SetUserAction(new B1PrimaryGeneratorAction(runAction));
	//G4cout<<"PROVA Action Init "<<fX0Scan<<G4endl;
	//G4cout<<"GaSetting "<<fGaSet<<G4endl;


  B1RunAction* runAction = new B1RunAction(fX0Scan, fZValue, fAbsHoleDiam, fTBR, fSourceSelect, fSensorChoice, fFileName);
  SetUserAction(runAction);
  
  B1EventAction* eventAction = new B1EventAction(runAction,fLightOutFlag);
  SetUserAction(eventAction);
	
  SetUserAction(new B1SteppingAction(eventAction, runAction, fAbsHoleDiam,fGaSet, fSourceSelect));
	
//	B1PrimaryGeneratorAction* primAction= new B1PrimaryGeneratorAction(eventAction, TRUE, fSrSourceFlag, TRUE, fTBR, fSrSourceFlag); // Y, Sr, PrintDist, TBR sorge
	B1PrimaryGeneratorAction* primAction= new B1PrimaryGeneratorAction(fDetector,eventAction,  fTBR, fSourceSelect, fSourceDiameter, fSourceThickness,fGaSet, fCaseDepth, fExtSourceFile);
	SetUserAction(primAction);
	SetUserAction(new B1StackingAction(runAction, eventAction, fSourceSelect));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
