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
// $Id: B1PrimaryGeneratorAction.hh 90623 2015-06-05 09:24:30Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.hh
/// \brief Definition of the B1PrimaryGeneratorAction class

#ifndef B1PrimaryGeneratorAction_h
#define B1PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "B1EventAction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"

class G4ParticleGun;
class G4Event;
class G4Tubs;
class B1DetectorConstruction;

/// The primary generator action class with particle gun.
///

class B1PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	//    B1PrimaryGeneratorAction(B1EventAction* eventAction);
	B1PrimaryGeneratorAction(B1DetectorConstruction *, B1EventAction *eventAction, G4double TBR = 1, G4int SourceSelect = 1, G4double SourceDiameter = 10, G4double SourceThickness = 7, G4int GaSetting = 1, G4double CaseDepth = -50, G4String ExtSourceFile = "");
	virtual ~B1PrimaryGeneratorAction();

	// method from the base class
	virtual void GeneratePrimaries(G4Event *);

	// method to access particle gun
	const G4ParticleGun *GetParticleGun() const { return fParticleGun; }

	G4double BetaDecaySpectrum(G4double Ek, G4double EndPoint);

	inline void SetRadius(const double Radius) { fRadiusMax = Radius; };
	inline void SetZ(const double Z) { fZ = Z; };

private:
	B1DetectorConstruction *fDetector;

	G4ParticleGun *fParticleGun; // pointer a to G4 gun class
	B1EventAction *evtPrimAction;
	//	G4Tubs * fSourceVol;

	G4double fZ;

	G4double fRadiusExt;
	G4double fRadiusInt;
	G4double fDZExt;
	G4double fDZInt;

	G4double fRadiusMax;
	G4double fRadiusMin;
	G4double fTBR;
	//	G4double fRatio;
	G4int fSourceSelect;
	G4double fSourceDiameter;
	G4double fSourceThickness;
	G4int fGaSet;
	G4bool FlatEle = false;
	G4bool FlatGamma = false;
	G4double fCaseDepth;
	//	G4LogicalVolume* sourceLV;
	//	G4VPhysicalVolume* sourcePV;
	G4VPrimaryGenerator *hepmcAscii;
	G4String fExtSourceFile;

	double spessoreFondoContenitoreFluoro = 10 * mm;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
