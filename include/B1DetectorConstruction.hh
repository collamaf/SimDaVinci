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
// $Id: B1DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file B1DetectorConstruction.hh
/// \brief Definition of the B1DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Region.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
public:
	B1DetectorConstruction(G4double, G4double, G4double, G4int, G4int ,G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double , G4double, G4double, G4int, G4bool, G4int, G4int, G4int,G4double, G4bool);
	virtual ~B1DetectorConstruction();
		
	virtual G4VPhysicalVolume* Construct();
	
	G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
	
	double ContainerOuterRadius;
	double ContainerMaxHeigth;

protected:
	G4LogicalVolume*  fScoringVolume;
	G4double fX0Scan;
	G4double fZValue;
	G4double fAbsHoleDiam;
	G4int fSourceSelect;
	G4int fAbsorberMaterial;
	G4double fPterDiameter;
	G4double fPterThickness;
	G4double fSourceDiameter;
	G4double fSourceThickness;
	G4double fAbsorberThickness;
	G4double fCaseDepth;
	G4double fLateralCaseThickness;
	G4double fBackCaseThickness;
	G4double fHorsesShoeLateralThickness;
	G4double fHorsesShoeBackThickness;
	G4int fHousingCase;
	G4bool fScintFlag;
	G4int fGaSet;
	G4int fApparatusMat;
	G4int fPosAbsorber;
	G4double fAbsCenter;
	G4bool fSecondShieldFlag;
	
	G4Region* sorgente = new G4Region("SourceReg");
	G4Region* ABSRegion = new G4Region("ABSRegion");
	G4Region* frontshieldreg = new G4Region("FrontShieldReg");
	G4Region* pterreg = new G4Region("PTERReg");
	G4Region* carrierreg = new G4Region("CarrierReg");


	
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


