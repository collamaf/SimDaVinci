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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class
///
///
///

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4StepLimiter.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction(G4double x0, G4double ZValue, G4double CuDiam, G4int SourceSelect, G4int AbsorberMaterial,G4double PterDiameter, G4double PterThickness,G4double SourceDiameter,G4double SourceThickness, G4double AbsorberThickness, G4double ProbeCaseDepth, G4double ProbeCaseLateralThickness, G4double ProbeCaseBackThickness, G4double HSLateralThickness, G4double HSBackThickness)
: G4VUserDetectorConstruction(),
fScoringVolume(0), fX0Scan(x0), fZValue(ZValue), fCuDiam(CuDiam), fSourceSelect(SourceSelect), fAbsorberMaterial(AbsorberMaterial), fPterDiameter(PterDiameter), fPterThickness(PterThickness), fSourceDiameter(SourceDiameter), fSourceThickness(SourceThickness), fAbsorberThickness(AbsorberThickness),fCaseDepth(ProbeCaseDepth),fLateralCaseThickness(ProbeCaseLateralThickness), fBackCaseThickness(ProbeCaseBackThickness), fHorsesShoeLateralThickness(HSLateralThickness),fHorsesShoeBackThickness(HSBackThickness)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	
	// Option to switch on/off checking of volumes overlaps
	//
	G4bool checkOverlaps = false;
	
	//
	// World
	//
	G4double world_sizeXY = 0.5*m;
	G4double world_sizeZ  = 0.5*m;
	G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
	//	G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
	
	G4Box* solidWorld =
	new G4Box("World",                       //its name
						0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
	
	G4LogicalVolume* logicWorld =
	new G4LogicalVolume(solidWorld,          //its solid
											world_mat,           //its material
											"World");            //its name
	
	G4VPhysicalVolume* physWorld =
	new G4PVPlacement(0,                     //no rotation
										G4ThreeVector(),       //at (0,0,0)
										logicWorld,            //its logical volume
										"World",               //its name
										0,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	//###################################################################
	//###################################################
	// Definitions of materials
	//##########################
	
	G4double z, a, density;
	G4String name, symbol;
	G4int ncomponents, natoms;
	
	a = 1.01*g/mole;
	G4Element* elH = new G4Element (name="Hydrogen", symbol="H", z=1.,a );
	a = 12.01*g/mole;
	G4Element* elC = new G4Element (name="Carbon", symbol="C", z=6.,a );
	a = 16.00*g/mole;
	G4Element* elO = new G4Element (name="Oxygen", symbol="O", z=8.,a );
	a = 14.00*g/mole;
	G4Element* elN = new G4Element (name="Nitrogen", symbol="N", z=7.,a );
	
	
	density = 4.000*g/cm3; //4 for MT9V011, 2.43 for MT9V115
	//if (fSensorChoice==2) density=2.43;
	G4Material* FrontShield = new G4Material (name="FrontShield", density, ncomponents=3);
	FrontShield->AddElement (elH, natoms=30);
	FrontShield->AddElement (elC, natoms=20);
	FrontShield->AddElement (elO, natoms=2);
	
	G4double densityAlu = 2.600*g/cm3;
	//	G4NistManager* man = G4NistManager::Instance();
	G4NistManager::Instance()->BuildMaterialWithNewDensity("MyAlu","G4_Al",densityAlu);
	
	/* ref from TestEm7
	 G4NistManager::Instance()->
	 BuildMaterialWithNewDensity("Water_1.05","G4_WATER",1.05*g/cm3);
	 */
	
	//G4int TypeSourceFlag=fSourceSelect;
	
	//###################################################
	// AGAR AGAR Source - AgarAgar should be C14 H24 O9
	//##########################
	
	G4double Agardensity = 1.030*g/cm3;
	G4Material* AgarAgar = new G4Material (name="AgarAgar", Agardensity, ncomponents=3);
	AgarAgar->AddElement (elH, natoms=24);
	AgarAgar->AddElement (elC, natoms=14);
	AgarAgar->AddElement (elO, natoms=9);
	
	
	//###################################################
	// ABS material - ABS should be C15 H17 N
	//##########################
	G4double ABSdensity = 0.7*g/cm3;
	G4Material* ABS = new G4Material (name="ABS", ABSdensity, ncomponents=3);
	ABS->AddElement (elH, natoms=17);
	ABS->AddElement (elC, natoms=15);
	ABS->AddElement (elN, natoms=1);
	
	//###################################################
	// P-Terphenyl Material
	//##########################
	
	G4double PTerphenyldensity = 1.23*g/cm3;
	G4Material* PTerphenyl= new G4Material (name="PTerphenyl", PTerphenyldensity, ncomponents=2);
	PTerphenyl->AddElement (elC, natoms=18);
	PTerphenyl->AddElement (elH, natoms=14);
	
	//###################################################
	// Delrin Material      C H2 O
	//##########################
	
	G4double Delrindensity = 1.41*g/cm3;
	G4Material* Delrin= new G4Material (name="Delrin", Delrindensity, ncomponents=3);
	Delrin->AddElement (elC, natoms=1);
	Delrin->AddElement (elH, natoms=2);
	Delrin->AddElement (elO, natoms=1);

	
	
	//############ MATERIAL ASSIGNMENT
	G4Material* SourceExtY_mat = AgarAgar;
	G4Material* ABSaround_mat = ABS;
	G4Material* ABSbehind_mat = ABS;
	G4Material* GaContainer_mat=nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	//G4Material* GaWater_mat=nist->FindOrBuildMaterial("G4_WATER");
	G4Material* SourceExtGa_mat=nist->FindOrBuildMaterial("G4_WATER");
	G4Material* SourceSR_mat = nist->FindOrBuildMaterial("MyAlu"); //G4_Al
	G4Material* FrontShield_mat = nist->FindOrBuildMaterial("MyAlu");
	G4Material* shapeDummy_mat = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* Pter_mat = PTerphenyl;
	G4Material* PVC_mat= nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material* Delrin_mat=Delrin;
	G4Material* shapeCo_mat = nist->FindOrBuildMaterial("G4_Cu");
	//G4Material* ProbeCase_mat = nist->FindOrBuildMaterial("MyAlu");
	//G4Material* CylinderB_mat = nist->FindOrBuildMaterial("MyAlu");
	
	G4Material* PlasticCase_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material* HorsesShoe_mat= nist->FindOrBuildMaterial("G4_Pb");
	G4Material* CaseInner_mat=nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material* MiddleCase_mat=CaseInner_mat;
	G4Material* TopCase_mat=CaseInner_mat;
	
	HorsesShoe_mat=world_mat;
	CaseInner_mat=world_mat;

	
	
	//Before Boolean
	/*
	 G4Material* CaseExt_mat=nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	 G4Material* CaseMetal_mat=nist->FindOrBuildMaterial("G4_Pb");
	 G4Material* CaseInner_mat=nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	 
	 G4Material* AroundBackMiddleCase_mat=CaseExt_mat;
	 G4Material* AroundTopCase_mat=CaseExt_mat;
	 G4Material* EndCase_mat=CaseExt_mat;
	 G4Material* ExtMiddleCase_mat=CaseExt_mat;
	 
	 G4Material* MiddleCase_mat=CaseInner_mat;
	 G4Material* TopCase_mat=CaseInner_mat;
	 G4Material* BackMiddleCase_mat=CaseMetal_mat;
	 G4Material* AroundMiddleCase_mat=CaseMetal_mat;
	 */
	
#if 0
	G4Material* TopCase_mat = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* AroundTopCase_mat = nist->FindOrBuildMaterial("MyAlu");
	//G4Material* ExternalCase_mat = nist->FindOrBuildMaterial("MyAlu");
	G4Material* MiddleCase_mat = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* AroundMiddleCase_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");  //Understand what type of material
	G4Material* ExtMiddleCase_mat = nist->FindOrBuildMaterial("MyAlu");
	G4Material* BackMiddleCase_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");  //Understand what type of material
	G4Material* AroundBackMiddleCase_mat = nist->FindOrBuildMaterial("MyAlu");
	G4Material* EndCase_mat = nist->FindOrBuildMaterial("MyAlu");
	
#endif
	
	// EndCase, AroundBackMiddleCase and ExtMiddleCase are made by the same material
	
  // AroundMiddleCase and BackMiddleCase are made by the same material
	
	
	
	//###################################################################
	//###################################################
	// Definitions of dimensions and sizes
	//##########################
	
	//### ExtY SOURCE
	G4double RminSourceExtY = 0.*mm;
	G4double RmaxSourceExtY = 10.5*mm; //10.48 per Rosa, 6.65 per PG
	G4double DzSourceExtY= 4.5*mm; //4.4 per Rosa, 5.5 per PG
	G4double SPhiSourceExtY = 0.*deg;
	G4double DPhiSourceExtY = 360.*deg;
	//###
	
	//### ABS
	G4double RminABSaround = RmaxSourceExtY;
	G4double RmaxABSaround = 12.5*mm;
	G4double DzABSaround= DzSourceExtY;
	G4double SPhiABSaround = 0.*deg;
	G4double DPhiABSaround = 360.*deg;
	
	G4double RminABSbehind = 0.*mm;
	G4double RmaxABSbehind = RmaxABSaround;
	G4double DzABSbehind= 3*mm;
	G4double SPhiABSbehind = 0.*deg;
	G4double DPhiABSbehind = 360.*deg;
	//###
	
	
	//### Sr Source
	G4double RminSourceSR = 0.*mm;
	G4double RmaxSourceSR = 12.5*mm; //physical dimensions same for PG/RM sources, the active one differs
	G4double DzSourceSR= 3*mm;
	G4double SPhiSourceSR = 0.*deg;
	G4double DPhiSourceSR = 360.*deg;
	//###
	
	//### Ga Source Container
	G4double rSourceExtGa = fSourceDiameter*0.5*mm;
	G4double RContainer = 50*mm;
	G4double dzSourceExtGa= fSourceThickness*mm;
	G4double DzContainer= 20*mm;
	G4double SPhiSourceExtGa = 0.*deg;
	G4double DPhiSourceExtGa = 360.*deg;
	
	
	//### Filter (FrontShield)
	G4double Z_FrontShield= 0*mm;
	G4double FrontShield_outer_r=(12.0/2.0)*mm;
	G4double FrontShield_sizeZ=5*um;
	G4double FrontShield_start_angle=0.*deg;
	G4double FrontShield_spanning_angle=360.0*deg;
	//###
	
	//### Copper Collimator
	G4double RminCo = fabs(fCuDiam)/2.*mm;
	G4double RmaxCo = 18.*mm;
	G4double DzCo= fAbsorberThickness*mm;
	G4double SPhiCo = 0.*deg;
	G4double DPhiCo = 360.*deg;
	//###
	
	//### Dummy
	G4double RminDummy = 0.*mm;
	G4double RmaxDummy = 18.*mm;
	if (fCuDiam>=0) RmaxDummy =RmaxCo;
	else RmaxDummy=RmaxSourceSR;
	G4double DzDummy= 1.e-5*mm;
	G4double SPhiDummy = 0.*deg;
	G4double DPhiDummy = 360.*deg;
	G4double zDummy;
	//###
	
	
	//### Pter
	
	G4double Pter_Diam=fPterDiameter*mm;
	G4double PVC_outer_r=12.0*mm/2.;
	G4double Pter_sizeZ=fPterThickness*mm;
	G4double Pter_start_angle=0.*deg;
	G4double Pter_spanning_angle=360.0*deg;
	G4double Pter_Posz=0.*mm;
	G4double Pter_ZScan=fZValue*mm;
	
	G4double PVC_inner_r= PVC_outer_r - 2*mm;

	
	//### Probe Case
	
	G4double CaseDepth = fCaseDepth*mm;
	//G4double SPhiCase = 0.*deg;
	//G4double DPhiCase = 360.*deg;
	G4double LateralCaseThickness = fLateralCaseThickness*mm;
	G4double BackCaseThickness = fBackCaseThickness*mm;
	//G4double ProbeCase_Posz=0.*mm;
	G4double TopCaseDepth=2.*mm;
	G4double SPhiTopCase = 0.*deg;
	G4double DPhiTopCase = 360.*deg;
	G4double SPhiAroundTopCase = 0.*deg;
	G4double DPhiAroundTopCase = 360.*deg;
	//G4double Case_r = PVC_outer_r +0.1*mm;
	G4double HorsesShoeLateralThickness = fHorsesShoeLateralThickness * mm;
	G4double HorsesShoeBackThickness = fHorsesShoeBackThickness * mm;

	
	
	//##########################
	//###################################################
	
	
	
	//###################################################################
	//###################################################
	// Definitions of volumes
	//##########################
	//###################################################################
	
	
	//###################################################
	// ExtY Source
	//##########################
	G4ThreeVector posSourceExtY = G4ThreeVector(0, 0, -DzSourceExtY*0.5);
	
	G4Tubs* solidSourceExtY =
	new G4Tubs("SourceExtY",                       //its name
						 RminSourceExtY,
						 RmaxSourceExtY,
						 0.5*DzSourceExtY,
						 SPhiSourceExtY,
						 DPhiSourceExtY);     //its size
	
	G4LogicalVolume* logicSourceExtY =
	new G4LogicalVolume(solidSourceExtY,          //its solid
											SourceExtY_mat,           //its material
											"SourceExtY");            //its name
	
	if(fSourceSelect==3) { //I place the ExtY source if I am not asking for Sr source

		G4cout<<"GEOMETRY DEBUG - Z thickness of solidSourceExtY= "<<DzSourceExtY/mm<<", Z pos= "<<-DzSourceExtY*0.5<<G4endl;
		
		G4cout<<"GEOMETRY DEBUG - ExtYTOC Source has been placed!!"<<G4endl;
		
		new G4PVPlacement(0,                     //no rotation
											posSourceExtY,       //at (0,0,0)
											logicSourceExtY,            //its logical volume
											"SourceExtY",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		//G4Region* sorgente = new G4Region("SourceReg");
		logicSourceExtY->SetRegion(sorgente);
		sorgente->AddRootLogicalVolume(logicSourceExtY);
	}
	//################################################### END ExtY SOURCE
	
	
	//###################################################
	// ABS carrier around ExtY Source
	//##########################
	G4ThreeVector posABSaround = G4ThreeVector(0, 0, -DzABSaround*0.5);
	
	G4Tubs* solidABSaround =
	new G4Tubs("ABSaround",                       //its name
						 RminABSaround,
						 RmaxABSaround,
						 0.5*DzABSaround,
						 SPhiABSaround,
						 DPhiABSaround);     //its size
	
	G4LogicalVolume* logicABSaround =
	new G4LogicalVolume(solidABSaround,          //its solid
											ABSaround_mat,           //its material
											"ABSaround");            //its name
	
	if(fSourceSelect==3) {  //I place the ABS carrier of the ExtY source if I am not asking for Sr source
		G4cout<<"GEOMETRY DEBUG - Z thickness of solidABSaround= "<<DzABSaround/mm<<", Z pos= "<<-DzABSaround*0.5<<G4endl;
		
		G4cout<<"GEOMETRY DEBUG - ExtYTOC Source has been placed!!"<<G4endl;
		
		new G4PVPlacement(0,                     //no rotation
											posABSaround,       //at (0,0,0)
											logicABSaround,            //its logical volume
											"ABSaround",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
	}
	//G4Region* ABSRegion = new G4Region("ABSRegion");
	logicABSaround->SetRegion(ABSRegion);
	ABSRegion->AddRootLogicalVolume(logicABSaround);
	
	//################################################### END ABS AROUND
	
	//###################################################
	// ABS carrier behind ExtY Source
	//##########################
	G4ThreeVector posABSbehind = G4ThreeVector(0, 0, -DzABSbehind*0.5- DzABSaround);
	
	G4Tubs* solidABSbehind =
	new G4Tubs("ABSbehind",                       //its name
						 RminABSbehind,
						 RmaxABSbehind,
						 0.5*DzABSbehind,
						 SPhiABSbehind,
						 DPhiABSbehind);     //its size
	
	G4LogicalVolume* logicABSbehind =
	new G4LogicalVolume(solidABSbehind,          //its solid
											ABSbehind_mat,           //its material
											"ABSbehind");            //its name
	
	if(fSourceSelect==3) { //I place the ABS carrier of the ExtY source if I am not asking for Sr source
		G4cout<<"GEOMETRY DEBUG - Z thickness of solidABSbehind= "<<DzABSbehind/mm<<", Z pos= "<<-DzABSbehind*0.5- DzABSaround<<G4endl;
		
		G4cout<<"GEOMETRY DEBUG - ExtYTOC Source has been placed!!"<<G4endl;
		
		new G4PVPlacement(0,                     //no rotation
											posABSbehind,          //at (0,0,0)
											logicABSbehind,        //its logical volume
											"ABSbehind",           //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
	}
	logicABSbehind->SetRegion(ABSRegion);
	ABSRegion->AddRootLogicalVolume(logicABSbehind);
	
	//################################################### END ABS BEHIND
	
	
	//###################################################
	// Sr90 lab Source
	//##########################
	G4ThreeVector posSourceSR = G4ThreeVector(0, 0, -DzSourceSR*0.5);
	
	G4cout<<"GEOMETRY DEBUG - Z thickness of solidSourceSR= "<<DzSourceSR/mm<<", Z pos= "<<-DzSourceSR*0.5<<G4endl;
	
	G4Tubs* solidSourceSR =
	new G4Tubs("SourceSR",                       //its name
						 RminSourceSR,
						 RmaxSourceSR,
						 0.5*DzSourceSR,
						 SPhiSourceSR,
						 DPhiSourceSR);     //its size
	
	G4LogicalVolume* logicSourceSR =
	new G4LogicalVolume(solidSourceSR,          //its solid
											SourceSR_mat,           //its material
											"SourceSR");            //its name
	
	if(fSourceSelect==1 || fSourceSelect==2) { //If i requested the Sr source
		G4cout<<"GEOMETRY DEBUG - Sr Source has been placed!!"<<G4endl;
		
		new G4PVPlacement(0,                     //no rotation
											posSourceSR,       //at (0,0,0)
											logicSourceSR,            //its logical volume
											"SourceSR",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		//G4Region* sorgente = new G4Region("SourceReg");
		logicSourceSR->SetRegion(sorgente);
		sorgente->AddRootLogicalVolume(logicSourceSR);
	}
	//################################################### END SR SOURCE
	
	
	//###################################################
	//Ga Source Container
	//##########################
	
	G4ThreeVector posContainerExtGa = G4ThreeVector(0, 0, -DzContainer*0.5);
	G4ThreeVector posExtGa = G4ThreeVector(0, 0, -dzSourceExtGa*0.5);

	
	//G4cout<<"GEOMETRY DEBUG - Z thickness of solidSourceSR= "<<DzSourceSR/mm<<", Z pos= "<<-DzSourceSR*0.5<<G4endl;
	
	//G4cout<<"GEOMETRY DEBUG Source Diameter "<<fSourceDiameter<<G4endl;
	//G4cout<<"GEOMETRY DEBUG Source thickness "<<fSourceThickness<<G4endl;
	G4cout<<"GEOMETRY DEBUG Absorber thickness "<<fAbsorberThickness<<G4endl;



	
	G4VSolid* GaCylinder =
	new G4Tubs("GaCylinder",                       //its name
						 0.,
						 RContainer,
						 0.5*DzContainer,
						 SPhiSourceExtGa,
						 DPhiSourceExtGa);     //its size
	
	
	G4VSolid* SourceExtGa =
	new G4Tubs("SourceExtGa",                       //its name
						 0.,
						 rSourceExtGa,
						 0.5*dzSourceExtGa,
						 SPhiSourceExtGa,
						 DPhiSourceExtGa);     //its size
	
	
	G4VSolid* GaContainer=
	new G4SubtractionSolid ("GaContainer",      //GaContainer=GaCylinder-GaSource
													GaCylinder,
													SourceExtGa,
													0,
													G4ThreeVector(0.,0.,DzContainer*0.5-dzSourceExtGa/2.));
	
	
	
	G4LogicalVolume* logicGaContainer =
	new G4LogicalVolume(GaContainer,               //its solid
											GaContainer_mat,           //its material
											"GaContainer");            //its name
	
	G4LogicalVolume* logicSourceExtGa =
	new G4LogicalVolume(SourceExtGa,               //its solid
											SourceExtGa_mat,           //its material
											"SourceExtGa");            //its name
	
	
	
	if(fSourceSelect==4) { //If i requested the Sr source
		//G4cout<<"GEOMETRY DEBUG - Sr Source has been placed!!"<<G4endl;
		
		new G4PVPlacement(0,                     //no rotation
											posContainerExtGa,       //at (0,0,0)
											logicGaContainer,            //its logical volume
											"GaContainer",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		new G4PVPlacement(0,                     //no rotation
											posExtGa,       //at (0,0,0)
											logicSourceExtGa,            //its logical volume
											"SourceExtGa",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		/*
		//G4Region* Container = new G4Region("ContainerReg");
		logicGaContainer->SetRegion(Container);
		Container->AddRootLogicalVolume(logicGaContainer);
		
		//G4Region* Water = new G4Region("WaterReg");
		logicGaWater->SetRegion(Water);
		Water->AddRootLogicalVolume(logicGaWater);
		 */
	}
	
	
	//################################################### END Ga SOURCE

	
	
	//###################################################
	//Copper Collimator
	//##########################
	
	G4ThreeVector posCo = G4ThreeVector(0, 0, DzCo*0.5);
	
	//G4cout<<"GEOMETRY DEBUG - Z thickness of solidShapeCo= "<<DzCo/mm<<", Z pos= "<<posCo.z()<<G4endl;
	
	G4Tubs* solidShapeCo =
	new G4Tubs("CuCollimator",                       //its name
						 RminCo,
						 RmaxCo,
						 0.5*DzCo,
						 SPhiCo,
						 DPhiCo);     //its size
	
	if(fAbsorberMaterial==1){
		shapeCo_mat = nist->FindOrBuildMaterial("G4_Cu");
	}else if(fAbsorberMaterial==2){
		shapeCo_mat = nist->FindOrBuildMaterial("G4_Pb");
	}else if(fAbsorberMaterial==3){
		shapeCo_mat = nist->FindOrBuildMaterial("MyAlu");
	}else if(fAbsorberMaterial==4){
		shapeCo_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	}
	
	G4cout<<"GEOMETRY DEBUG - absorber mat="<<shapeCo_mat<<G4endl;
	
	G4LogicalVolume* logicShapeCo =
	new G4LogicalVolume(solidShapeCo,          //its solid
											shapeCo_mat,           //its material
											"CuCollimator");            //its name
	
	if (fCuDiam>=0) {
		G4cout<<"GEOMETRY DEBUG - Copper collimator has been placed!!"<<G4endl;
		
		new G4PVPlacement(0,                     //no rotation
											posCo,       //at (0,0,0)
											logicShapeCo,            //its logical volume
											"CuCollimator",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		//		G4Region* sorgente = new G4Region("SourceReg");
		logicShapeCo->SetRegion(sorgente);
		sorgente->AddRootLogicalVolume(logicShapeCo);
		
	}
	
	//################################################### END OF COPPER COLLIMATOR
	
	
	
	//###################################################
	//Dummy volume for scoring what exits source
	//##########################
	
	if (fCuDiam<0) {
		zDummy=DzDummy*0.5;
	} else {
		zDummy=DzDummy*0.5+DzCo;
	}
	G4ThreeVector posDummy = G4ThreeVector(0, 0, zDummy);
	
	G4cout<<"GEOMETRY DEBUG - Z thickness of solidShapeDummy= "<<DzDummy/mm<<", Z pos= "<<zDummy<<G4endl;
	
	G4Tubs* solidShapeDummy =
	new G4Tubs("Dummy",                       //its name
						 RminDummy,
						 RmaxDummy,
						 0.5*DzDummy,
						 SPhiDummy,
						 DPhiDummy);     //its size
	
	G4LogicalVolume* logicShapeDummy =
	new G4LogicalVolume(solidShapeDummy,          //its solid
											shapeDummy_mat,           //its material
											"Dummy");            //its name
	
	G4cout<<"GEOMETRY DEBUG - Dummy volume has been placed!!"<<G4endl;
	
	new G4PVPlacement(0,                     //no rotation
										posDummy,       //at (0,0,0)
										logicShapeDummy,            //its logical volume
										"Dummy",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	//		G4Region* sorgente = new G4Region("SourceReg");
	logicShapeDummy->SetRegion(sorgente);
	sorgente->AddRootLogicalVolume(logicShapeDummy);
	
	//################################################### END OF DUMMY VOLUME
	
	
	
	
	//###################################################
	//Frontal part of the Probe
	//##########################
	
	
	
	
	//###################################################
	//Electron Filter FrontShield
	//##########################
	
	/*
	if (fSensorChoice==2) fFilterFlag=1; //Sensor 2 is always with filter
	if (fSensorChoice==3) fFilterFlag=0; //Sensor 3 is always with filter
																			 //	if (fFilterFlag==1) {
	//FrontShield_sizeX = noX*PixelSize*ScaleFactor;
	//FrontShield_sizeY = noY*PixelSize*ScaleFactor;
	*/
	
	
	
	Z_FrontShield = fZValue + FrontShield_sizeZ*0.5;
	
	/*
	if (fFilterFlag==0) { //if I do not want the filter, place it but make it thin and empty
		FrontShield_mat=world_mat;
	}*/
	
	G4ThreeVector posFilter = G4ThreeVector(fX0Scan, 0, Z_FrontShield);
	
	G4cout<<"GEOMETRY DEBUG - Z thickness of solidFrontShield= "<<FrontShield_sizeZ/mm<<", Z pos= "<<Z_FrontShield/mm<<G4endl;
	
	G4Tubs* solidFrontShield =
	new G4Tubs("FrontShield",                       //its name
						0.,FrontShield_outer_r,FrontShield_sizeZ*0.5,FrontShield_start_angle,FrontShield_spanning_angle);     //its size
	
	
	
	G4LogicalVolume* logicFrontShield =
	new G4LogicalVolume(solidFrontShield,          //its solid
											FrontShield_mat,           //its material
											"FrontShield");            //its name
	
	new G4PVPlacement(0,                     //no rotation
										posFilter,
										logicFrontShield,            //its logical volume
										"FrontShield",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
																					 //	G4Region* filtro = new G4Region("FrontShieldReg");
	
	logicFrontShield->SetRegion(frontshieldreg);
	frontshieldreg->AddRootLogicalVolume(logicFrontShield);

	
	//################################################### END OF FrontShield FILTER
	/*
	if(fSensorChoice==1) {
		Pter_ZScan=fZValue + FrontShield_sizeZ+Pter_sizeZ*0.5; //modified on 2017.11.21 by collamaf - Z distance does not take into account Cu thickness! is always from source top to possible resin
	}
	
	G4cout<<"GEOMETRY DEBUG - Z thickness of solidPter= "<<Pter_sizeZ/mm<<", Z pos= "<<Pter_ZScan/mm<<G4endl;
  //G4cout<<"GEOMETRY DEBUG - PterSizeX= "<<Pter_sizeX/mm<<", PterSizeY= "<<Pter_sizeY/mm<<", PterSizeZ= "<<pixZ/mm<<G4endl;
	
*/
	
	
	//###################################################
	// 	P-Terphenyl
	//##########################
	
	//G4cout<<"GEOMETRY DEBUG Pter Diameter "<<fPterDiameter<<G4endl;
	//G4cout<<"GEOMETRY DEBUG Pter thickness "<<fPterThickness<<G4endl;
	
	G4Tubs* solidPter =
	new G4Tubs("Pter",                                                                         //its name
						0.,Pter_Diam*0.5,Pter_sizeZ*0.5,Pter_start_angle,Pter_spanning_angle);                //its size
	
	G4LogicalVolume* logicPter =
	new G4LogicalVolume(solidPter,          //its solid
											Pter_mat,           //its material
											"Pter");            //its name
	
	
	// place detector-Pter in world

		Pter_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ*0.5;
		G4ThreeVector pos2 = G4ThreeVector(fX0Scan, 0, Pter_Posz);
		new G4PVPlacement(0,                     //no rotation
											pos2,
											logicPter,            //its logical volume
											"Pter",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
	
	
		//G4Region* Pterreg = new G4Region("PterReg");
		logicPter->SetRegion(pterreg);
		pterreg->AddRootLogicalVolume(logicPter);
	
		//Solid Si Pter
		fScoringVolume = logicPter;
	
	
	//###################################################
	// PVC around P-Terphenyl
	//##########################
	
	G4Tubs* solidPVC =
	new G4Tubs("PVC",                                                                         //its name
						 PVC_inner_r,PVC_outer_r,Pter_sizeZ*0.5,Pter_start_angle,Pter_spanning_angle);                //its size
	
	G4LogicalVolume* logicPVC =
	new G4LogicalVolume(solidPVC,          //its solid
											PVC_mat,              //its material
											"PVC");            //its name
	
	
	new G4PVPlacement(0,                     //no rotation
										pos2,
										logicPVC,            //its logical volume
										"PVC",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	/*
	G4Tubs* solidExternalCase =
	new G4Tubs("ExternalCase",                                                                         //its name
						 PVC_outer_r, Case_r , Pter_sizeZ*0.5,Pter_start_angle,Pter_spanning_angle);                //its size
	
	G4LogicalVolume* logicExternalCase =
	new G4LogicalVolume(solidExternalCase,          //its solid
											ExternalCase_mat,              //its material
											"ExternalCase");            //its name
	
	
	new G4PVPlacement(0,                     //no rotation
										pos2,
										logicExternalCase,            //its logical volume
										"ExternalCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	*/
	
	//###################################################
	// Delrin around P-Terphenyl
	//##########################
	

	
	G4Tubs* solidDelrin =
	new G4Tubs("Delrin",                                                                         //its name
						 Pter_Diam*0.5,PVC_inner_r,Pter_sizeZ*0.5,Pter_start_angle,Pter_spanning_angle);                //its size

	
	G4LogicalVolume* logicDelrin =
	new G4LogicalVolume(solidDelrin,          //its solid
											Delrin_mat,           //its material
											"Delrin");            //its name
	
	
		new G4PVPlacement(0,                     //no rotation
											pos2,
											logicDelrin,            //its logical volume
											"Delrin",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
	
	
	
	
	
	
	//################################################### END OF Frontal part of the Probe
	
	
	if(CaseDepth>0){
	
	
	//###################################################
	//Back Part of the Probe (Case)
	//##########################
	
	//###################################################
	// Top Probe Case
	//##########################
	
	
	G4Tubs* solidTopCase =
	new G4Tubs("TopCase",
						 0.,
						 PVC_outer_r - LateralCaseThickness,
						 TopCaseDepth*0.5,
						 SPhiTopCase,
						 DPhiTopCase);
	
	
	G4Tubs* solidAroundTopCase =
	new G4Tubs("AroundTopCase",
						 PVC_outer_r - LateralCaseThickness,
						 PVC_outer_r,
						 TopCaseDepth*0.5,
						 SPhiAroundTopCase,
						 DPhiAroundTopCase);
	
	
	G4LogicalVolume* logicTopCase =
	new G4LogicalVolume(solidTopCase,               //its solid
											TopCase_mat,           //its material
											"TopCase");            //its name
	/*
	G4LogicalVolume* logicAroundTopCase =
	new G4LogicalVolume(solidAroundTopCase,               //its solid
											AroundTopCase_mat,           //its material
											"AroundTopCase");            //its name
	*/
	
	
	G4double ProbeTopCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
	G4ThreeVector posTopCase = G4ThreeVector(fX0Scan, 0, ProbeTopCase_Posz);

	
	new G4PVPlacement(0,                     //no rotation
										posTopCase,       //at (0,0,0)
										logicTopCase,            //its logical volume
										"TopCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	/*
	new G4PVPlacement(0,                     //no rotation
										posTopCase,       //at (0,0,0)
										logicAroundTopCase,            //its logical volume
										"AroundTopCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	*/
	
	//###################################################
	// Middle Probe case
	//##########################
	
	G4double MiddleCaseDepth= (CaseDepth - HorsesShoeBackThickness - BackCaseThickness);
	
	G4Tubs* solidMiddleCase =
	new G4Tubs("MiddleCase",
						 0.,
						 PVC_outer_r - LateralCaseThickness - HorsesShoeLateralThickness,
						 MiddleCaseDepth * 0.5,
						 SPhiTopCase,
						 DPhiTopCase);
	
	
	G4Tubs* solidAroundMiddleCase =
	new G4Tubs("AroundMiddleCase",
						 PVC_outer_r - LateralCaseThickness - HorsesShoeLateralThickness,
						 PVC_outer_r - LateralCaseThickness,
						 MiddleCaseDepth * 0.5,
						 SPhiAroundTopCase,
						 DPhiAroundTopCase);
	
	G4Tubs* solidExtMiddleCase =
	new G4Tubs("AroundMiddleCase",
						 PVC_outer_r - LateralCaseThickness,
						 PVC_outer_r,
						 MiddleCaseDepth * 0.5,
						 SPhiAroundTopCase,
						 DPhiAroundTopCase);
	
	
	G4LogicalVolume* logicMiddleCase =
	new G4LogicalVolume(solidMiddleCase,               //its solid
											MiddleCase_mat,           //its material
											"MiddleCase");            //its name
	/*
	G4LogicalVolume* logicAroundMiddleCase =
	new G4LogicalVolume(solidAroundMiddleCase,               //its solid
											AroundMiddleCase_mat,           //its material
											"AroundMiddleCase");            //its name
	
	G4LogicalVolume* logicExtMiddleCase =
	new G4LogicalVolume(solidExtMiddleCase,               //its solid
											ExtMiddleCase_mat,           //its material
											"ExtMiddleCase");            //its name
	
	*/
	
	G4double ProbeMiddleCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth * 0.5;
	G4ThreeVector posMiddleCase = G4ThreeVector(fX0Scan, 0, ProbeMiddleCase_Posz);
	
	
	new G4PVPlacement(0,                     //no rotation
										posMiddleCase,       //at (0,0,0)
										logicMiddleCase,            //its logical volume
										"MiddleCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	/*
	new G4PVPlacement(0,                     //no rotation
										posMiddleCase,       //at (0,0,0)
										logicAroundMiddleCase,            //its logical volume
										"AroundMiddleCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	new G4PVPlacement(0,                     //no rotation
										posMiddleCase,       //at (0,0,0)
										logicExtMiddleCase,            //its logical volume
										"ExtMiddleCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	*/
	
	
	//###################################################
	// Back Middle Probe case
	//##########################
	
	
	
	G4Tubs* solidBackMiddleCase =
	new G4Tubs("BackMiddleCase",
						 0.,
						 PVC_outer_r - LateralCaseThickness,
						 HorsesShoeBackThickness * 0.5,
						 SPhiTopCase,
						 DPhiTopCase);
	
	
	G4Tubs* solidAroundBackMiddleCase =
	new G4Tubs("AroundBackMiddleCase",
						 PVC_outer_r - LateralCaseThickness,
						 PVC_outer_r,
						 HorsesShoeBackThickness * 0.5,
						 SPhiAroundTopCase,
						 DPhiAroundTopCase);
	
	
	/*
	G4LogicalVolume* logicBackMiddleCase =
	new G4LogicalVolume(solidBackMiddleCase,               //its solid
											BackMiddleCase_mat,           //its material
											"BackMiddleCase");            //its name
	
	G4LogicalVolume* logicAroundBackMiddleCase =
	new G4LogicalVolume(solidAroundBackMiddleCase,               //its solid
											AroundBackMiddleCase_mat,           //its material
											"ArounBackdMiddleCase");            //its name
	
	
	
	G4double ProbeBackMiddleCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth + HorsesShoeBackThickness * 0.5;
	G4ThreeVector posBackMiddleCase = G4ThreeVector(fX0Scan, 0, ProbeBackMiddleCase_Posz);
	
	
	new G4PVPlacement(0,                     //no rotation
										posBackMiddleCase,       //at (0,0,0)
										logicBackMiddleCase,            //its logical volume
										"BackMiddleCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	new G4PVPlacement(0,                     //no rotation
										posBackMiddleCase,       //at (0,0,0)
										logicAroundBackMiddleCase,            //its logical volume
										"AroundBackMiddleCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	

	
	*/
	
	//###################################################
	//  End Probe case
	//##########################
	
	
	
	G4Tubs* solidEndCase =
	new G4Tubs("EndCase",
						 0.,
						 PVC_outer_r,
						 BackCaseThickness * 0.5,
						 SPhiTopCase,
						 DPhiTopCase);
	
	
	/*
	
	G4LogicalVolume* logicEndCase =
	new G4LogicalVolume(solidEndCase,               //its solid
											EndCase_mat,           //its material
											"EndCase");            //its name
	

	
	
	G4double ProbeEndCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth  + HorsesShoeBackThickness + BackCaseThickness * 0.5;
	G4ThreeVector posEndCase = G4ThreeVector(fX0Scan, 0, ProbeEndCase_Posz);
	
	
	new G4PVPlacement(0,                     //no rotation
										posEndCase,       //at (0,0,0)
										logicEndCase,            //its logical volume
										"EndCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	*/


	//###################################################
	// G4Union Probe
	//##########################
	
	
	
	//###################################################
	// G4Union External Case
	//##########################
	
	
	G4VSolid* Union1
	= new G4UnionSolid("Union1", solidAroundTopCase , solidExtMiddleCase,0 , G4ThreeVector(0.,0.,(MiddleCaseDepth+TopCaseDepth)*0.5));
	//The G4ThreeVector(0.,0.,(MiddleCaseDepth+TopCaseDepth)*0.5) generates a translation of the second solid respect to the center of the first one. When I place the resulting solid I've to take in consideration as center of it the center of the first solid.
	
	
	
	
	G4VSolid* Union2
	= new G4UnionSolid("Union2", solidAroundBackMiddleCase , solidEndCase, 0 , G4ThreeVector(0.,0.,(BackCaseThickness + HorsesShoeBackThickness)*0.5));
	
	
	G4VSolid* PlasticCase
	= new G4UnionSolid("PlasticCase", Union1 , Union2 ,0, G4ThreeVector(0.,0.,MiddleCaseDepth + TopCaseDepth*0.5 + HorsesShoeBackThickness*0.5));
	
	G4LogicalVolume* logicPlasticCase =
	new G4LogicalVolume(PlasticCase,               //its solid
											PlasticCase_mat,           //its material
											"PlasticCase");            //its name
	
	
	
	G4double PlastiCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
	G4ThreeVector posPlasticCase = G4ThreeVector(fX0Scan, 0, PlastiCase_Posz);
	
	
	
	
	new G4PVPlacement(0,                     //no rotation
										posPlasticCase,        //at (0,0,0)
										logicPlasticCase,            //its logical volume
										"PlasticCase",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	//###################################################
	// G4Union HorsesShoe Probe
	//##########################
	
	
	G4VSolid* HorsesShoe
	= new G4UnionSolid("HorsesShoe", solidAroundMiddleCase , solidBackMiddleCase  ,0, G4ThreeVector(0.,0.,(MiddleCaseDepth +  HorsesShoeBackThickness)*0.5));
	
	
	G4LogicalVolume* logicHorsesShoe =
	new G4LogicalVolume(HorsesShoe,               //its solid
											HorsesShoe_mat,           //its material
											"HorsesShoe");            //its name
	
	G4double HorsesShoe_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth*0.5;
	G4ThreeVector posHorsesShoe = G4ThreeVector(fX0Scan, 0, HorsesShoe_Posz);
										 
	
	new G4PVPlacement(0,                     //no rotation
							      posHorsesShoe,        //at (0,0,0)
							      logicHorsesShoe,            //its logical volume
							      "HorsesShoe",               //its name
								    logicWorld,            //its mother  volume
								    false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
										 
	
	
	
	//###################################################
	// End of Probe
	//##########################
	
	}
	
	
	
	return physWorld;
	
	
	
	
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Voglio Pter_inner_r da terminale


