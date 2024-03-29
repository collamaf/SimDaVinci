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
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction(G4double x0, G4double ZValue, G4double AbsHoleDiam, G4int SourceSelect, G4int AbsorberMaterial, G4double PterDiameter, G4double PterThickness, G4double PVCDiameter, G4double PVCLateralThickness, G4double SourceDiameter, G4double SourceThickness, G4double AbsorberThickness, G4double ProbeCaseDepth, G4double ProbeCaseLateralThickness, G4double ProbeCaseBackThickness, G4double HSLateralThickness, G4double HSBackThickness, G4int HousingCase, G4bool ScintFlag, G4int GaSet, G4int ApparatusMat, G4bool SecondShieldFlag)
	: G4VUserDetectorConstruction(),
	  fScoringVolume(0), fX0Scan(x0 * mm), fZValue(ZValue), fAbsHoleDiam(AbsHoleDiam), fSourceSelect(SourceSelect), fAbsorberMaterial(AbsorberMaterial), fPterDiameter(PterDiameter), fPterThickness(PterThickness), fPVCDiameter(PVCDiameter), fPVCLateralThickness(PVCLateralThickness), fSourceDiameter(SourceDiameter), fSourceThickness(SourceThickness), fAbsorberThickness(AbsorberThickness), fCaseDepth(ProbeCaseDepth), fLateralCaseThickness(ProbeCaseLateralThickness), fBackCaseThickness(ProbeCaseBackThickness), fHorsesShoeLateralThickness(HSLateralThickness), fHorsesShoeBackThickness(HSBackThickness), fHousingCase(HousingCase), fScintFlag(ScintFlag), fGaSet(GaSet), fApparatusMat(ApparatusMat), fSecondShieldFlag(SecondShieldFlag)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B1DetectorConstruction::Construct()
{
	// Get nist material manager
	G4NistManager *nist = G4NistManager::Instance();

	// Option to switch on/off checking of volumes overlaps
	//
	G4bool checkOverlaps = true;
	G4bool secondShieldFlag = fSecondShieldFlag;
	//
	// World
	//
	G4double world_sizeXY = 1. * m;
	G4double world_sizeZ = 1. * m;
	G4Material *world_mat = nist->FindOrBuildMaterial("G4_AIR");
	//	G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

	G4Box *solidWorld =
		new G4Box("World",													  // its name
				  0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ); // its size

	G4LogicalVolume *logicWorld =
		new G4LogicalVolume(solidWorld, // its solid
							world_mat,	// its material
							"World");	// its name

	G4VPhysicalVolume *physWorld =
		new G4PVPlacement(0,			   // no rotation
						  G4ThreeVector(), // at (0,0,0)
						  logicWorld,	   // its logical volume
						  "World",		   // its name
						  0,			   // its mother  volume
						  false,		   // no boolean operation
						  0,			   // copy number
						  checkOverlaps);  // overlaps checking

#pragma mark Definition of Materials
	// ###################################################################
	// ###################################################
	//  Definitions of materials
	// ##########################

	G4double z, a;
	G4String name, symbol;
	G4int ncomponents, natoms;

	a = 1.01 * g / mole;
	G4Element *elH = new G4Element(name = "Hydrogen", symbol = "H", z = 1., a);
	a = 12.01 * g / mole;
	G4Element *elC = new G4Element(name = "Carbon", symbol = "C", z = 6., a);
	a = 16.00 * g / mole;
	G4Element *elO = new G4Element(name = "Oxygen", symbol = "O", z = 8., a);
	a = 14.00 * g / mole;
	G4Element *elN = new G4Element(name = "Nitrogen", symbol = "N", z = 7., a);

	G4double densityAlu = 2.600 * g / cm3;
	G4NistManager::Instance()->BuildMaterialWithNewDensity("MyAlu", "G4_Al", densityAlu);

	// ###################################################
	//  AGAR AGAR Source - AgarAgar should be C14 H24 O9
	// ##########################
	G4double Agardensity = 1.030 * g / cm3;
	G4Material *AgarAgar = new G4Material(name = "AgarAgar", Agardensity, ncomponents = 3);
	AgarAgar->AddElement(elH, natoms = 24);
	AgarAgar->AddElement(elC, natoms = 14);
	AgarAgar->AddElement(elO, natoms = 9);

	// ###################################################
	//  ABS material - ABS should be C15 H17 N
	// ##########################
	//	G4double ABSdensity = 0.9*g/cm3;
	G4double ABSdensity = 1.037 * g / cm3;
	//	G4double ABSdensity = 1.08*g/cm3;
	G4Material *ABS = new G4Material(name = "ABS", ABSdensity, ncomponents = 3);
	ABS->AddElement(elH, natoms = 17);
	ABS->AddElement(elC, natoms = 15);
	ABS->AddElement(elN, natoms = 1);

	// ###################################################
	//  P-Terphenyl Material
	// ##########################
	G4double PTerphenyldensity = 1.23 * g / cm3;
	G4Material *PTerphenyl = new G4Material(name = "PTerphenyl", PTerphenyldensity, ncomponents = 2);
	PTerphenyl->AddElement(elC, natoms = 18);
	PTerphenyl->AddElement(elH, natoms = 14);

	// ###################################################
	//  Delrin Material      C H2 O
	// ##########################
	G4double Delrindensity = 1.41 * g / cm3;
	G4Material *Delrin = new G4Material(name = "Delrin", Delrindensity, ncomponents = 3);
	Delrin->AddElement(elC, natoms = 1);
	Delrin->AddElement(elH, natoms = 2);
	Delrin->AddElement(elO, natoms = 1);

	// ###################################################
	//  GaContainer2 Material (ABS)
	// ##########################
	G4double GaContainer2Matdensity = 0.43 * g / cm3; // mean measured density
	G4Material *GaContainer2Mat = new G4Material(name = "GaContainer2Mat", GaContainer2Matdensity, ncomponents = 3);
	GaContainer2Mat->AddElement(elH, natoms = 17);
	GaContainer2Mat->AddElement(elC, natoms = 15);
	GaContainer2Mat->AddElement(elN, natoms = 1);

	// ###################################################
	//  PEEK Material
	//  C 114/150 = 76
	//  H 12/150  = 8
	//  O 24/150  = 16
	// ##########################
	G4Material *Peek = new G4Material("Peek", 1.31 * g / cm3, ncomponents = 3);
	Peek->AddMaterial(nist->FindOrBuildMaterial("G4_C"), 76 * perCent);
	Peek->AddMaterial(nist->FindOrBuildMaterial("G4_H"), 8 * perCent);
	Peek->AddMaterial(nist->FindOrBuildMaterial("G4_O"), 16 * perCent);

	// ###################################################
	//  KAPTON Material
	G4Material *Kapton = new G4Material(name = "Kapton", 1.42 * g / cm3, ncomponents = 4);
	Kapton->AddElement(elH, 0.0273);
	Kapton->AddElement(elC, 0.7213);
	Kapton->AddElement(elN, 0.0765);
	Kapton->AddElement(elO, 0.1749);

	// ##########################
	// ###################################################
	// ###################################################################

#pragma mark Assignment of Materials
	// ###################################################################
	// ###################################################
	//  Assignment of materials
	// ##########################

	G4Material *SourceExtY_mat = AgarAgar;
	G4Material *ABSaround_mat = ABS;
	G4Material *ABSbehind_mat = ABS;
	G4Material *GaContainer_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material *GaContainer2_mat = GaContainer2Mat;
	G4Material *GaContainer3_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material *ProbeContainer_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material *SourceExtGa_mat = nist->FindOrBuildMaterial("G4_WATER");
	G4Material *SourceSR_mat = nist->FindOrBuildMaterial("MyAlu");
	G4Material *SourceGenericIon_mat = nist->FindOrBuildMaterial("G4_WATER");
	G4Material *FrontShield_mat = nist->FindOrBuildMaterial("MyAlu");
	G4Material *shapeDummy_mat = nist->FindOrBuildMaterial("G4_AIR");
	G4Material *Pter_mat = PTerphenyl;
	G4Material *PVC_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material *Delrin_mat = Delrin;
	G4Material *Absorber_mat = nist->FindOrBuildMaterial("G4_Cu");
	G4Material *SiPM_mat = nist->FindOrBuildMaterial("G4_Si");
	G4Material *Table_mat = nist->FindOrBuildMaterial("MyAlu");

	G4NistManager::Instance()->BuildMaterialWithNewDensity("BlackAbsMaterial", "G4_POLYVINYL_CHLORIDE", 1.4 * g / cm3);
	G4NistManager::Instance()->BuildMaterialWithNewDensity("NUCLEOMEDMaterial", "G4_POLYVINYL_CHLORIDE", 1.21 * g / cm3);

	// Materials for NL probe case
	G4Material *PlasticCase_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material *HorsesShoe_mat = nist->FindOrBuildMaterial("G4_Pb");
	G4Material *CaseInner_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	if (fHousingCase == 2) // 1- InnerPlastic: Pb+Pl || 2- InnerAir: Pb+Air || 3- TotalAir: Air+Air (norm)
		CaseInner_mat = world_mat;
	else if (fHousingCase == 3)
	{ // main.cc default
		HorsesShoe_mat = world_mat;
		CaseInner_mat = world_mat;
	}
	G4Material *MiddleCase_mat = CaseInner_mat;
	G4Material *TopCase_mat = CaseInner_mat;

	if (fSourceSelect == 6 || fSourceSelect == 7 || fSourceSelect == 10)
		SourceSR_mat = world_mat;

	// To modify GaSet 2/3 materials according to flags
	if (fGaSet == 3)
	{
		if (fApparatusMat == 1)
		{
			ProbeContainer_mat = world_mat;
			GaContainer2_mat = world_mat;
			GaContainer3_mat = world_mat;
		}
		else if (fApparatusMat == 2)
		{
			ProbeContainer_mat = HorsesShoe_mat;
			GaContainer2_mat = HorsesShoe_mat;
			GaContainer3_mat = HorsesShoe_mat;
		}
	}

	// ###################################################################
	// ###################################################
	//  Optics characteristics part
	// ##########################
	//
	//  ------------ Generate & Add Material Properties Table ------------
	//
	G4double photonEnergy[] =
		{2.96 * eV};

	const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

	G4double refractiveIndex1[] =
		{1.65}; // da misteriosa mail del 30.9.2013, confermato da silvio il 07-01-2019

	assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

	G4double absorption[] =
		{20 * mm}; // da elsarticle CMT, è lunghezza di attenuazione "fisica" (non "efficace")

	assert(sizeof(absorption) == sizeof(photonEnergy));

	G4double scintilFast[] =
		{1.00};

	assert(sizeof(scintilFast) == sizeof(photonEnergy));

	G4double scintilSlow[] =
		{0.0};

	assert(sizeof(scintilSlow) == sizeof(photonEnergy));

	G4MaterialPropertiesTable *materialTablePter = new G4MaterialPropertiesTable();
	//	materialTablePter->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
	//	->SetSpline(true);
	//	materialTablePter->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
	//	->SetSpline(true);
	//	materialTablePter->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
	//	->SetSpline(true);
	//	materialTablePter->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
	//	->SetSpline(true);

	//	materialTablePter->AddConstProperty("SCINTILLATIONYIELD",0.1*28000./MeV); //33k da nostro papero, 28k da papero recente elsa CMT
	//	materialTablePter->AddConstProperty("RESOLUTIONSCALE",1.0);
	////	materialTablePter->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
	//	materialTablePter->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
	//	materialTablePter->AddConstProperty("YIELDRATIO",1);
	//
	//
	//	G4double refractiveIndexDelrin = 1.48;
	//	G4MaterialPropertiesTable* materialTableDelrin = new G4MaterialPropertiesTable();
	//	materialTableDelrin->AddConstProperty("RINDEX", refractiveIndexDelrin);
	//	materialTableDelrin->AddConstProperty("REFLECTIVITY",1);
	//
	//
	//	G4cout << "PTERP G4MaterialPropertiesTable" << G4endl;
	//	materialTablePter->DumpTable();
	//	materialTableDelrin->DumpTable();
	//

	/*
	 G4OpticalSurface* wrapper = new G4OpticalSurface("wrapper");
	 new G4LogicalBorderSurface("wrapper", slab, expHall_phys, wrapper);
	 wrapper->SetType(dielectric_metal);
	 wrapper->SetFinish(polished);
	 wrapper->SetModel(glisur);
	 const G4int NUM = 2;
	 G4double pp[NUM] = {2.0*eV, 3.5*eV};
	 G4double reflectivity[NUM] = {1., 1.};
	 G4double efficiency[NUM] = {0.0, 0.0};
	 G4MaterialPropertiesTable* wrapperProperty = new G4MaterialPropertiesTable();
	 wrapperProperty->AddProperty("REFLECTIVITY", pp, reflectivity, NUM);
	 wrapperProperty->AddProperty("EFFICIENCY", pp, efficiency, NUM);
	 wrapper->SetMaterialPropertiesTable(wrapperProperty);


	 */

	//	if (fScintFlag) {
	//		PTerphenyl->SetMaterialPropertiesTable(materialTablePter); //to toggle scintillation
	//		Delrin->SetMaterialPropertiesTable(materialTableDelrin);
	//	}
	// ##########################
	// ###################################################
	// ###################################################################

#pragma mark Definitions of Dimensions and Sizes

	// ###################################################################
	// ###################################################
	//  Definitions of dimensions and sizes
	// ##########################

	G4double Ang0 = 0. * deg;
	G4double Ang2Pi = 360. * deg;

	// ### ExtY SOURCE
	G4double RminSourceExtY = 0. * mm;
	G4double RmaxSourceExtY = 10.5 * mm;
	G4double DzSourceExtY = 4.5 * mm;
	// ###

	// ### ABS structure for Agar sources
	G4double RminABSaround = RmaxSourceExtY;
	G4double RmaxABSaround = 12.5 * mm;
	G4double DzABSaround = DzSourceExtY;

	G4double RminABSbehind = 0. * mm;
	G4double RmaxABSbehind = RmaxABSaround;
	G4double DzABSbehind = 3 * mm;
	// ###

	// ### Sr Source
	G4double RminSourceSR = 0. * mm;
	G4double RmaxSourceSR = 12.5 * mm;
	G4double DzSourceSR = 3 * mm;
	// ###

	// ### Cu Source
	G4double RminSourceGenericIon = 0. * mm;
	G4double RmaxSourceGenericIon = fSourceDiameter * 0.5 * mm;
	G4double DzSourceGenericIon = fSourceThickness * mm;
	// ###

	// ### Ga Source Container
	G4double rSourceExtGa = fSourceDiameter * 0.5 * mm;
	G4double RContainer = 50 * mm;
	G4double dzSourceExtGa = fSourceThickness * mm;
	G4double DzContainer = 20 * mm;
	// ###

	// ### FrontShield (Al)
	G4double Z_FrontShield = 0 * mm;
	G4double Z_FrontShieldBis = 0 * mm;
	G4double FrontShield_outer_r = (12.0 / 2.0) * mm;
	G4double FrontShield_sizeZ = 15 * um;
	G4double FrontShieldBis_sizeZ = !secondShieldFlag ? 0.0 * mm : 10 * um;
	G4ThreeVector posFrontShield;
	G4ThreeVector posFrontShieldBis;
	// ###

	// ###

	// ### SiPm
	G4double DxSiPm = 3. * mm;
	G4double DySiPm = 3. * mm;
	G4double DzSiPm = 300.e-6 * mm;

	// ### Table
	G4double DxTable = 500 * mm;
	G4double DyTable = 500 * mm;
	G4double DzTable = 50 * mm;

	// ### Pter
	G4double Pter_Diam = fPterDiameter * mm;
	G4double PVC_outer_r = fPVCDiameter * mm / 2.;
	G4double Pter_sizeZ = fPterThickness * mm;
	G4double Pter_Posz = 0. * mm;
	G4double Pter_ZScan = fZValue * mm;
	G4double PVC_inner_r = PVC_outer_r - fPVCLateralThickness * mm;

	// ### Absorber
	G4double RminAbs = fabs(fAbsHoleDiam) / 2. * mm;
	G4double RmaxAbs = 30 * mm;
	G4double DzAbs = fAbsorberThickness * mm;
	if (fGaSet == 3)
	{
		RmaxAbs = 22 / 2. * mm;
	};
	if (fCaseDepth > 0)
	{
		RmaxAbs = PVC_outer_r;
	};

	// ### Probe Case
	G4double CaseDepth = fCaseDepth * mm;
	G4double LateralCaseThickness = fLateralCaseThickness * mm;
	G4double BackCaseThickness = fBackCaseThickness * mm;
	// G4double ProbeCase_Posz=0.*mm;
	G4double TopCaseDepth = 2. * mm;
	// G4double Case_r = PVC_outer_r +0.1*mm;
	G4double HorsesShoeLateralThickness = fHorsesShoeLateralThickness * mm;
	G4double HorsesShoeBackThickness = fHorsesShoeBackThickness * mm;

	// ### GaSet 3
	G4double D_CylABCD3 = 50.8 * mm;
	G4double H_CylA3 = 34.7 * mm;
	G4double d_CylB3 = 10 * mm;
	G4double H_CylB3 = 7.5 * mm;
	G4double d_CylC3 = 22 * mm;
	G4double H_CylC3 = 3.3 * mm;
	G4double d_CylD3 = 41.3 * mm;
	G4double H_CylD3 = 13.8 * mm;

	G4double D_CylEF3 = 41.3 * mm;
	G4double d_CylE3 = 26 * mm;
	G4double H_CylE3 = 5.45 * mm;
	G4double d_CylF3 = 13 * mm;
	G4double H_CylF3 = 8 * mm;
	G4double D_CylG3 = 51 * mm;
	G4double d_CylG3 = 13 * mm;
	G4double H_CylG3 = 32 * mm;
	ContainerOuterRadius = D_CylG3 / 2.;
	ContainerMaxHeigth = H_CylA3;

	// ### Dummy Exit Sorg
	G4double RminDummyExitSorg = 0. * mm;
	G4double RmaxDummyExitSorg = 18. * mm;
	if (fGaSet == 3)
		RmaxDummyExitSorg = d_CylG3 / 2. * mm;
	G4double DzDummyExitSorg = 1.e-5 * mm;
	G4ThreeVector posDummyExitSorg;
	// ###

	// ### Dummy Exit Abs
	G4double RminDummyExitAbs = 0. * mm;
	G4double RmaxDummyExitAbs = 18. * mm;
	if (fGaSet == 3)
		RmaxDummyExitAbs = RmaxAbs;
	G4double DzDummyExitAbs = 1.e-5 * mm;
	G4ThreeVector posDummyExitAbs;
	// ###

	// ### Dummy Enter Pter
	G4double RminDummyEnterProbe = 0. * mm;
	G4double RmaxDummyEnterProbe = 18. * mm;
	//	G4double RmaxDummyEnterProbe = PVC_outer_r;
	if (fGaSet == 3)
		RmaxDummyEnterProbe = d_CylG3 / 2. * mm;
	G4double DzDummyEnterProbe = 1.e-5 * mm;
	G4ThreeVector posDummyEnterProbe;
	// ###

	G4ThreeVector posAbs;

	G4double FExtLiquidLateralHeigth = fSourceThickness * mm;
	G4double FExtLiquidBottomThickness = 1 * cm;
	// ##########################
	// ###################################################

#pragma mark Definitions of Volumes
	// ###################################################################
	// ###################################################
	//  Definitions of volumes
	// ##########################
	// ###################################################################

	// ###########################
	// Dummy exit sorg
	// ###########################
	G4Tubs *solidShapeDummyExitSorg =
		new G4Tubs("DummyExitSorg", // its name
				   RminDummyExitSorg,
				   RmaxDummyExitSorg,
				   0.5 * DzDummyExitSorg,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicShapeDummyExitSorg =
		new G4LogicalVolume(solidShapeDummyExitSorg, // its solid
							shapeDummy_mat,			 // its material
							"DummyExitSorg");		 // its name

	// ###########################
	// Dummy Enter Abs
	// ###########################
	G4Tubs *solidShapeDummyExitAbs =
		new G4Tubs("DummyExitAbs", // its name
				   RminDummyExitAbs,
				   RmaxDummyExitAbs,
				   0.5 * DzDummyExitAbs,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicShapeDummyExitAbs =
		new G4LogicalVolume(solidShapeDummyExitAbs, // its solid
							shapeDummy_mat,			// its material
							"DummyExitAbs");		// its name

	// ###########################
	// Dummy Enter Pter
	// ###########################
	G4Tubs *solidShapeDummyEnterProbe =
		new G4Tubs("DummyEnterProbe", // its name
				   RminDummyEnterProbe,
				   RmaxDummyEnterProbe,
				   0.5 * DzDummyEnterProbe,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicShapeDummyEnterProbe =
		new G4LogicalVolume(solidShapeDummyEnterProbe, // its solid
							shapeDummy_mat,			   // its material
							"DummyEnterProbe");		   // its name

	// ###################################################
	//  ExtY Source
	// ##########################
	G4ThreeVector posSourceExtY = G4ThreeVector(0, 0, -DzSourceExtY * 0.5 - DzDummyExitSorg);

	G4Tubs *solidSourceExtY =
		new G4Tubs("Source", // its name
				   RminSourceExtY,
				   RmaxSourceExtY,
				   0.5 * DzSourceExtY,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicSourceExtY =
		new G4LogicalVolume(solidSourceExtY, // its solid
							SourceExtY_mat,	 // its material
							"Source");		 // its name

	// ################################################### END ExtY SOURCE

	// ###################################################
	//  ABS carrier around ExtY Source
	// ##########################
	G4ThreeVector posABSaround = G4ThreeVector(0, 0, -DzABSaround * 0.5 - DzDummyExitSorg);

	G4Tubs *solidABSaround =
		new G4Tubs("ABSaround", // its name
				   RminABSaround,
				   RmaxABSaround,
				   0.5 * DzABSaround,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicABSaround =
		new G4LogicalVolume(solidABSaround, // its solid
							ABSaround_mat,	// its material
							"ABSaround");	// its name

	// ################################################### END ABS AROUND

	// ###################################################
	//  ABS carrier behind ExtY Source
	// ##########################
	G4ThreeVector posABSbehind = G4ThreeVector(0, 0, -DzABSbehind * 0.5 - DzABSaround - DzDummyExitSorg);

	G4Tubs *solidABSbehind =
		new G4Tubs("ABSbehind", // its name
				   RminABSbehind,
				   RmaxABSbehind,
				   0.5 * DzABSbehind,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicABSbehind =
		new G4LogicalVolume(solidABSbehind, // its solid
							ABSbehind_mat,	// its material
							"ABSbehind");	// its name

	// ################################################### END ABS BEHIND

	// ###################################################
	//  Sr90 lab Source
	// ##########################
	G4ThreeVector posSourceSR = G4ThreeVector(0, 0, -DzSourceSR * 0.5 - DzDummyExitSorg);

	G4cout << "GEOMETRY DEBUG - Z thickness of solidSourceSR= " << DzSourceSR / mm << ", Z pos= " << posSourceSR.z() / mm << G4endl;

	G4Tubs *solidSourceSR =
		new G4Tubs("Source", // its name
				   RminSourceSR,
				   RmaxSourceSR,
				   0.5 * DzSourceSR,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicSourceSR =
		new G4LogicalVolume(solidSourceSR, // its solid
							SourceSR_mat,  // its material
							"Source");	   // its name

	// ################################################### END SR SOURCE

	// ###################################################
	//  Ga68 puntforme
	// ##########################
	G4double DzSourceGaPoint = 3 * cm;
	G4ThreeVector posSourceGaPoint = G4ThreeVector(0, 0, -DzSourceGaPoint * 0.5 - DzDummyExitSorg);
	G4Tubs *solidSourceGaPoint =
		new G4Tubs("Source", // its name
				   RminSourceSR,
				   RmaxSourceSR,
				   0.5 * DzSourceGaPoint,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicSourceGaPoint =
		new G4LogicalVolume(solidSourceGaPoint, // its solid
							ABSbehind_mat,		// its material
							"Source");			// its name

	// ################################################### END SR SOURCE

	// ###################################################
	//  Cu67/64-F18 volume Sources for pure MC tests (not exp meas.)
	// ##########################
	G4ThreeVector posSourceGenericIon = G4ThreeVector(0, 0, -DzSourceGenericIon * 0.5 - DzDummyExitSorg);

	G4Tubs *solidSourceGenericIon =
		new G4Tubs("Source", // its name
				   RminSourceGenericIon,
				   RmaxSourceGenericIon,
				   0.5 * DzSourceGenericIon,
				   Ang0,
				   Ang2Pi); // its size

	G4LogicalVolume *logicSourceGenericIon =
		new G4LogicalVolume(solidSourceGenericIon, // its solid
							SourceGenericIon_mat,  // its material
							"Source");			   // its name

	// ################################################### END SR SOURCE

	// ###################################################################################
	// ###################################################
	//  Probe's Solids
	// ##########################

	// ###########################
	// Front Shield
	// ###########################
	G4Tubs *solidFrontShield =
		new G4Tubs("FrontShield",													// its name
				   0., FrontShield_outer_r, FrontShield_sizeZ * 0.5, Ang0, Ang2Pi); // its size

	G4LogicalVolume *logicFrontShield =
		new G4LogicalVolume(solidFrontShield, // its solid
							FrontShield_mat,  // its material
							//											world_mat,           //its material
							"FrontShield"); // its name

	G4VPhysicalVolume *physFrontShield;

	// ###########################
	// Front Shield Bis (scotch to cover hole) - addedd on 2021.01.21 by collamaf
	// ###########################
	G4double FrontShieldBis_sizeZBIS = fSecondShieldFlag ? FrontShieldBis_sizeZ : 1 * mm;
	G4Tubs *solidFrontShieldBis =
		new G4Tubs("FrontShieldBis",													  // its name
				   0., FrontShield_outer_r, FrontShieldBis_sizeZBIS * 0.5, Ang0, Ang2Pi); // its size

	G4LogicalVolume *logicFrontShieldBis =
		new G4LogicalVolume(solidFrontShieldBis, // its solid
							Kapton,				 // its material
							"FrontShieldBis");	 // its name

	G4VPhysicalVolume *physFrontShieldBis;

	// ################ Pter
	G4Tubs *solidPter =
		new G4Tubs("Pter", // its name
				   0.,
				   Pter_Diam * 0.5,
				   Pter_sizeZ * 0.5,
				   Ang0, Ang2Pi); // its size

	G4LogicalVolume *logicPter =
		new G4LogicalVolume(solidPter, // its solid
							Pter_mat,  // its material
							"Pter");   // its name

	G4VPhysicalVolume *physPter;

	// ################ SiPM
	G4Box *solidSiPm =
		new G4Box("SiPm", // its name
				  DxSiPm / 2.,
				  DySiPm / 2.,
				  DzSiPm / 2.); // its size

	G4LogicalVolume *logicSiPm =
		new G4LogicalVolume(solidSiPm, // its solid
							SiPM_mat,  // its material
							"SiPm");   // its name

	G4VPhysicalVolume *physSiPm;

	// ################ Table
	G4Box *solidTable =
		new G4Box("Table", // its name
				  DxTable / 2.,
				  DyTable / 2.,
				  DzTable / 2.); // its size

	G4LogicalVolume *logicTable =
		new G4LogicalVolume(solidTable, // its solid
							Table_mat,	// its material
							"Table");	// its name

	// ################ PVC Around delrin
	G4Tubs *solidPVC =
		new G4Tubs("PVC",													  // its name
				   PVC_inner_r, PVC_outer_r, Pter_sizeZ * 0.5, Ang0, Ang2Pi); // its size

	G4LogicalVolume *logicPVC =
		new G4LogicalVolume(solidPVC, // its solid
							PVC_mat,  // its material
							"PVC");	  // its name

	G4VPhysicalVolume *physPVC;

	// ################ Delrin Around Pter
	G4Tubs *solidDelrin =
		new G4Tubs("Delrin", Pter_Diam * 0.5, PVC_inner_r, Pter_sizeZ * 0.5, Ang0, Ang2Pi); // its size

	G4LogicalVolume *logicDelrin =
		new G4LogicalVolume(solidDelrin, // its solid
							Delrin_mat,	 // its material
							"Delrin");	 // its name

	G4VPhysicalVolume *physDelrin;

	// ################ Absorber
	G4Tubs *solidAbsorber = new G4Tubs("Absorber",
									   RminAbs,
									   RmaxAbs,
									   0.5 * DzAbs,
									   Ang0,
									   Ang2Pi); // its size

	if (fAbsorberMaterial == 1)
	{
		Absorber_mat = nist->FindOrBuildMaterial("G4_Cu");
	}
	else if (fAbsorberMaterial == 2)
	{
		Absorber_mat = nist->FindOrBuildMaterial("G4_Pb");
	}
	else if (fAbsorberMaterial == 3)
	{
		Absorber_mat = nist->FindOrBuildMaterial("MyAlu");
	}
	else if (fAbsorberMaterial == 4)
	{
		Absorber_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	}
	else if (fAbsorberMaterial == 5)
	{
		Absorber_mat = ABS;
	}
	else if (fAbsorberMaterial == 6)
	{ // Assorbitori neri, probabilmente PVC ma densità 1.4
		Absorber_mat = nist->FindOrBuildMaterial("BlackAbsMaterial");
	}
	else if (fAbsorberMaterial == 7)
	{ // Material stampa 3d NUCLEOMED densità 1.4
		Absorber_mat = nist->FindOrBuildMaterial("NUCLEOMEDMaterial");
	}
	else if (fAbsorberMaterial == 8)
	{ // PEEK
		Absorber_mat = Peek;
	}

	G4LogicalVolume *logicAbsorber =
		new G4LogicalVolume(solidAbsorber, // its solid
							Absorber_mat,  // its material
							"Absorber");   // its name

	G4VPhysicalVolume *physAbsorber;

	// ########## New case (Ga Container made with PVC)

	// ###################################################
	//  "CATAFALCO" - Ga Source Container
	// ##########################

	G4ThreeVector posContainerExtGa = G4ThreeVector(0, 0, -H_CylA3 * 0.5 - DzDummyExitSorg);
	//		G4ThreeVector posExtGa3 = G4ThreeVector(0, 0, -H_CylB3*0.5-DzDummy3/2);
	G4ThreeVector posExtGa3 = G4ThreeVector(0, 0, -H_CylB3 + dzSourceExtGa * 0.5 - DzDummyExitSorg);

	G4VSolid *CylinderA =
		new G4Tubs("CylinderA", // its name
				   0.,
				   D_CylABCD3 * 0.5,
				   H_CylA3 * 0.5,
				   Ang0,
				   Ang2Pi); // its size

	G4VSolid *CylinderB =
		new G4Tubs("CylinderB", // its name
				   0.,
				   d_CylB3 * 0.5,
				   H_CylB3 * 0.5,
				   Ang0,
				   Ang2Pi); // its size

	G4VSolid *SourceGa =
		new G4Tubs("SourceGa",
				   0.,
				   d_CylB3 * 0.5,
				   (fSourceThickness * mm) * 0.5,
				   Ang0,
				   Ang2Pi);

	G4VSolid *BackGaContainer =
		new G4SubtractionSolid("BackGaContainer", // GaContainer=GaCylinder-GaSource
							   CylinderA,
							   CylinderB,
							   0,
							   G4ThreeVector(0., 0., H_CylA3 * 0.5 - H_CylB3 * 0.5));

	G4LogicalVolume *logicSourceExtGa2 =
		new G4LogicalVolume(SourceGa,		 // its solid
							SourceExtGa_mat, // its material
							"Source");		 // its name

	//	if (fSourceSelect!=11) new G4PVPlacement(0,                     //no rotation
	//																					 posExtGa3,             //at (0,0,0)
	//																					 logicSourceExtGa2,            //its logical volume
	//																					 "Source",               //its name
	//																					 logicWorld,            //its mother  volume
	//																					 false,                 //no boolean operation
	//																					 0,                     //copy number
	//																					 checkOverlaps);        //overlaps checking

	G4VSolid *CylinderC =
		new G4Tubs("CylinderC", // its name
				   d_CylC3 * 0.5,
				   D_CylABCD3 * 0.5,
				   H_CylC3 * 0.5,
				   Ang0,
				   Ang2Pi); // its size

	G4VSolid *CylinderD =
		new G4Tubs("CylinderD", // its name
				   d_CylD3 * 0.5,
				   D_CylABCD3 * 0.5,
				   H_CylD3 * 0.5,
				   Ang0,
				   Ang2Pi); // its size

	G4VSolid *UnionCD = new G4UnionSolid("UnionCD",
										 CylinderC,
										 CylinderD,
										 0,
										 G4ThreeVector(0., 0., (H_CylC3 + H_CylD3) * 0.5));
	// The G4ThreeVector(0.,0.,(MiddleCaseDepth+TopCaseDepth)*0.5) generates a translation of the second solid respect to the center of the first one. When I place the resulting solid I've to take in consideration as center of it the center of the first solid.

	G4VSolid *GaContainer = new G4UnionSolid("GaContainer",
											 BackGaContainer,
											 UnionCD,
											 0,
											 G4ThreeVector(0., 0., (H_CylA3 + H_CylC3) * 0.5));

	G4LogicalVolume *logicGaContainer =
		new G4LogicalVolume(GaContainer,	  // its solid
							GaContainer3_mat, // its material
							"GaContainer");	  // its name

	// ###################################################
	//  "CATAFALCO" - Probe Container
	// ##########################

	G4VSolid *CylinderE =
		new G4Tubs("CylinderE", // its name
				   d_CylE3 * 0.5,
				   D_CylEF3 * 0.5,
				   H_CylE3 * 0.5,
				   Ang0,
				   Ang2Pi); // its size

	G4VSolid *CylinderF =
		new G4Tubs("CylinderF", // its name
				   d_CylF3 * 0.5,
				   D_CylEF3 * 0.5,
				   H_CylF3 * 0.5,
				   Ang0,
				   Ang2Pi); // its size

	G4VSolid *CylinderG =
		new G4Tubs("CylinderG", // its name
				   d_CylG3 * 0.5,
				   D_CylG3 * 0.5,
				   H_CylG3 * 0.5,
				   Ang0,
				   Ang2Pi); // its size

	G4VSolid *UnionEF = new G4UnionSolid("UnionEF",
										 CylinderF,
										 CylinderE,
										 0,
										 G4ThreeVector(0., 0., (H_CylE3 + H_CylF3) * 0.5));

	G4RotationMatrix *rm = new G4RotationMatrix();
	rm->rotateY(180. * deg);

	G4VSolid *ProbeContainer = new G4UnionSolid("ProbeContainer",
												CylinderG,
												UnionEF,
												0,
												G4ThreeVector(0., 0., (H_CylG3 + H_CylF3) * 0.5));

	G4LogicalVolume *logicProbeContainer =
		new G4LogicalVolume(ProbeContainer,		// its solid
							ProbeContainer_mat, // its material
							"ProbeContainer");	// its name

	G4ThreeVector ProbeContainerPos = G4ThreeVector(0, 0, 3.31 * cm); // H_C+(H_G/2)+H_F+H_E

	// ###################################################
	//  "CATAFALCO" -  18F Liquid layer around GaSet3 (only placed for source 11)
	// ##########################

	G4Tubs *FExtLiquidSide =
		new G4Tubs("FExtLiquidSide",
				   D_CylG3 * 0.5,
				   D_CylG3 * 0.5 + fSourceDiameter * 0.5,
				   //						 (H_CylA3+H_CylG3+H_CylC3+H_CylD3)*0.5,
				   FExtLiquidLateralHeigth * 0.5,
				   Ang0,
				   Ang2Pi);

	//	G4LogicalVolume* logicFExtLiquidSide =
	//	new G4LogicalVolume(FExtLiquidSide,               //its solid
	//											SourceExtGa_mat,           //its material
	////											world_mat,           //its material
	//											"FExtLiquidSide");            //its name

	G4Tubs *FExtLiquidBack =
		new G4Tubs("FExtLiquidBack",
				   0,
				   D_CylG3 * 0.5 + fSourceDiameter * 0.5,
				   FExtLiquidBottomThickness * 0.5,
				   Ang0,
				   Ang2Pi);

	G4VSolid *FExtLiquidTot = new G4UnionSolid("FExtLiquidTot",
											   FExtLiquidSide,
											   FExtLiquidBack,
											   0,
											   G4ThreeVector(0., 0., -FExtLiquidLateralHeigth * 0.5 - FExtLiquidBottomThickness * 0.5));

	G4LogicalVolume *logicFExtLiquidTot =
		new G4LogicalVolume(FExtLiquidTot,	 // its solid
							SourceExtGa_mat, // its material
							//											world_mat,           //its material
							"FExtLiquidTot"); // its name

	if (fGaSet == 1)
	{
		posAbs = G4ThreeVector(0, 0, fAbsorberThickness / 2.);
	}
	else if (fGaSet == 3)
	{
		posAbs = G4ThreeVector(0, 0, fAbsorberThickness / 2.);
	}

#pragma mark Placement of Volumes

	if (fAbsHoleDiam >= 0)
	{
		G4cout << "GEOMETRY DEBUG - Absorber has been placed!!" << G4endl;

		physAbsorber = new G4PVPlacement(0,				 // no rotation
										 posAbs,		 // at (0,0,0)
										 logicAbsorber,	 // its logical volume
										 "Absorber",	 // its name
										 logicWorld,	 // its mother  volume
										 false,			 // no boolean operation
										 0,				 // copy number
										 checkOverlaps); // overlaps checking

		logicAbsorber->SetRegion(sorgente);
		sorgente->AddRootLogicalVolume(logicAbsorber);
	}
	// ################################################### END OF ABSORBER PLACEMENT

	// ###################################################
	// Front shield
	// ##########################

	Z_FrontShield = Pter_ZScan + FrontShield_sizeZ * 0.5 + FrontShieldBis_sizeZ;

	posFrontShield = G4ThreeVector(fX0Scan, 0, Z_FrontShield);
	G4cout << "GEOMETRY DEBUG - Z thickness of solidFrontShield= " << FrontShield_sizeZ / mm << ", Z pos= " << posFrontShield.z() / mm << G4endl;

	physFrontShield = new G4PVPlacement(0, // no rotation
										posFrontShield,
										logicFrontShield, // its logical volume
										"FrontShield",	  // its name
										logicWorld,		  // its mother  volume
										false,			  // no boolean operation
										0,				  // copy number
										checkOverlaps);	  // overlaps checking
	logicFrontShield->SetRegion(frontshieldreg);
	frontshieldreg->AddRootLogicalVolume(logicFrontShield);

	// ###################################################
	// Front shield BIS
	// ##########################

	Z_FrontShieldBis = Pter_ZScan + FrontShieldBis_sizeZ * 0.5;

	posFrontShieldBis = G4ThreeVector(fX0Scan, 0, Z_FrontShieldBis);
	G4cout << "GEOMETRY DEBUG - Z thickness of solidFrontShieldBis= " << FrontShieldBis_sizeZ / mm << ", Z pos= " << posFrontShieldBis.z() / mm << G4endl;

	if (secondShieldFlag)
		physFrontShieldBis = new G4PVPlacement(0, // no rotation
											   posFrontShieldBis,
											   logicFrontShieldBis, // its logical volume
											   "FrontShieldBis",	// its name
											   logicWorld,			// its mother  volume
											   false,				// no boolean operation
											   0,					// copy number
											   checkOverlaps);		// overlaps checking
	logicFrontShieldBis->SetRegion(frontshieldreg);
	frontshieldreg->AddRootLogicalVolume(logicFrontShieldBis);

	// ###################################################
	//  	P-Terphenyl
	// ##########################

	Pter_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ * 0.5 + FrontShieldBis_sizeZ;

	G4ThreeVector posPter = G4ThreeVector(fX0Scan, 0, Pter_Posz);
	G4cout << "GEOMETRY DEBUG - Z thickness of Pterp= " << Pter_sizeZ / mm << ", Z pos= " << posPter.z() / mm << G4endl;

	physPter = new G4PVPlacement(0, // no rotation
								 posPter,
								 logicPter,		 // its logical volume
								 "Pter",		 // its name
								 logicWorld,	 // its mother  volume
								 false,			 // no boolean operation
								 0,				 // copy number
								 checkOverlaps); // overlaps checking

	logicPter->SetRegion(pterreg);
	pterreg->AddRootLogicalVolume(logicPter);

	fScoringVolume = logicPter;

	// ###################################################
	//  SiPm volume behind PTER
	// ##########################

	G4ThreeVector posSiPm = G4ThreeVector(fX0Scan, 0, Pter_Posz + Pter_sizeZ / 2. + DzSiPm / 2. + FrontShieldBis_sizeZ);

	physSiPm = new G4PVPlacement(0,				 // no rotation
								 posSiPm,		 // at (0,0,0)
								 logicSiPm,		 // its logical volume
								 "SiPm",		 // its name
								 logicWorld,	 // its mother  volume
								 false,			 // no boolean operation
								 0,				 // copy number
								 checkOverlaps); // overlaps checking

	// ###################################################
	//  PVC around P-Terphenyl
	// ##########################
	physPVC = new G4PVPlacement(0, // no rotation
								posPter,
								logicPVC,		// its logical volume
								"PVC",			// its name
								logicWorld,		// its mother  volume
								false,			// no boolean operation
								0,				// copy number
								checkOverlaps); // overlaps checking

	// ###################################################
	//  Delrin around P-Terphenyl
	// ##########################

	physDelrin = new G4PVPlacement(0, // no rotation
								   posPter,
								   logicDelrin, // its logical volume
								   "Delrin",	// its name
								   logicWorld,	// its mother  volume
								   false,		// no boolean operation
								   0,			// copy number
								   checkOverlaps);

	if (CaseDepth > 0)
	{
#pragma mark Robotic Probe
		// ################ NL probe case for laparoscopic surgery
		// ################ Top Probe Case
		G4Tubs *solidTopCase =
			new G4Tubs("solidTopCase",
					   0.,
					   PVC_outer_r - LateralCaseThickness,
					   TopCaseDepth * 0.5,
					   Ang0,
					   Ang2Pi);

		G4Tubs *solidAroundTopCase =
			new G4Tubs("AroundTopCase",
					   PVC_outer_r - LateralCaseThickness,
					   PVC_outer_r,
					   TopCaseDepth * 0.5,
					   Ang0,
					   Ang2Pi);

		// ################ Middle Case
		G4double MiddleCaseDepth = (CaseDepth - HorsesShoeBackThickness - BackCaseThickness);

		G4Tubs *solidMiddleCase =
			new G4Tubs("MiddleCase",
					   0.,
					   PVC_outer_r - LateralCaseThickness - HorsesShoeLateralThickness,
					   MiddleCaseDepth * 0.5,
					   Ang0,
					   Ang2Pi);

		G4Tubs *solidAroundMiddleCase =
			new G4Tubs("AroundMiddleCase",
					   PVC_outer_r - LateralCaseThickness - HorsesShoeLateralThickness,
					   PVC_outer_r - LateralCaseThickness,
					   MiddleCaseDepth * 0.5,
					   Ang0,
					   Ang2Pi);

		G4Tubs *solidExtMiddleCase =
			//		new G4Tubs("AroundMiddleCase", //corrected on 2019.01.10, probably was a cut and paste error
			new G4Tubs("ExtMiddleCase",
					   PVC_outer_r - LateralCaseThickness,
					   PVC_outer_r,
					   MiddleCaseDepth * 0.5,
					   Ang0,
					   Ang2Pi);

		G4LogicalVolume *logicMiddleCase =
			new G4LogicalVolume(solidMiddleCase, // its solid
								MiddleCase_mat,	 // its material
								"MiddleCase");	 // its name

		// ################ Back Middle Case
		G4Tubs *solidBackMiddleCase =
			new G4Tubs("BackMiddleCase",
					   0.,
					   PVC_outer_r - LateralCaseThickness,
					   HorsesShoeBackThickness * 0.5,
					   Ang0,
					   Ang2Pi);

		G4Tubs *solidAroundBackMiddleCase =
			new G4Tubs("AroundBackMiddleCase",
					   PVC_outer_r - LateralCaseThickness,
					   PVC_outer_r,
					   HorsesShoeBackThickness * 0.5,
					   Ang0,
					   Ang2Pi);

		// ################ End Case
		G4Tubs *solidEndCase =
			new G4Tubs("EndCase",
					   0.,
					   PVC_outer_r,
					   BackCaseThickness * 0.5,
					   Ang0,
					   Ang2Pi);

		// ###################################################
		//  Middle Probe case
		// ##########################
		G4double ProbeMiddleCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth * 0.5;
		G4ThreeVector posMiddleCase = G4ThreeVector(fX0Scan, 0, ProbeMiddleCase_Posz);

		new G4PVPlacement(0,			   // no rotation
						  posMiddleCase,   // at (0,0,0)
						  logicMiddleCase, // its logical volume
						  "MiddleCase",	   // its name
						  logicWorld,	   // its mother  volume
						  false,		   // no boolean operation
						  0,			   // copy number
						  checkOverlaps);  // overlaps checking

		// ###################################################
		//  External Plastic Case
		// ##########################

		G4VSolid *PlasticCasePartialUnion1 = new G4UnionSolid("PlasticCasePartialUnion1", solidAroundTopCase, solidExtMiddleCase, 0, G4ThreeVector(0., 0., (MiddleCaseDepth + TopCaseDepth) * 0.5));
		// The G4ThreeVector(0.,0.,(MiddleCaseDepth+TopCaseDepth)*0.5) generates a translation of the second solid respect to the center of the first one. When I place the resulting solid I've to take in consideration as center of it the center of the first solid.

		G4VSolid *PlasticCasePartialUnion2 = new G4UnionSolid("PlasticCasePartialUnion2", solidAroundBackMiddleCase, solidEndCase, 0, G4ThreeVector(0., 0., (BackCaseThickness + HorsesShoeBackThickness) * 0.5));

		G4VSolid *PlasticCase = new G4UnionSolid("PlasticCase", PlasticCasePartialUnion1, PlasticCasePartialUnion2, 0, G4ThreeVector(0., 0., MiddleCaseDepth + TopCaseDepth * 0.5 + HorsesShoeBackThickness * 0.5));

		G4LogicalVolume *logicPlasticCase =
			new G4LogicalVolume(PlasticCase,	 // its solid
								PlasticCase_mat, // its material
								"PlasticCase");	 // its name

		G4double PlasticCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth * 0.5;
		G4ThreeVector posPlasticCase = G4ThreeVector(fX0Scan, 0, PlasticCase_Posz);

		new G4PVPlacement(0,				// no rotation
						  posPlasticCase,	// at (0,0,0)
						  logicPlasticCase, // its logical volume
						  "PlasticCase",	// its name
						  logicWorld,		// its mother  volume
						  false,			// no boolean operation
						  0,				// copy number
						  checkOverlaps);	// overlaps checking

		// ###################################################
		//  G4Union HorsesShoe Probe
		// ##########################

		G4VSolid *HorsesShoe = new G4UnionSolid("HorsesShoe", solidAroundMiddleCase, solidBackMiddleCase, 0, G4ThreeVector(0., 0., (MiddleCaseDepth + HorsesShoeBackThickness) * 0.5));

		G4LogicalVolume *logicHorsesShoe =
			new G4LogicalVolume(HorsesShoe,		// its solid
								HorsesShoe_mat, // its material
								"HorsesShoe");	// its name

		// G4double HorsesShoe_Posz = DzAbs*0.5+ Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth*0.5;
		G4double HorsesShoe_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth * 0.5;

		G4ThreeVector posHorsesShoe = G4ThreeVector(fX0Scan, 0, HorsesShoe_Posz);

		new G4PVPlacement(0,			   // no rotation
						  posHorsesShoe,   // at (0,0,0)
						  logicHorsesShoe, // its logical volume
						  "HorsesShoe",	   // its name
						  logicWorld,	   // its mother  volume
						  false,		   // no boolean operation
						  0,			   // copy number
						  checkOverlaps);  // overlaps checking

		// ###################################################
		//  G4Union TopCase Probe
		// ###################################################

		G4VSolid *TopCase =
			new G4SubtractionSolid("TopCase",
								   solidTopCase,
								   solidSiPm,
								   0,
								   G4ThreeVector(0., 0., TopCaseDepth * 0.5 - DzSiPm * 0.5));

		G4LogicalVolume *logicTopCase =
			new G4LogicalVolume(TopCase,	 // its solid
								TopCase_mat, // its material
								"TopCase");	 // its name

		// G4double ProbeTopCase_Posz = DzAbs*0.5+ Pter_ZScan + DzDummy2 + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
		G4double ProbeTopCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth * 0.5;

		G4ThreeVector posTopCase = G4ThreeVector(fX0Scan, 0, ProbeTopCase_Posz);

		G4RotationMatrix *RotMatrix = new G4RotationMatrix();
		RotMatrix->rotateY(180 * deg);

		new G4PVPlacement(RotMatrix,	  // rotation
						  posTopCase,	  // at (0,0,0)
						  logicTopCase,	  // its logical volume
						  "TopCase",	  // its name
						  logicWorld,	  // its mother  volume
						  false,		  // no boolean operation
						  0,			  // copy number
						  checkOverlaps); // overlaps checking

	} // finish placing robotic probe
	else if (CaseDepth < 0)
	{
#pragma mark Open Surgery Probe
		G4double AluCaseDepth = fabs(CaseDepth);

		G4double Alu_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth - 1.15) * 0.5 * mm + FrontShieldBis_sizeZ;
		//		G4double Alu_Posz = DzDummy2 + Pter_ZScan  + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth-1.15)*0.5*mm;
		G4ThreeVector posAlu = G4ThreeVector(0, 0, Alu_Posz);

		G4double BackAlu_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth - 1.15) * mm + 1.15 * 0.5 * mm + FrontShieldBis_sizeZ;
		//		G4double BackAlu_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth-1.15)*mm + 1.15*0.5*mm;
		G4ThreeVector posBackAlu = G4ThreeVector(0, 0, BackAlu_Posz);

		G4Tubs *solidAlu =
			new G4Tubs("AluCase",
					   9.7 * 0.5 * mm,
					   12. * 0.5 * mm,
					   (AluCaseDepth - 1.15) * 0.5 * mm,
					   Ang0,
					   Ang2Pi);

		G4Tubs *solidBackAlu =
			new G4Tubs("BackAluCase",
					   0.,
					   12. * 0.5 * mm,
					   1.15 * 0.5 * mm,
					   Ang0,
					   Ang2Pi);

		G4LogicalVolume *logicAlu =
			new G4LogicalVolume(solidAlu,
								FrontShield_mat,
								"AluCase");

		G4LogicalVolume *logicBackAlu =
			new G4LogicalVolume(solidBackAlu,
								FrontShield_mat,
								"BackAluCase");

		new G4PVPlacement(0,
						  posAlu,
						  logicAlu,
						  "AluCase",
						  logicWorld,
						  false,
						  0,
						  checkOverlaps);

		new G4PVPlacement(0,
						  posBackAlu,
						  logicBackAlu,
						  "BackAluCase",
						  logicWorld,
						  false,
						  0,
						  checkOverlaps);
	} // finish placing open surgery probe

	// ##################################################################
	// ################################################
	// ######################## PLACEMENTS
#pragma mark Placement of Sources etc
	//	G4ThreeVector posPter = G4ThreeVector(fX0Scan, 0, Pter_Posz);

	if (fGaSet == 1)
	{ // no "catafalco"

		if (fSourceSelect == 8 || fSourceSelect == 9 || fSourceSelect < 0)
		{ // If I requested the Cu67/F18/GenericIon source (or the flat electron one for efficiencies)
			G4cout << "GEOMETRY DEBUG - Cu/F Source has been placed!!" << G4endl;

			new G4PVPlacement(0,					 // no rotation
							  posSourceGenericIon,	 // at (0,0,0)
							  logicSourceGenericIon, // its logical volume
							  "Source",				 // its name
							  logicWorld,			 // its mother  volume
							  false,				 // no boolean operation
							  0,					 // copy number
							  checkOverlaps);		 // overlaps checking

			logicSourceGenericIon->SetRegion(sorgente);
			sorgente->AddRootLogicalVolume(logicSourceGenericIon);
		}

		if (fSourceSelect == 3)
		{ // Y extended source (AGAR AGAR)

			G4cout << "GEOMETRY DEBUG - Z thickness of solidSourceExtY= " << DzSourceExtY / mm << ", Z pos= " << posSourceExtY.z() / mm << G4endl;
			G4cout << "GEOMETRY DEBUG - ExtY Source has been placed!!" << G4endl;

			new G4PVPlacement(0, // no rotation
							  posSourceExtY,
							  logicSourceExtY, // its logical volume
							  "Source",		   // its name
							  logicWorld,	   // its mother  volume
							  false,		   // no boolean operation
							  0,			   // copy number
							  checkOverlaps);  // overlaps checking

			logicSourceExtY->SetRegion(sorgente);
			sorgente->AddRootLogicalVolume(logicSourceExtY);

			G4cout << "GEOMETRY DEBUG - Z thickness of solidABSaround= " << DzABSaround / mm << ", Z pos= " << posABSaround.z() / mm << G4endl;
			G4cout << "GEOMETRY DEBUG - ExtYTOC Source has been placed!!" << G4endl;

			new G4PVPlacement(0,			  // no rotation
							  posABSaround,	  // at (0,0,0)
							  logicABSaround, // its logical volume
							  "ABSaround",	  // its name
							  logicWorld,	  // its mother  volume
							  false,		  // no boolean operation
							  0,			  // copy number
							  checkOverlaps); // overlaps checking

			logicABSaround->SetRegion(ABSRegion);
			ABSRegion->AddRootLogicalVolume(logicABSaround);

			G4cout << "GEOMETRY DEBUG - Z thickness of solidABSbehind= " << DzABSbehind / mm << ", Z pos= " << posABSbehind.z() / mm << G4endl;
			G4cout << "GEOMETRY DEBUG - ExtYTOC Source has been placed!!" << G4endl;

			new G4PVPlacement(0,			  // no rotation
							  posABSbehind,	  // at (0,0,0)
							  logicABSbehind, // its logical volume
							  "ABSbehind",	  // its name
							  logicWorld,	  // its mother  volume
							  false,		  // no boolean operation
							  0,			  // copy number
							  checkOverlaps); // overlaps checking

			logicABSbehind->SetRegion(ABSRegion);
			ABSRegion->AddRootLogicalVolume(logicABSbehind);

			G4cout << "GEOMETRY DEBUG - Z thickness of solidABSbehind= " << DzABSbehind / mm << ", Z pos= " << posABSbehind.z() / mm << G4endl;
			G4cout << "GEOMETRY DEBUG - ExtYTOC Source has been placed!!" << G4endl;

			new G4PVPlacement(0,			  // no rotation
							  posABSbehind,	  // at (0,0,0)
							  logicABSbehind, // its logical volume
							  "ABSbehind",	  // its name
							  logicWorld,	  // its mother  volume
							  false,		  // no boolean operation
							  0,			  // copy number
							  checkOverlaps); // overlaps checking

			logicABSbehind->SetRegion(ABSRegion);
			ABSRegion->AddRootLogicalVolume(logicABSbehind);
		}
		if (fSourceSelect == 1 || fSourceSelect == 2 || fSourceSelect == 6 || fSourceSelect == 13)
		{ // If I requested the Sr source (or the flat electron one for efficiencies)
			G4cout << "GEOMETRY DEBUG - Sr(-like) Source has been placed!!" << G4endl;

			new G4PVPlacement(0,			  // no rotation
							  posSourceSR,	  // at (0,0,0)
							  logicSourceSR,  // its logical volume
							  "Source",		  // its name
							  logicWorld,	  // its mother  volume
							  false,		  // no boolean operation
							  0,			  // copy number
							  checkOverlaps); // overlaps checking

			logicSourceSR->SetRegion(sorgente);
			sorgente->AddRootLogicalVolume(logicSourceSR);
		}

		if (fSourceSelect == 12)
		{ // If I requested the GA point source
			G4cout << "GEOMETRY DEBUG - GA PointLike Source has been placed!!" << G4endl;

			new G4PVPlacement(0,				  // no rotation
							  posSourceGaPoint,	  // at (0,0,0)
							  logicSourceGaPoint, // its logical volume
							  "Source",			  // its name
							  logicWorld,		  // its mother  volume
							  false,			  // no boolean operation
							  0,				  // copy number
							  checkOverlaps);	  // overlaps checking

			logicSourceGaPoint->SetRegion(sorgente);
			sorgente->AddRootLogicalVolume(logicSourceGaPoint);
		}
	}
	else if (fGaSet == 3)
	{ // with "catafalco"

		new G4PVPlacement(rm,				   // no rotation
						  ProbeContainerPos,   // at (0,0,0)
						  logicProbeContainer, // its logical volume
						  "ProbeContainer",	   // its name
						  logicWorld,		   // its mother  volume
						  false,			   // no boolean operation
						  0,				   // copy number
						  checkOverlaps);

		new G4PVPlacement(0,				 // no rotation
						  posContainerExtGa, // at (0,0,0)
						  logicGaContainer,	 // its logical volume
						  "GaContainer",	 // its name
						  logicWorld,		 // its mother  volume
						  false,			 // no boolean operation
						  0,				 // copy number
						  checkOverlaps);

		if (fSourceSelect == 11)
		{ // liquid around "catafalco"
			G4ThreeVector posFLiquidBack = posPter + G4ThreeVector(0, 0, -(H_CylA3 + H_CylG3 + H_CylC3 + H_CylD3) * 0.5 - 5 * mm);

			new G4PVPlacement(0, // no rotation
							  G4ThreeVector(0, 0,
											FExtLiquidLateralHeigth * 0.5 - H_CylA3 - DzDummyExitSorg),
							  logicFExtLiquidTot, // its logical volume
							  "FExtLiquidTot",	  // its name
							  logicWorld,		  // its mother  volume
							  false,			  // no boolean operation
							  0,				  // copy number
							  checkOverlaps);	  // overlaps checking
		}
		else
		{
			new G4PVPlacement(0,				 // no rotation
							  posExtGa3,		 // at (0,0,0)
							  logicSourceExtGa2, // its logical volume
							  "Source",			 // its name
							  logicWorld,		 // its mother  volume
							  false,			 // no boolean operation
							  0,				 // copy number
							  checkOverlaps);	 // overlaps checking
		}
		// ###################################################
		//  Table
		// ##########################

		G4ThreeVector posTable = G4ThreeVector(0, 0, -H_CylA3 - DzTable / 2. - DzDummyExitSorg);
		//		if (fSourceSelect==11) posTable = G4ThreeVector(0, 0, -H_CylA3 - DzTable/2.-DzDummyExitSorg-FExtLiquidBottomThickness-1*um); //per fare spazio allo strato di F18 sotto il contenitore
		if (fSourceSelect == 11)
			posTable = G4ThreeVector(0, 0, -H_CylA3 - DzTable / 2. - DzDummyExitSorg - FExtLiquidBottomThickness); // per fare spazio allo strato di F18 sotto il contenitore
		new G4PVPlacement(0,																					   // no rotation
						  posTable,																				   // at (0,0,0)
						  logicTable,																			   // its logical volume
						  "Table",																				   // its name
						  logicWorld,																			   // its mother  volume
						  false,																				   // no boolean operation
						  0,																					   // copy number
						  checkOverlaps);																		   // overlaps checking
		//
		//		//###################################################
		//		// PVC around P-Terphenyl
		//		//##########################
		//
		//		physPVC = new G4PVPlacement(0,                     //no rotation
		//																posPter,
		//																logicPVC,            //its logical volume
		//																"PVC",               //its name
		//																logicWorld,            //its mother  volume
		//																false,                 //no boolean operation
		//																0,                     //copy number
		//																checkOverlaps);        //overlaps checking
		//
		//
		//		//###################################################
		//		// Delrin around P-Terphenyl
		//		//##########################
		//
		//		physDelrin = new G4PVPlacement(0,                     //no rotation
		//																	 posPter,
		//																	 logicDelrin,            //its logical volume
		//																	 "Delrin",               //its name
		//																	 logicWorld,            //its mother  volume
		//																	 false,                 //no boolean operation
		//																	 0,                     //copy number
		//																	 checkOverlaps);        //overlaps checking
		//
	} // end of GaSet3

	if (fGaSet == 1)
	{
		posDummyExitSorg = G4ThreeVector(0, 0, -DzDummyExitSorg / 2.);
		posDummyExitAbs = G4ThreeVector(0, 0, fAbsorberThickness + DzDummyExitAbs / 2.);
		G4cout << "CIAONE debug: prima: " << fAbsorberThickness + DzDummyExitAbs / 2. << " dopo " << posAbs.z() + fAbsorberThickness / 2. + DzDummyExitAbs / 2. << G4endl;
		posDummyEnterProbe = G4ThreeVector(0, 0, posFrontShield.z() - FrontShield_sizeZ / 2. - DzDummyEnterProbe / 2. - FrontShieldBis_sizeZ);
	}
	if (fGaSet == 3)
	{
		posDummyExitSorg = G4ThreeVector(0, 0, -DzDummyExitSorg / 2.);
		//		posDummyExitAbs = G4ThreeVector(0, 0, posAbs.z() +fAbsorberThickness/2.+ DzDummyExitAbs/2.);
		posDummyExitAbs = G4ThreeVector(0, 0, fAbsorberThickness + DzDummyExitAbs / 2.);
		G4cout << "CIAONE geodebug: PosAbsDummy prima: " << posDummyExitAbs.z() << G4endl;
		posDummyEnterProbe = G4ThreeVector(0, 0, posFrontShield.z() - FrontShield_sizeZ / 2. - DzDummyEnterProbe / 2.);
	}

	if (fSourceSelect != 5 && fSourceSelect != 7 && fSourceSelect != 10)
		new G4PVPlacement(0, posDummyExitSorg, logicShapeDummyExitSorg, "DummyExitSorg", logicWorld, false, 0, checkOverlaps); // overlaps checking
	if (fAbsHoleDiam >= 0)
		new G4PVPlacement(0, posDummyExitAbs, logicShapeDummyExitAbs, "DummyExitAbs", logicWorld, false, 0, checkOverlaps);		 // overlaps checking
	new G4PVPlacement(0, posDummyEnterProbe, logicShapeDummyEnterProbe, "DummyEnterProbe", logicWorld, false, 0, checkOverlaps); // overlaps checking

	logicShapeDummyExitSorg->SetRegion(sorgente);
	sorgente->AddRootLogicalVolume(logicShapeDummyExitSorg);

	logicShapeDummyExitAbs->SetRegion(sorgente);
	sorgente->AddRootLogicalVolume(logicShapeDummyExitAbs);

	logicShapeDummyEnterProbe->SetRegion(sorgente);
	sorgente->AddRootLogicalVolume(logicShapeDummyEnterProbe);

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
