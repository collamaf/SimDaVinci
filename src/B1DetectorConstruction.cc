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

B1DetectorConstruction::B1DetectorConstruction(G4double x0, G4double ZValue, G4double CuDiam, G4int SourceSelect, G4int AbsorberMaterial,G4double PterDiameter, G4double PterThickness,G4double SourceDiameter,G4double SourceThickness, G4double AbsorberThickness, G4double ProbeCaseDepth, G4double ProbeCaseLateralThickness, G4double ProbeCaseBackThickness, G4double HSLateralThickness, G4double HSBackThickness, G4int HousingCase, G4bool ScintFlag, G4int GaSet, G4int ApparatusMat,G4int PosAbsorber,G4double AbsCenter)
: G4VUserDetectorConstruction(),
fScoringVolume(0), fX0Scan(x0), fZValue(ZValue), fCuDiam(CuDiam), fSourceSelect(SourceSelect), fAbsorberMaterial(AbsorberMaterial), fPterDiameter(PterDiameter), fPterThickness(PterThickness), fSourceDiameter(SourceDiameter), fSourceThickness(SourceThickness), fAbsorberThickness(AbsorberThickness),fCaseDepth(ProbeCaseDepth),fLateralCaseThickness(ProbeCaseLateralThickness), fBackCaseThickness(ProbeCaseBackThickness), fHorsesShoeLateralThickness(HSLateralThickness),fHorsesShoeBackThickness(HSBackThickness), fHousingCase(HousingCase), fScintFlag(ScintFlag), fGaSet(GaSet), fApparatusMat (ApparatusMat), fPosAbsorber (PosAbsorber), fAbsCenter (AbsCenter)
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
	G4bool checkOverlaps = true;
	
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
	
	
	//###################################################
	// GaContainer2 Material
	//##########################
	
	G4double GaContainer2Matdensity = 0.43*g/cm3;
	G4Material* GaContainer2Mat = new G4Material (name="GaContainer2Mat", GaContainer2Matdensity, ncomponents=3);
	GaContainer2Mat->AddElement (elH, natoms=17);
	GaContainer2Mat->AddElement (elC, natoms=15);
	GaContainer2Mat->AddElement (elN, natoms=1);
	
	//############ MATERIAL ASSIGNMENT
	G4Material* SourceExtY_mat = AgarAgar;
	G4Material* ABSaround_mat = ABS;
	G4Material* ABSbehind_mat = ABS;
	G4Material* GaContainer_mat=nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material* GaContainer2_mat = GaContainer2Mat;
	G4Material* GaContainer3_mat =nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material* ProbeContainer_mat=nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	//G4Material* GaWater_mat=nist->FindOrBuildMaterial("G4_WATER");
	G4Material* SourceExtGa_mat=nist->FindOrBuildMaterial("G4_WATER");
	G4Material* SourceSR_mat = nist->FindOrBuildMaterial("MyAlu"); //G4_Al
	G4Material* FrontShield_mat = nist->FindOrBuildMaterial("MyAlu");
	G4Material* shapeDummy_mat = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* Pter_mat = PTerphenyl;
	G4Material* PVC_mat= nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material* Delrin_mat=Delrin;
	G4Material* shapeCo_mat = nist->FindOrBuildMaterial("G4_Cu");
	G4Material* Absorber_mat = nist->FindOrBuildMaterial("G4_Cu");
	G4Material* SiPM_mat = nist->FindOrBuildMaterial("G4_Si");
	G4Material* Table_mat = nist->FindOrBuildMaterial("MyAlu");


	
	G4Material* PlasticCase_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Material* HorsesShoe_mat= nist->FindOrBuildMaterial("G4_Pb");
	G4Material* CaseInner_mat=nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	// To make them of air
	if (fHousingCase==2) //1- InnerPlastic (norm): Pb+Pl || 2- InnerAir: Pb+Air || 3- TotalAir: Air+Air
		CaseInner_mat=world_mat;
	else if (fHousingCase==3) {
		HorsesShoe_mat=world_mat;
		CaseInner_mat=world_mat;
	}
	//
	if (fSourceSelect==6 || fSourceSelect==7) SourceSR_mat=world_mat;

	
	// To make -GaSet 2  expreimental environment of air or of Pb
	if (fApparatusMat==2 && fGaSet==2){
		ProbeContainer_mat = world_mat;
		GaContainer2_mat = world_mat;
	}else if (fApparatusMat==3 && fGaSet==2){
		ProbeContainer_mat = HorsesShoe_mat;
		GaContainer2_mat = HorsesShoe_mat;
	}
	
	
	// To make -GaSet 3  expreimental environment of air or of Pb
	if (fApparatusMat==2 && fGaSet==3){
		ProbeContainer_mat = world_mat;
		GaContainer3_mat = world_mat;
	}else if (fApparatusMat==3 && fGaSet==3){
		ProbeContainer_mat = HorsesShoe_mat;
		GaContainer3_mat = HorsesShoe_mat;
	}
	
	
	
	G4Material* MiddleCase_mat=CaseInner_mat;
	G4Material* TopCase_mat=CaseInner_mat;
	


	
	
	
	//###################################################
	//###################################################
	// Optics characteristics part
	//##########################
	//
	// ------------ Generate & Add Material Properties Table ------------
	//
	G4double photonEnergy[] =
	{ 2.96*eV};
	
	const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);
	
	//
	// Water
	//
	G4double refractiveIndex1[] =
	{ 1.65}; //da misteriosa mail del 30.9.2013
	
	assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));
	
	G4double absorption[] =
	{20*mm }; //da elsarticle CMT
	
	assert(sizeof(absorption) == sizeof(photonEnergy));
	
	G4double scintilFast[] =
	{ 1.00};
	
	assert(sizeof(scintilFast) == sizeof(photonEnergy));
	
	G4double scintilSlow[] =
	{ 0.0 };
	
	assert(sizeof(scintilSlow) == sizeof(photonEnergy));
	
	
	
	G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
	myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
	->SetSpline(true);
	myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
	->SetSpline(true);
	myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
	->SetSpline(true);
	myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
	->SetSpline(true);
	
	myMPT1->AddConstProperty("SCINTILLATIONYIELD",28000./MeV); //33k da nostro papero, 28k da papero recente elsa CMT
	myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
	myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
	myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
	myMPT1->AddConstProperty("YIELDRATIO",1);

	
	
	G4cout << "PTERP G4MaterialPropertiesTable" << G4endl;
	myMPT1->DumpTable();
	
	if (fScintFlag)
	PTerphenyl->SetMaterialPropertiesTable(myMPT1); //to toggle scintillation
	
	//##########################
	//###################################################

	
	
	
	
	
	
	
	
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
	G4double FrontShield_sizeZ=15*um;
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
	
	//### Absorber
	G4double RminAbs = fabs(fCuDiam)/2.*mm;
	G4double RmaxAbs = 0*mm;
	G4double DzAbs= fAbsorberThickness*mm;
	if(fPosAbsorber==1 && fGaSet==2){
		RmaxAbs = 21/2.*mm;
	}else if (fPosAbsorber==1 && fGaSet==3){
		RmaxAbs = 22/2.*mm;
	}else if (fPosAbsorber==2 && fGaSet==2){
		RmaxAbs = 26/2.*mm;
	}else if (fPosAbsorber==2 && fGaSet==3){
		RmaxAbs = 26/2.*mm;
	};
	G4double SPhiAbs = 0.*deg;
	G4double DPhiAbs = 360.*deg;
	G4double ZCenterAbs=fAbsCenter*mm;
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
	
	
	//### Dummy2
	G4double RminDummy2 = 0.*mm;
	G4double RmaxDummy2 = 10.5*mm;
	 if (fCuDiam>=0 && fPosAbsorber==2) {
		RmaxDummy2= 13. *mm;
	}else if(fCuDiam>=0 && fPosAbsorber==1 && fGaSet==2) {
		RmaxDummy2 = 10.5*mm;
	}else if(fCuDiam>=0 && fPosAbsorber==1 && fGaSet==3) {
		RmaxDummy2 = 11.*mm;
	}else if(fCuDiam<0 && fPosAbsorber==1 && fGaSet==3) {
		RmaxDummy2 = 11.*mm;
	}else RmaxDummy2=10.5*mm;
	/*
	if (fGaSet==3) {
		RmaxDummy2= 11.*mm;
	}else RmaxDummy2 = 10.5*mm;
	 */
	G4double DzDummy2= 1.e-5*mm;
	G4double SPhiDummy2 = 0.*deg;
	G4double DPhiDummy2 = 360.*deg;
	G4double zDummy2;
	//###
	
	//### Dummy3
	G4double RminDummy3 = 0.*mm;
	G4double RmaxDummy3 = 5.*mm;
	G4double DzDummy3= 7.5*mm-fSourceThickness*mm;
	G4double SPhiDummy3 = 0.*deg;
	G4double DPhiDummy3 = 360.*deg;
	G4double zDummy3;
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

	
	//### GaSet 2

	G4double D_CylABCD = 70*mm;
	G4double H_CylA = 19*mm;
	G4double d_CylB = 10*mm;
	G4double H_CylB = 7*mm;
	G4double d_CylC = 21*mm;
	G4double H_CylC = 2*mm;
	G4double d_CylD = 48*mm;
	G4double H_CylD = 14*mm;
	
	G4double D_CylEF = 48*mm;
	G4double d_CylE = 26*mm;
	G4double H_CylE = 6*mm;
	G4double d_CylF = 13*mm;
	G4double H_CylF = 8*mm;
	G4double D_CylG = 51*mm;
	G4double d_CylG = 13*mm;
	G4double H_CylG = 32*mm;

	G4double SPhiCyl = 0.*deg;
	G4double DPhiCyl = 360.*deg;
	
	//### GaSet 3
	
	G4double D_CylABCD3 = 50.8*mm;
	G4double H_CylA3 = 34.7*mm;
	G4double d_CylB3 = 10*mm;
	G4double H_CylB3 = 7.5*mm;
	G4double d_CylC3 = 22*mm;
	G4double H_CylC3 = 3.3*mm;
	G4double d_CylD3 = 41.3*mm;
	G4double H_CylD3 = 13.8*mm;
	
	G4double D_CylEF3 = 41.3*mm;
	G4double d_CylE3 = 26*mm;
	G4double H_CylE3 = 5.45*mm;
	G4double d_CylF3 = 13*mm;
	G4double H_CylF3 = 8*mm;
	G4double D_CylG3 = 51*mm;
	G4double d_CylG3 = 13*mm;
	G4double H_CylG3 = 32*mm;
	
	G4double SPhiCyl3 = 0.*deg;
	G4double DPhiCyl3 = 360.*deg;
	
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
	
	if(fSourceSelect==1 || fSourceSelect==2 || fSourceSelect==5 || fSourceSelect==6) { //If i requested the Sr source
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
	// Probe's Solids
	//##########################
	
	
	//################ Dummy

	
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
	

	
	//###########################
	//Dummy2
	//###########################
	
	
	G4Tubs* solidShapeDummy2 =
	new G4Tubs("Dummy2",                       //its name
						 RminDummy2,
						 RmaxDummy2,
						 0.5*DzDummy2,
						 SPhiDummy2,
						 DPhiDummy2);     //its size
	
	G4LogicalVolume* logicShapeDummy2 =
	new G4LogicalVolume(solidShapeDummy2,          //its solid
											shapeDummy_mat,           //its material
											"Dummy2");            //its name
	
	
	//###########################
	//Dummy3
	//###########################
	
	
	G4Tubs* solidShapeDummy3 =
	new G4Tubs("Dummy3",                       //its name
						 RminDummy3,
						 RmaxDummy3,
						 0.5*DzDummy3,
						 SPhiDummy3,
						 DPhiDummy3);     //its size
	
	G4LogicalVolume* logicShapeDummy3 =
	new G4LogicalVolume(solidShapeDummy3,          //its solid
											shapeDummy_mat,           //its material
											"Dummy3");            //its name
	
	
	//################ Front-Shield

	
	
	G4Tubs* solidFrontShield =
	new G4Tubs("FrontShield",                       //its name
						 0.,FrontShield_outer_r,FrontShield_sizeZ*0.5,FrontShield_start_angle,FrontShield_spanning_angle);     //its size
	
	
	
	G4LogicalVolume* logicFrontShield =
	new G4LogicalVolume(solidFrontShield,          //its solid
											FrontShield_mat,           //its material
											"FrontShield");            //its name
	
	
	//################ Pter

	
	G4Tubs* solidPter =
	new G4Tubs("Pter",                                                                         //its name
						 0.,Pter_Diam*0.5,Pter_sizeZ*0.5,Pter_start_angle,Pter_spanning_angle);                //its size
	
	G4LogicalVolume* logicPter =
	new G4LogicalVolume(solidPter,          //its solid
											Pter_mat,           //its material
											"Pter");            //its name
	
	
	//################ SiPM

	
	G4double DxSiPm = 3.*mm;
	G4double DySiPm= 3.*mm;
	G4double DzSiPm = 300.e-6*mm;
	
	G4Box* solidSiPm =
	new G4Box("SiPm",                       //its name
						DxSiPm/2.,
						DySiPm/2.,
						DzSiPm/2.);     //its size
	
	G4LogicalVolume* logicSiPm =
	new G4LogicalVolume(solidSiPm,          //its solid
											SiPM_mat,           //its material
											"SiPm");            //its name
	
	
	//################ Table
	
	
	G4double DxTable = 500*mm;
	G4double DyTable= 500*mm;
	G4double DzTable = 50*mm;
	
	G4Box* solidTable =
	new G4Box("Table",                       //its name
						DxTable/2.,
						DyTable/2.,
						DzTable/2.);     //its size
	
	G4LogicalVolume* logicTable =
	new G4LogicalVolume(solidTable,          //its solid
											Table_mat,           //its material
											"Table");            //its name
	
	
	
	//################ PVC Around Pter

	
	G4Tubs* solidPVC =
	new G4Tubs("PVC",                                                                         //its name
						 PVC_inner_r,PVC_outer_r,Pter_sizeZ*0.5,Pter_start_angle,Pter_spanning_angle);                //its size
	
	G4LogicalVolume* logicPVC =
	new G4LogicalVolume(solidPVC,          //its solid
											PVC_mat,              //its material
											"PVC");            //its name
	
	
	
	
	//################ Delrin Around Pter

	
	G4Tubs* solidDelrin =
	new G4Tubs("Delrin",                                                                         //its name
						 Pter_Diam*0.5,PVC_inner_r,Pter_sizeZ*0.5,Pter_start_angle,Pter_spanning_angle);                //its size
	
	
	G4LogicalVolume* logicDelrin =
	new G4LogicalVolume(solidDelrin,          //its solid
											Delrin_mat,           //its material
											"Delrin");            //its name
	
	
	
	
	
	//################################################### Old case (Ga Container = Hole in a table)
	
	if(fGaSet==1){
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
	
		/*
		//###################################################
		// Table
		//##########################
		
		
		 G4ThreeVector posTable = G4ThreeVector(0, 0, -3.0*mm - DzTable/2.);
		 
		 
		 new G4PVPlacement(0,                     //no rotation
		 posTable,       //at (0,0,0)
		 logicTable,            //its logical volume
		 "Table",               //its name
		 logicWorld,            //its mother  volume
		 false,                 //no boolean operation
		 0,                     //copy number
		 checkOverlaps);        //overlaps checking
		 */
		 
	
	
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
	
	
		if (fCuDiam>=0){
	Z_FrontShield = DzCo + Pter_ZScan + FrontShield_sizeZ*0.5;
		} else{
			Z_FrontShield = Pter_ZScan + FrontShield_sizeZ*0.5;
		}
	
	/*
	if (fFilterFlag==0) { //if I do not want the filter, place it but make it thin and empty
		FrontShield_mat=world_mat;
	}*/
	
	G4ThreeVector posFilter = G4ThreeVector(fX0Scan, 0, Z_FrontShield);
	
	G4cout<<"GEOMETRY DEBUG - Z thickness of solidFrontShield= "<<FrontShield_sizeZ/mm<<", Z pos= "<<Z_FrontShield/mm<<G4endl;
	
	
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
	

	
	// place detector-Pter in world
		
		if (fCuDiam>=0){
			Pter_Posz = DzCo+ Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ*0.5;
		} else{
			Pter_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ*0.5;
		}

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
	// SiPm volume behind PTER
	//##########################
	
		
	G4ThreeVector posSiPm = G4ThreeVector(fX0Scan, 0, Pter_Posz + Pter_sizeZ/2. + DzSiPm/2.);

	
	new G4PVPlacement(0,                     //no rotation
										posSiPm,       //at (0,0,0)
										logicSiPm,            //its logical volume
										"SiPm",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
		
	
	
	
	//###################################################
	// PVC around P-Terphenyl
	//##########################
	
	
	new G4PVPlacement(0,                     //no rotation
										pos2,
										logicPVC,            //its logical volume
										"PVC",               //its name
										logicWorld,            //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	//###################################################
	// Delrin around P-Terphenyl
	//##########################
	

	
	
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
	
		
		
		
		//################ Top Prob Case
		
		
		G4Tubs* solidTopCase =
		new G4Tubs("solidTopCase",
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
		
		
		//################ Middle Case
		
		
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
		
		
		
		
		//################ Back Middle Case
		
		
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
		
		
		
		
		//################ End Case
		
		
		G4Tubs* solidEndCase =
		new G4Tubs("EndCase",
							 0.,
							 PVC_outer_r,
							 BackCaseThickness * 0.5,
							 SPhiTopCase,
							 DPhiTopCase);
	
	
	//###################################################
	// Middle Probe case
	//##########################
	
	
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
	// G4Union TopCase Probe
	//###################################################
	
		
	
	G4VSolid* TopCase=
	new G4SubtractionSolid ("TopCase",
													solidTopCase,
													solidSiPm,
													0,
													G4ThreeVector(0.,0.,TopCaseDepth*0.5-DzSiPm*0.5));
		
		
		
		G4LogicalVolume* logicTopCase =
		new G4LogicalVolume(TopCase,          //its solid
												TopCase_mat,           //its material
												"TopCase");            //its name
		
		
		
		G4double ProbeTopCase_Posz = Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
		G4ThreeVector posTopCase = G4ThreeVector(fX0Scan, 0, ProbeTopCase_Posz);
		
		G4RotationMatrix *rm = new G4RotationMatrix();
		rm->rotateY(180*deg);
		
		
		new G4PVPlacement(rm,                    // rotation
											posTopCase,            //at (0,0,0)
											logicTopCase,          //its logical volume
											"TopCase",             //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		

	//###################################################
	// End of Probe
	//##########################
	
	}
	
		
		
	}else if(fGaSet==2){      	//########## New case (Ga Container made with 3D Printer)

		
		//###################################################
		//Ga Source Container
		//##########################
		
		G4ThreeVector posContainerExtGa2 = G4ThreeVector(0, 0, -H_CylA*0.5);
		G4ThreeVector posExtGa2 = G4ThreeVector(0, 0, -H_CylB*0.5);
		
		
		G4VSolid* CylinderA =
		new G4Tubs("CylinderA",                       //its name
							 0.,
							 D_CylABCD*0.5,
							 H_CylA*0.5,
							 SPhiCyl,
							 DPhiCyl);     //its size
		
		
		G4VSolid* CylinderB =
		new G4Tubs("CylinderB",                       //its name
							 0.,
							 d_CylB*0.5,
							 H_CylB*0.5,
							 SPhiCyl,
							 DPhiCyl);     //its size
		
		
		G4VSolid* BackGaContainer=
		new G4SubtractionSolid ("BackGaContainer",      //GaContainer=GaCylinder-GaSource
														CylinderA,
														CylinderB,
														0,
														G4ThreeVector(0.,0.,H_CylA*0.5-H_CylB*0.5));
		
		
		

		
		G4LogicalVolume* logicSourceExtGa2 =
		new G4LogicalVolume(CylinderB,               //its solid
												SourceExtGa_mat,           //its material
												"SourceExtGa");            //its name
		
		
	
			new G4PVPlacement(0,                     //no rotation
												posExtGa2,       //at (0,0,0)
												logicSourceExtGa2,            //its logical volume
												"SourceExtGa",               //its name
												logicWorld,            //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
		
		
		G4VSolid* CylinderC =
		new G4Tubs("CylinderC",                       //its name
							 d_CylC*0.5,
							 D_CylABCD*0.5,
							 H_CylC*0.5,
							 SPhiCyl,
							 DPhiCyl);     //its size
		
		G4VSolid* CylinderD =
		new G4Tubs("CylinderD",                       //its name
							 d_CylD*0.5,
							 D_CylABCD*0.5,
							 H_CylD*0.5,
							 SPhiCyl,
							 DPhiCyl);     //its size
		
		G4VSolid* UnionCD
		= new G4UnionSolid("UnionCD",
											 CylinderC,
											 CylinderD,
											 0,
											 G4ThreeVector(0.,0.,(H_CylC+H_CylD)*0.5));
		//The G4ThreeVector(0.,0.,(MiddleCaseDepth+TopCaseDepth)*0.5) generates a translation of the second solid respect to the center of the first one. When I place the resulting solid I've to take in consideration as center of it the center of the first solid.
		
		
		
		G4VSolid* GaContainer
		= new G4UnionSolid("GaContainer",
											 BackGaContainer,
											 UnionCD,
											 0,
											 G4ThreeVector(0.,0.,(H_CylA + H_CylC)*0.5));
		
		
		G4LogicalVolume* logicGaContainer2 =
		new G4LogicalVolume(GaContainer,               //its solid
												GaContainer2_mat,           //its material
												"GaContainer");            //its name
		
		new G4PVPlacement(0,                     //no rotation
											posContainerExtGa2,       //at (0,0,0)
											logicGaContainer2,            //its logical volume
											"GaContainer",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);
		
		
		//###################################################
		//Probe Container
		//##########################
		
		G4VSolid* CylinderE =
		new G4Tubs("CylinderE",                       //its name
							 d_CylE*0.5,
							 D_CylEF*0.5,
							 H_CylE*0.5,
							 SPhiCyl,
							 DPhiCyl);     //its size
		
		G4VSolid* CylinderF =
		new G4Tubs("CylinderF",                       //its name
							 d_CylF*0.5,
							 D_CylEF*0.5,
							 H_CylF*0.5,
							 SPhiCyl,
							 DPhiCyl);     //its size
		
		G4VSolid* CylinderG =
		new G4Tubs("CylinderG",                       //its name
							 d_CylG*0.5,
							 D_CylG*0.5,
							 H_CylG*0.5,
							 SPhiCyl,
							 DPhiCyl);     //its size
		
		
		
		
	
		G4VSolid* UnionEF
		= new G4UnionSolid("UnionEF",
											 CylinderF,
											 CylinderE,
											 0,
											 G4ThreeVector(0.,0.,(H_CylE+H_CylF)*0.5));
		
		G4RotationMatrix* rm = new G4RotationMatrix();
		rm->rotateY(180.*deg);
		
		G4VSolid* ProbeContainer
		= new G4UnionSolid("ProbeContainer",
											 CylinderG,
											 UnionEF,
											 0,
											 G4ThreeVector(0.,0.,(H_CylG+H_CylF)*0.5));
		
		G4LogicalVolume* logicProbeContainer =
		new G4LogicalVolume(ProbeContainer,               //its solid
												ProbeContainer_mat,           //its material
												"ProbeContainer");            //its name
		
		
		G4ThreeVector ProbeContainerPos= G4ThreeVector(0,0,3.2*cm); //H_C+(H_G/2)+H_F+H_E
		
		new G4PVPlacement(rm,                     //no rotation
											ProbeContainerPos,       //at (0,0,0)
											logicProbeContainer,            //its logical volume
											"ProbeContainer",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);
		
		
		
		//###################################################
		//Absrober
		//##########################
		
		
		G4ThreeVector posAbs = G4ThreeVector(0, 0, ZCenterAbs);
		
		//G4cout<<"GEOMETRY DEBUG - Z thickness of solidShapeCo= "<<DzCo/mm<<", Z pos= "<<posCo.z()<<G4endl;
		
		G4Tubs* solidAbsorber =
		new G4Tubs("Absorber",                       //its name
							 RminAbs,
							 RmaxAbs,
							 0.5*DzAbs,
							 SPhiAbs,
							 DPhiAbs);     //its size
		
		if(fAbsorberMaterial==1){
			Absorber_mat = nist->FindOrBuildMaterial("G4_Cu");
		}else if(fAbsorberMaterial==2){
			Absorber_mat = nist->FindOrBuildMaterial("G4_Pb");
		}else if(fAbsorberMaterial==3){
			Absorber_mat = nist->FindOrBuildMaterial("MyAlu");
		}else if(fAbsorberMaterial==4){
			Absorber_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
		}
		
		G4cout<<"GEOMETRY DEBUG - absorber mat="<<Absorber_mat<<G4endl;
		
		G4LogicalVolume* logicAbsorber =
		new G4LogicalVolume(solidAbsorber,          //its solid
												Absorber_mat,           //its material
												"Absorber");            //its name
		
		if (fCuDiam>=0) {
			G4cout<<"GEOMETRY DEBUG - Copper collimator has been placed!!"<<G4endl;
			
			new G4PVPlacement(0,                     //no rotation
												posAbs,       //at (0,0,0)
												logicAbsorber,            //its logical volume
												"Absorber",               //its name
												logicWorld,            //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			//		G4Region* sorgente = new G4Region("SourceReg");
			logicAbsorber->SetRegion(sorgente);
			sorgente->AddRootLogicalVolume(logicAbsorber);
			
		}
		
		//################################################### END OF COPPER COLLIMATOR
		
	
		
		
		
		//###################################################
		//Dummy volume for scoring what exit from the source
		//##########################
		
		
		if (fCuDiam<0) {           // No Absorber
			zDummy2=DzDummy2*0.5;
		} else {                   // With Absorber (if fCuDiam>0 the absorber is drilled in the midle)
			zDummy2=DzDummy2*0.5+DzAbs*0.5+ ZCenterAbs;   //N.B. ZCenterAbs is the position of the center of the absorber
		}
		
		G4ThreeVector posDummy2 = G4ThreeVector(0, 0, zDummy2);
		
		
		
		new G4PVPlacement(0,                     //no rotation
											posDummy2,       //at (0,0,0)
											logicShapeDummy2,            //its logical volume
											"Dummy2",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		//		G4Region* sorgente = new G4Region("SourceReg");
		logicShapeDummy2->SetRegion(sorgente);
		sorgente->AddRootLogicalVolume(logicShapeDummy2);
		
		
		
		//###################################################
		//Electron Filter FrontShield
		//##########################
		
		
		
		//Z_FrontShield = DzAbs*0.5 + Pter_ZScan + DzDummy2 + FrontShield_sizeZ*0.5;
		Z_FrontShield = DzDummy2 + Pter_ZScan + FrontShield_sizeZ*0.5;
		
		G4ThreeVector posFilter = G4ThreeVector(fX0Scan, 0, Z_FrontShield);
		
		
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
		
		
		
		//###################################################
		// 	P-Terphenyl
		//##########################
		
		//Pter_Posz = DzAbs*0.5+ Pter_ZScan + DzDummy2 + FrontShield_sizeZ + Pter_sizeZ*0.5;
		Pter_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ*0.5;
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
		// SiPm volume behind PTER
		//##########################
		
		
		G4ThreeVector posSiPm = G4ThreeVector(fX0Scan, 0, Pter_Posz + Pter_sizeZ/2. + DzSiPm/2.);

		
		new G4PVPlacement(0,                     //no rotation
											posSiPm,       //at (0,0,0)
											logicSiPm,            //its logical volume
											"SiPm",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		//###################################################
		// Table
		//##########################
		
		/*
		G4ThreeVector posTable = G4ThreeVector(0, 0, -H_CylA - DzTable/2.);
		
		
		new G4PVPlacement(0,                     //no rotation
											posTable,       //at (0,0,0)
											logicTable,            //its logical volume
											"Table",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		*/
		
		//###################################################
		// PVC around P-Terphenyl
		//##########################
		
		
		new G4PVPlacement(0,                     //no rotation
											pos2,
											logicPVC,            //its logical volume
											"PVC",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking

		
		//###################################################
		// Delrin around P-Terphenyl
		//##########################
		
		
		
		
		new G4PVPlacement(0,                     //no rotation
											pos2,
											logicDelrin,            //its logical volume
											"Delrin",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		if(CaseDepth>0){

			
			//################ Top Prob Case
			
			
			G4Tubs* solidTopCase =
			new G4Tubs("solidTopCase",
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
			
			
			//################ Middle Case
			
			
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
			
			
			
			
			//################ Back Middle Case
			
			
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
			
			
			
			
			//################ End Case
			
			
			G4Tubs* solidEndCase =
			new G4Tubs("EndCase",
								 0.,
								 PVC_outer_r,
								 BackCaseThickness * 0.5,
								 SPhiTopCase,
								 DPhiTopCase);
			
			
			
			
		
		//###################################################
		// Middle Probe case
		//##########################
		
			//G4double ProbeMiddleCase_Posz = DzAbs*0.5+ Pter_ZScan + DzDummy2 + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth * 0.5;

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
		
		
		
		//G4double PlastiCase_Posz = DzAbs*0.5+ Pter_ZScan + DzDummy2 + DzDummy2 + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
		G4double PlastiCase_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
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
		
		//G4double HorsesShoe_Posz = DzAbs*0.5+ Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth*0.5;
		G4double HorsesShoe_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth*0.5;
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
		// G4Union TopCase Probe
		//###################################################
		
		
		
		G4VSolid* TopCase=
		new G4SubtractionSolid ("TopCase",
														solidTopCase,
														solidSiPm,
														0,
														G4ThreeVector(0.,0.,TopCaseDepth*0.5-DzSiPm*0.5));
		
		
		
		G4LogicalVolume* logicTopCase =
		new G4LogicalVolume(TopCase,          //its solid
												TopCase_mat,           //its material
												"TopCase");            //its name
		
		
		//G4double ProbeTopCase_Posz = DzAbs*0.5+ Pter_ZScan + DzDummy2 + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
		G4double ProbeTopCase_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
		G4ThreeVector posTopCase = G4ThreeVector(fX0Scan, 0, ProbeTopCase_Posz);
		
		G4RotationMatrix *rm1 = new G4RotationMatrix();
		rm1->rotateY(180*deg);
		
		
		new G4PVPlacement(rm1,                    // rotation
											posTopCase,            //at (0,0,0)
											logicTopCase,          //its logical volume
											"TopCase",             //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		
		
		//###################################################
		// End of Probe
		//##########################
		
		}else if (CaseDepth<0){
		
			
		G4double AluCaseDepth=fabs(CaseDepth);
		/*
		G4double Alu_Posz = DzAbs*0.5+ Pter_ZScan + DzDummy2 + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth-1.15)*0.5*mm;
		G4ThreeVector posAlu = G4ThreeVector(0, 0, Alu_Posz);
		
		G4double BackAlu_Posz = DzAbs*0.5+ Pter_ZScan + DzDummy2 + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth-1.15)*mm + 1.15*0.5*mm;
		G4ThreeVector posBackAlu = G4ThreeVector(0, 0, BackAlu_Posz);
		*/
		G4double Alu_Posz = DzDummy2 + Pter_ZScan  + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth-1.15)*0.5*mm;
		G4ThreeVector posAlu = G4ThreeVector(0, 0, Alu_Posz);
			
		G4double BackAlu_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth-1.15)*mm + 1.15*0.5*mm;
		G4ThreeVector posBackAlu = G4ThreeVector(0, 0, BackAlu_Posz);
			
			
		G4Tubs* solidAlu =
		new G4Tubs("AluCase",
							 9.7*0.5*mm,
							 12.*0.5*mm,
							 (AluCaseDepth-1.15)*0.5*mm,
							 SPhiTopCase,
							 DPhiTopCase);
		
		
		G4Tubs* solidBackAlu =
		new G4Tubs("BackAluCase",
							 0.,
							 12.*0.5*mm,
							 1.15*0.5*mm,
							 SPhiTopCase,
							 DPhiTopCase);
		
		G4LogicalVolume* logicsolidAlu =
		new G4LogicalVolume(solidAlu,          //its solid
												FrontShield_mat,           //its material
												"AluCase");            //its name
		
		G4LogicalVolume* logicsolidBackAlu =
		new G4LogicalVolume(solidBackAlu,          //its solid
												FrontShield_mat,           //its material
												"BackAluCase");            //its name
		
		new G4PVPlacement(0,                    // rotation
											posAlu,            //at (0,0,0)
											logicsolidAlu,          //its logical volume
											"AluCase",             //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		new G4PVPlacement(0,                    // rotation
											posBackAlu,            //at (0,0,0)
											logicsolidBackAlu,          //its logical volume
											"BackAluCase",             //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		}
		
		
	}else if(fGaSet==3){ // end of GaSet 2
		//########## New case (Ga Container made with PVC)
		
		
		//###################################################
		//Ga Source Container
		//##########################
		
		G4ThreeVector posContainerExtGa2 = G4ThreeVector(0, 0, -H_CylA3*0.5);
		G4ThreeVector posExtGa3 = G4ThreeVector(0, 0, -H_CylB3*0.5-DzDummy3/2); //centro la sorgente in -0.4
		
		
		G4VSolid* CylinderA =
		new G4Tubs("CylinderA",                       //its name
							 0.,
							 D_CylABCD3*0.5,
							 H_CylA3*0.5,
							 SPhiCyl3,
							 DPhiCyl3);     //its size
		
		
		G4VSolid* CylinderB =
		new G4Tubs("CylinderB",                       //its name
							 0.,
							 d_CylB3*0.5,
							 H_CylB3*0.5,
							 SPhiCyl3,
							 DPhiCyl3);     //its size
		
		G4VSolid* SourceGaDOTATOC =
		new G4Tubs("SourceGa-DOTATOC",
							 0.,
							 d_CylB3*0.5,
							 (fSourceThickness*mm)*0.5,
							 SPhiCyl3,
							 DPhiCyl3);
		
		
		G4VSolid* BackGaContainer=
		new G4SubtractionSolid ("BackGaContainer",      //GaContainer=GaCylinder-GaSource
														CylinderA,
														CylinderB,
														0,
														G4ThreeVector(0.,0.,H_CylA3*0.5-H_CylB3*0.5));
		
		
		
		
		
		G4LogicalVolume* logicSourceExtGa2 =
		new G4LogicalVolume(SourceGaDOTATOC,               //its solid
												SourceExtGa_mat,           //its material
												"SourceExtGa");            //its name
		
		
		
		new G4PVPlacement(0,                     //no rotation
											posExtGa3,             //at (0,0,0)
											logicSourceExtGa2,            //its logical volume
											"SourceExtGa",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		G4VSolid* CylinderC =
		new G4Tubs("CylinderC",                       //its name
							 d_CylC3*0.5,
							 D_CylABCD3*0.5,
							 H_CylC3*0.5,
							 SPhiCyl3,
							 DPhiCyl3);     //its size
		
		G4VSolid* CylinderD =
		new G4Tubs("CylinderD",                       //its name
							 d_CylD3*0.5,
							 D_CylABCD3*0.5,
							 H_CylD3*0.5,
							 SPhiCyl3,
							 DPhiCyl3);     //its size
		
		G4VSolid* UnionCD
		= new G4UnionSolid("UnionCD",
											 CylinderC,
											 CylinderD,
											 0,
											 G4ThreeVector(0.,0.,(H_CylC3+H_CylD3)*0.5));
		//The G4ThreeVector(0.,0.,(MiddleCaseDepth+TopCaseDepth)*0.5) generates a translation of the second solid respect to the center of the first one. When I place the resulting solid I've to take in consideration as center of it the center of the first solid.
		
		
		
		G4VSolid* GaContainer
		= new G4UnionSolid("GaContainer",
											 BackGaContainer,
											 UnionCD,
											 0,
											 G4ThreeVector(0.,0.,(H_CylA3 + H_CylC3)*0.5));
		
		
		G4LogicalVolume* logicGaContainer2 =
		new G4LogicalVolume(GaContainer,               //its solid
												GaContainer3_mat,           //its material
												"GaContainer");            //its name
		
		new G4PVPlacement(0,                     //no rotation
											posContainerExtGa2,       //at (0,0,0)
											logicGaContainer2,            //its logical volume
											"GaContainer",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);
		
		
		//###################################################
		//Probe Container
		//##########################
		
		G4VSolid* CylinderE =
		new G4Tubs("CylinderE",                       //its name
							 d_CylE3*0.5,
							 D_CylEF3*0.5,
							 H_CylE3*0.5,
							 SPhiCyl3,
							 DPhiCyl3);     //its size
		
		G4VSolid* CylinderF =
		new G4Tubs("CylinderF",                       //its name
							 d_CylF3*0.5,
							 D_CylEF3*0.5,
							 H_CylF3*0.5,
							 SPhiCyl3,
							 DPhiCyl3);     //its size
		
		G4VSolid* CylinderG =
		new G4Tubs("CylinderG",                       //its name
							 d_CylG3*0.5,
							 D_CylG3*0.5,
							 H_CylG3*0.5,
							 SPhiCyl3,
							 DPhiCyl3);     //its size
		
		
		
		
		
		G4VSolid* UnionEF
		= new G4UnionSolid("UnionEF",
											 CylinderF,
											 CylinderE,
											 0,
											 G4ThreeVector(0.,0.,(H_CylE3+H_CylF3)*0.5));
		
		G4RotationMatrix* rm = new G4RotationMatrix();
		rm->rotateY(180.*deg);
		
		G4VSolid* ProbeContainer
		= new G4UnionSolid("ProbeContainer",
											 CylinderG,
											 UnionEF,
											 0,
											 G4ThreeVector(0.,0.,(H_CylG3+H_CylF3)*0.5));
		
		G4LogicalVolume* logicProbeContainer =
		new G4LogicalVolume(ProbeContainer,               //its solid
												ProbeContainer_mat,           //its material
												"ProbeContainer");            //its name
		
		
		G4ThreeVector ProbeContainerPos= G4ThreeVector(0,0,3.31*cm);  //H_C+(H_G/2)+H_F+H_E
		
		new G4PVPlacement(rm,                     //no rotation
											ProbeContainerPos,       //at (0,0,0)
											logicProbeContainer,            //its logical volume
											"ProbeContainer",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);
		
		
		
		//###################################################
		//Absrober
		//##########################
		
		G4ThreeVector posAbs = G4ThreeVector(0, 0, ZCenterAbs);
		
		//G4cout<<"GEOMETRY DEBUG - Z thickness of solidShapeCo= "<<DzCo/mm<<", Z pos= "<<posCo.z()<<G4endl;
		
		G4Tubs* solidAbsorber =
		new G4Tubs("Absorber",                       //its name
							 RminAbs,
							 RmaxAbs,
							 0.5*DzAbs,
							 SPhiAbs,
							 DPhiAbs);     //its size
		
		if(fAbsorberMaterial==1){
			Absorber_mat = nist->FindOrBuildMaterial("G4_Cu");
		}else if(fAbsorberMaterial==2){
			Absorber_mat = nist->FindOrBuildMaterial("G4_Pb");
		}else if(fAbsorberMaterial==3){
			Absorber_mat = nist->FindOrBuildMaterial("MyAlu");
		}else if(fAbsorberMaterial==4){
			Absorber_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
		}
		
		G4cout<<"GEOMETRY DEBUG - absorber mat="<<Absorber_mat<<G4endl;
		
		G4LogicalVolume* logicAbsorber =
		new G4LogicalVolume(solidAbsorber,          //its solid
												Absorber_mat,           //its material
												"Absorber");            //its name
		
		if (fCuDiam>=0) {
			G4cout<<"GEOMETRY DEBUG - Copper collimator has been placed!!"<<G4endl;
			
			new G4PVPlacement(0,                     //no rotation
												posAbs,       //at (0,0,0)
												logicAbsorber,            //its logical volume
												"Absorber",               //its name
												logicWorld,            //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			//		G4Region* sorgente = new G4Region("SourceReg");
			logicAbsorber->SetRegion(sorgente);
			sorgente->AddRootLogicalVolume(logicAbsorber);
			
		}
		
		//################################################### END OF COPPER COLLIMATOR
		
		
		
		
		
		//###################################################
		//Dummy volume for scoring what exit from the source
		//##########################
		
		if (fCuDiam<0) {           // No Absorber
			zDummy2=DzDummy2*0.5;
		} else {                   // With Absorber (if fCuDiam>0 the absorber is drilled in the midle)
			zDummy2=DzDummy2*0.5+DzAbs*0.5+ ZCenterAbs;   //N.B. Pter_ZScan is the position of the center of the absorber
		}
		
		G4ThreeVector posDummy2 = G4ThreeVector(0, 0, zDummy2);
		
		
		
		new G4PVPlacement(0,                     //no rotation
											posDummy2,       //at (0,0,0)
											logicShapeDummy2,            //its logical volume
											"Dummy2",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		//		G4Region* sorgente = new G4Region("SourceReg");
		logicShapeDummy2->SetRegion(sorgente);
		sorgente->AddRootLogicalVolume(logicShapeDummy2);
		
		zDummy3=-DzDummy3*0.5;

		
		G4ThreeVector posDummy3 = G4ThreeVector(0, 0, zDummy3);
		
		
		
		new G4PVPlacement(0,                     //no rotation
											posDummy3,       //at (0,0,0)
											logicShapeDummy3,            //its logical volume
											"Dummy3",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		//		G4Region* sorgente = new G4Region("SourceReg");
		logicShapeDummy3->SetRegion(sorgente);
		sorgente->AddRootLogicalVolume(logicShapeDummy3);
		
		
		
		//###################################################
		//Electron Filter FrontShield
		//##########################
		
		
		
		Z_FrontShield = DzDummy2 + Pter_ZScan + FrontShield_sizeZ*0.5;
		
		G4ThreeVector posFilter = G4ThreeVector(fX0Scan, 0, Z_FrontShield);
		
		
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
		
		
		
		//###################################################
		// 	P-Terphenyl
		//##########################
		
		
		Pter_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ*0.5;
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
		// SiPm volume behind PTER
		//##########################
		
		
		G4ThreeVector posSiPm = G4ThreeVector(fX0Scan, 0, Pter_Posz + Pter_sizeZ/2. + DzSiPm/2.);
		
		
		new G4PVPlacement(0,                     //no rotation
											posSiPm,       //at (0,0,0)
											logicSiPm,            //its logical volume
											"SiPm",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		//###################################################
		// Table
		//##########################
		
		/*
		 G4ThreeVector posTable = G4ThreeVector(0, 0, -H_CylA - DzTable/2.);
		 
		 
		 new G4PVPlacement(0,                     //no rotation
		 posTable,       //at (0,0,0)
		 logicTable,            //its logical volume
		 "Table",               //its name
		 logicWorld,            //its mother  volume
		 false,                 //no boolean operation
		 0,                     //copy number
		 checkOverlaps);        //overlaps checking
		 
		 
		 */
		
		//###################################################
		// PVC around P-Terphenyl
		//##########################
		
		
		new G4PVPlacement(0,                     //no rotation
											pos2,
											logicPVC,            //its logical volume
											"PVC",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		//###################################################
		// Delrin around P-Terphenyl
		//##########################
		
		
		
		
		new G4PVPlacement(0,                     //no rotation
											pos2,
											logicDelrin,            //its logical volume
											"Delrin",               //its name
											logicWorld,            //its mother  volume
											false,                 //no boolean operation
											0,                     //copy number
											checkOverlaps);        //overlaps checking
		
		
		if(CaseDepth>0){
			
			
			//################ Top Prob Case
			
			
			G4Tubs* solidTopCase =
			new G4Tubs("solidTopCase",
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
			
			
			//################ Middle Case
			
			
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
			
			
			
			
			//################ Back Middle Case
			
			
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
			
			
			
			
			//################ End Case
			
			
			G4Tubs* solidEndCase =
			new G4Tubs("EndCase",
								 0.,
								 PVC_outer_r,
								 BackCaseThickness * 0.5,
								 SPhiTopCase,
								 DPhiTopCase);
			
			
			
			
			
			//###################################################
			// Middle Probe case
			//##########################
			
			
			G4double ProbeMiddleCase_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth * 0.5;
			G4ThreeVector posMiddleCase = G4ThreeVector(fX0Scan, 0, ProbeMiddleCase_Posz);
			
			
			new G4PVPlacement(0,                     //no rotation
												posMiddleCase,       //at (0,0,0)
												logicMiddleCase,            //its logical volume
												"MiddleCase",               //its name
												logicWorld,            //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			
			
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
			
			
			
			G4double PlastiCase_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
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
			
			G4double HorsesShoe_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth + MiddleCaseDepth*0.5;
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
			// G4Union TopCase Probe
			//###################################################
			
			
			
			G4VSolid* TopCase=
			new G4SubtractionSolid ("TopCase",
															solidTopCase,
															solidSiPm,
															0,
															G4ThreeVector(0.,0.,TopCaseDepth*0.5-DzSiPm*0.5));
			
			
			
			G4LogicalVolume* logicTopCase =
			new G4LogicalVolume(TopCase,          //its solid
													TopCase_mat,           //its material
													"TopCase");            //its name
			
			
			
			G4double ProbeTopCase_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + TopCaseDepth*0.5;
			G4ThreeVector posTopCase = G4ThreeVector(fX0Scan, 0, ProbeTopCase_Posz);
			
			G4RotationMatrix *rm1 = new G4RotationMatrix();
			rm1->rotateY(180*deg);
			
			
			new G4PVPlacement(rm1,                    // rotation
												posTopCase,            //at (0,0,0)
												logicTopCase,          //its logical volume
												"TopCase",             //its name
												logicWorld,            //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			
			
			
			//###################################################
			// End of Probe
			//##########################
			
		}else if (CaseDepth<0){
			
			
			G4double AluCaseDepth=fabs(CaseDepth);
			
			G4double Alu_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth-1.15)*0.5*mm;  //1.15 BackCase Thickness
			G4ThreeVector posAlu = G4ThreeVector(0, 0, Alu_Posz);
			
			G4double BackAlu_Posz = DzDummy2 + Pter_ZScan + FrontShield_sizeZ + Pter_sizeZ + (AluCaseDepth-1.15)*mm + 1.15*0.5*mm;
			G4ThreeVector posBackAlu = G4ThreeVector(0, 0, BackAlu_Posz);
			
			G4Tubs* solidAlu =
			new G4Tubs("AluCase",
								 9.7*0.5*mm,
								 12.*0.5*mm,
								 (AluCaseDepth-1.15)*0.5*mm,
								 SPhiTopCase,
								 DPhiTopCase);
			
			
			G4Tubs* solidBackAlu =
			new G4Tubs("BackAluCase",
								 0.,
								 12.*0.5*mm,
								 1.15*0.5*mm,
								 SPhiTopCase,
								 DPhiTopCase);
			
			G4LogicalVolume* logicsolidAlu =
			new G4LogicalVolume(solidAlu,          //its solid
													FrontShield_mat,           //its material
													"AluCase");            //its name
			
			G4LogicalVolume* logicsolidBackAlu =
			new G4LogicalVolume(solidBackAlu,          //its solid
													FrontShield_mat,           //its material
													"BackAluCase");            //its name
			
			new G4PVPlacement(0,                    // rotation
												posAlu,            //at (0,0,0)
												logicsolidAlu,          //its logical volume
												"AluCase",             //its name
												logicWorld,            //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			
			new G4PVPlacement(0,                    // rotation
												posBackAlu,            //at (0,0,0)
												logicsolidBackAlu,          //its logical volume
												"BackAluCase",             //its name
												logicWorld,            //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
		} // end of GaSet3
		
	}
	
	
	return physWorld;
	
	
	
	
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Voglio Pter_inner_r da terminale


