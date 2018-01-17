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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction(G4double x0, G4double ZValue, G4double CuDiam, G4int FilterFlag, G4bool SrSourceFlag, G4int SensorChoice)
: G4VUserDetectorConstruction(),
fScoringVolume(0), fX0Scan(x0), fZValue(ZValue), fCuDiam(CuDiam), fFilterFlag(FilterFlag), fSrSourceFlag(SrSourceFlag), fSensorChoice(SensorChoice)
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
	G4double world_sizeXY = 0.2*m;
	G4double world_sizeZ  = 0.2*m;
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
	
	
	//###################################################
	
	
	
	//materials definition
	
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

	
	density = 4.000*g/cm3;
	G4Material* Resin = new G4Material (name="Resin", density, ncomponents=3);
	Resin->AddElement (elH, natoms=30);
	Resin->AddElement (elC, natoms=20);
	Resin->AddElement (elO, natoms=2);
	
	G4double densityAlu = 2.600*g/cm3;
//	G4NistManager* man = G4NistManager::Instance();
	G4NistManager::Instance()->BuildMaterialWithNewDensity("MyAlu","G4_Al",densityAlu);

	/* ref from TestEm7
	G4NistManager::Instance()->
	BuildMaterialWithNewDensity("Water_1.05","G4_WATER",1.05*g/cm3);
	*/

	//	G4bool DOTASource=FALSE;
	G4bool SrSource=fSrSourceFlag;
	
	
	//###################################################
	// AGAR AGAR Source - AgarAgar should be C14 H24 O9
	//##########################
	
	G4double Agardensity = 1.030*g/cm3;
	G4Material* AgarAgar = new G4Material (name="AgarAgar", Agardensity, ncomponents=3);
	AgarAgar->AddElement (elH, natoms=24);
	AgarAgar->AddElement (elC, natoms=14);
	AgarAgar->AddElement (elO, natoms=9);
	
	G4double RminSourceDOTA = 0.*mm;
	G4double RmaxSourceDOTA = 10.5*mm;
	G4double DzSourceDOTA= 4.5*mm;
	G4double SPhiSourceDOTA = 0.*deg;
	G4double DPhiSourceDOTA = 360.*deg;
	//	G4Material* SourceDOTA_mat = nist->FindOrBuildMaterial("G4_Al");
	G4Material* SourceDOTA_mat = AgarAgar;
	G4ThreeVector posSourceDOTA = G4ThreeVector(0, 0, -DzSourceDOTA*0.5);
	
	
	G4Tubs* solidSourceDOTA =
	new G4Tubs("SourceDOTA",                       //its name
			   RminSourceDOTA,
			   RmaxSourceDOTA,
			   0.5*DzSourceDOTA,
			   SPhiSourceDOTA,
			   DPhiSourceDOTA);     //its size
	
	G4LogicalVolume* logicSourceDOTA =
	new G4LogicalVolume(solidSourceDOTA,          //its solid
						SourceDOTA_mat,           //its material
						"SourceDOTA");            //its name
	
	// G4VPhysicalVolume* physSr=
	
	if(!SrSource) { //I place the DOTA source if I am not asking for Sr source
		G4cout<<"GEOMETRY DEBUG - Z thickness of solidSourceDOTA= "<<DzSourceDOTA/mm<<", Z pos= "<<-DzSourceDOTA*0.5<<G4endl;

		G4cout<<"GEOMETRY DEBUG - DOTATOC Source has been placed!!"<<G4endl;

		new G4PVPlacement(0,                     //no rotation
						  posSourceDOTA,       //at (0,0,0)
						  logicSourceDOTA,            //its logical volume
						  "SourceDOTA",               //its name
						  logicWorld,            //its mother  volume
						  false,                 //no boolean operation
						  0,                     //copy number
						  checkOverlaps);        //overlaps checking
		

	//G4Region* sorgente = new G4Region("SourceReg");
	logicSourceDOTA->SetRegion(sorgente);
	sorgente->AddRootLogicalVolume(logicSourceDOTA);
	}
	//###################################################
	
	
	//###################################################
	// ABS carrier around DOTA Source - ABS should be C15 H17 N
	//##########################
	
	G4double ABSdensity = 0.7*g/cm3;
	G4Material* ABS = new G4Material (name="ABS", ABSdensity, ncomponents=3);
	ABS->AddElement (elH, natoms=17);
	ABS->AddElement (elC, natoms=15);
	ABS->AddElement (elN, natoms=1);
	
	
	G4double RminABSaround = 10.5*mm;
	G4double RmaxABSaround = 12.5*mm;
	G4double DzABSaround= DzSourceDOTA;
	G4double SPhiABSaround = 0.*deg;
	G4double DPhiABSaround = 360.*deg;
	//	G4Material* Source_mat = nist->FindOrBuildMaterial("G4_Al");
	G4Material* ABSaround_mat = ABS;
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
	
	// G4VPhysicalVolume* physSr=
		if(!SrSource) {  //I place the ABS carrier of the DOTA source if I am not asking for Sr source
			G4cout<<"GEOMETRY DEBUG - Z thickness of solidABSaround= "<<DzABSaround/mm<<", Z pos= "<<-DzABSaround*0.5<<G4endl;

			G4cout<<"GEOMETRY DEBUG - DOTATOC Source has been placed!!"<<G4endl;

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
	
	//###################################################
	
	//###################################################
	// ABS carrier behind DOTA Source - ABS should be C15 H17 N
	//##########################
	
	G4double RminABSbehind = 0.*mm;
	G4double RmaxABSbehind = 12.5*mm;
	G4double DzABSbehind= 3*mm;
	G4double SPhiABSbehind = 0.*deg;
	G4double DPhiABSbehind = 360.*deg;
	G4Material* ABSbehind_mat = ABS;
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
	
	// G4VPhysicalVolume* physSr=
		if(!SrSource) { //I place the ABS carrier of the DOTA source if I am not asking for Sr source
			G4cout<<"GEOMETRY DEBUG - Z thickness of solidABSbehind= "<<DzABSbehind/mm<<", Z pos= "<<-DzABSbehind*0.5- DzABSaround<<G4endl;

			G4cout<<"GEOMETRY DEBUG - DOTATOC Source has been placed!!"<<G4endl;

			new G4PVPlacement(0,                     //no rotation
					  posABSbehind,       //at (0,0,0)
					  logicABSbehind,            //its logical volume
					  "ABSbehind",               //its name
					  logicWorld,            //its mother  volume
					  false,                 //no boolean operation
					  0,                     //copy number
					  checkOverlaps);        //overlaps checking
	
		}
	logicABSbehind->SetRegion(ABSRegion);
	ABSRegion->AddRootLogicalVolume(logicABSbehind);
	
	
	
	
	//###################################################
	// Sr90 lab Source -
	//##########################
	
	
	G4double RminSourceSR = 0.*mm;
	G4double RmaxSourceSR = 12.5*mm; //12.5mm for PG source, 8mm for RM
	G4double DzSourceSR= 3*mm;
	G4double SPhiSourceSR = 0.*deg;
	G4double DPhiSourceSR = 360.*deg;
	G4Material* SourceSR_mat = nist->FindOrBuildMaterial("MyAlu"); //G4_Al
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
	
	// G4VPhysicalVolume* physSr=
	
	if(SrSource) { //If i requested the Sr source
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
	//###################################################

	
	//###################################################
	//Electron Filter Resin
	//##########################

	G4double Resin_sizeX=0*mm;
	G4double Resin_sizeY = 0*mm;
	G4double Resin_sizeZ  = 0.*mm;
	G4double Z_resin= 0*mm;
//	cmos_ZScan=fZValue + Resin_sizeZ + Cmos_sizeZ*0.5; //modified on 2017.11.21 by collamaf - Z distance does not take into account Cu thickness! is always from source top

	if (fSensorChoice==2) fFilterFlag=1; //Sensor 2 is always with filter
	//	Resin=nist->FindOrBuildMaterial("G4_AIR");
	if (fFilterFlag==1) {
		
		Resin_sizeX = 2.69*mm;
		Resin_sizeY = 3.58*mm;
		Resin_sizeZ  = 0.520*mm;
		if (fCuDiam>0) Z_resin= fZValue + Resin_sizeZ*0.5; //corrected 2017.12.12 by collamaf - no need to add 1mm of copper since distance is always from source top
		else Z_resin= fZValue + Resin_sizeZ*0.5;
		G4ThreeVector pos4 = G4ThreeVector(fX0Scan, 0,Z_resin);
		
		G4cout<<"GEOMETRY DEBUG - Z thickness of solidResin= "<<Resin_sizeZ/mm<<", Z pos= "<<Z_resin/mm<<G4endl;

		G4Box* solidResin =
		new G4Box("Resin",                       //its name
				  0.5*Resin_sizeX,0.5*Resin_sizeY,0.5*Resin_sizeZ);     //its size

		G4LogicalVolume* logicResin =
		new G4LogicalVolume(solidResin,          //its solid
							Resin,           //its material
							"Resin");            //its name

		new G4PVPlacement(0,                     //no rotation
						  pos4,       //at (0,0,0)
						  logicResin,            //its logical volume
						  "Resin",               //its name
						  logicWorld,            //its mother  volume
						  false,                 //no boolean operation
						  0,                     //copy number
						  checkOverlaps);        //overlaps checking
		//	G4Region* filtro = new G4Region("ResinReg");
		logicResin->SetRegion(filtro);
		filtro->AddRootLogicalVolume(logicResin);
		
	}
	
	//###################################################
	
	
	
	//###################################################
	//Copper Collimator
	//##########################
	
	
	G4double RminCo = fabs(fCuDiam)/2.;
	G4double RmaxCo = 18.*mm;
	G4double DzCo= 1.*mm;
	G4double SPhiCo = 0.*deg;
	G4double DPhiCo = 360.*deg;
	G4double zCo=DzCo*0.5;
	G4Material* shapeCo_mat = nist->FindOrBuildMaterial("G4_Cu");
	G4ThreeVector posCo = G4ThreeVector(0, 0, zCo);
	
	G4cout<<"GEOMETRY DEBUG - Z thickness of solidShapeCo= "<<DzCo/mm<<", Z pos= "<<zCo<<G4endl;
	
	
	G4Tubs* solidShapeCo =
	new G4Tubs("CuCollimator",                       //its name
			   RminCo,
			   RmaxCo,
			   0.5*DzCo,
			   SPhiCo,
			   DPhiCo);     //its size
	
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

	
	//###################################################
	

	
	//###################################################
	//Dummy volume for scoring
	//##########################
	
	
	G4double RminDummy = 0;
	G4double RmaxDummy = 18.*mm;
	if (fCuDiam>=0) RmaxDummy =RmaxCo;
	else RmaxDummy=RmaxSourceSR;
	G4double DzDummy= 1.e-5*mm;
	G4double SPhiDummy = 0.*deg;
	G4double DPhiDummy = 360.*deg;
	G4double zDummy=DzDummy*0.5;
	
	if (fCuDiam>0) {
		zDummy=DzDummy*0.5+DzCo;
	} else {
		zDummy=DzDummy*0.5;
	}
	
	G4Material* shapeDummy_mat = nist->FindOrBuildMaterial("G4_AIR");
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
	
	
	//	if (fCuDiam>=0) {
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
	
	
	//	}
	
	
	//###################################################
	
	
	
	
	
	//###################################################
	// CMOS Si sensoor PIXEL
	//
	
	//pixel dimension
	
	G4int ScaleFactor=10; //set to 1 for full simulation, 10 for quick view
	
	//defaults geom values are for MTV011 Sensor (1)
	G4double PixelSize=5.6*um;
	G4double PixelThickness=4.5*um;
	G4double gapX =0.01*um*ScaleFactor;
	G4double gapY =0.01*um*ScaleFactor;
	G4int noX = 480;
	G4int noY = 640;
	G4int copyNo=0;
	G4Material* pix_mat = nist->FindOrBuildMaterial("G4_Si");
	
	if (fSensorChoice==2) {
		PixelSize=1.75*um;
		PixelThickness=129*um;
		noX=488;
		noY=648;
	}
	
	G4double pixX =PixelSize*ScaleFactor;
	G4double pixY =PixelSize*ScaleFactor;
	G4double pixZ =PixelThickness;
	noX/=ScaleFactor;
	noY/=ScaleFactor;
	
	
	//CMOS dimension
	
	G4double Cmos_sizeX = (pixX+gapX)*noX-gapX;
	G4double Cmos_sizeY = (pixY+gapY)*noY-gapY;
	G4double Cmos_sizeZ  = pixZ;
	
	
	//G4double cmos_Z= Z_resin+ Resin_sizeZ*0.5 + Cmos_sizeZ*0.5;
	//	G4Material* Cmos_mat = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* Cmos_mat = nist->FindOrBuildMaterial("G4_Si");
	//	G4double cmos_ZScan=DzSr*0.5+DzCo+5.5*mm;
	
	G4double cmos_ZScan;
	/*
	if(fCuDiam>=0)
	cmos_ZScan=fZValue + Resin_sizeZ + Cmos_sizeZ*0.5 + 1*mm;
	else cmos_ZScan=fZValue + Resin_sizeZ + Cmos_sizeZ*0.5;
	 */
	cmos_ZScan=fZValue + Resin_sizeZ + Cmos_sizeZ*0.5; //modified on 2017.11.21 by collamaf - Z distance does not take into account Cu thickness! is always from source top to possible resin
	G4ThreeVector pos2 = G4ThreeVector(fX0Scan, 0, cmos_ZScan);
	
	G4cout<<"GEOMETRY DEBUG - Z thickness of solidCmos= "<<Cmos_sizeZ/mm<<", Z pos= "<<cmos_ZScan/mm<<G4endl;
	
	//CMOS
	
	G4Box* solidCmos =
	new G4Box("CMOS",                       //its name
			  0.5*Cmos_sizeX,0.5*Cmos_sizeY,0.5*Cmos_sizeZ);     //its size
	
	G4LogicalVolume* logicCmos =
	new G4LogicalVolume(solidCmos,          //its solid
						Cmos_mat,           //its material
						"CMOS");            //its name
	
	
	//pixel
	
	G4Box* solidPix =
	new G4Box("Pix",                       //its name
			  0.5*pixX,0.5*pixY,0.5*pixZ);     //its size
	
	G4LogicalVolume* logicPix =
	new G4LogicalVolume(solidPix,          //its solid
						pix_mat,           //its material
						"Pix");            //its name
	
	
	//G4Region* cmosreg = new G4Region("CMOSReg");
	logicCmos->SetRegion(cmosreg);
	cmosreg->AddRootLogicalVolume(logicCmos);
	
	logicPix->SetRegion(cmosreg);
	cmosreg->AddRootLogicalVolume(logicPix);
	
	//placement of the pixel in CMOS
	
	for (G4int ix = 1; ix <= noX ; ix++){
		for (G4int iy = 1; iy <= noY ; iy++){
			G4ThreeVector posPixX = G4ThreeVector((-0.5*Cmos_sizeX+ix*(pixX+gapX)-0.5*pixX-gapX),(-0.5*Cmos_sizeY+iy*(pixY+gapY)-0.5*pixY-gapY) ,0);
			copyNo++;
			
			new G4PVPlacement(0,                     //no rotation
							  posPixX,       //at (0,0,0)
							  logicPix,            //its logical volume
							  "CMOS",               //its name
							  logicCmos,            //its mother  volume
							  false,                 //no boolean operation
							  copyNo,                     //copy number
							  0);        //overlaps checking
		}
	}
	
	// place detector-CMOS in world
	
	new G4PVPlacement(0,                     //no rotation
					  pos2,
					  logicCmos,            //its logical volume
					  "CMOS",               //its name
					  logicWorld,            //its mother  volume
					  false,                 //no boolean operation
					  0,                     //copy number
					  checkOverlaps);        //overlaps checking
	
	
	//
	// CMOS-carrier PVC
	//
	
	G4double carrier_sizeX = 50.*mm;
	G4double carrier_sizeY = 70.*mm;
	G4double carrier_sizeZ  = 2.*mm;
	//G4double carrier_Z= Z_resin+ Resin_sizeZ*0.5 + Cmos_sizeZ+0.5*carrier_sizeZ;
	G4double carrier_Z=cmos_ZScan +0.5*carrier_sizeZ + Cmos_sizeZ;
	G4Material* carrier_mat = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4ThreeVector posCarrier = G4ThreeVector(fX0Scan, 0, carrier_Z);
	
	
	G4cout<<"GEOMETRY DEBUG - Z thickness of solidCarrier= "<<carrier_sizeZ/mm<<", Z pos= "<<carrier_Z/mm<<G4endl;
	
	
	G4Box* solidCarrier =
	new G4Box("Carrier",                       //its name
			  0.5*carrier_sizeX,0.5*carrier_sizeY,0.5*carrier_sizeZ);     //its size
	
	G4LogicalVolume*  logicCarrier=
	new G4LogicalVolume(solidCarrier,          //its solid
						carrier_mat,           //its material
						"Carrier");            //its name
	
	
	new G4PVPlacement(0,                     //no rotation
					  posCarrier,       //at (0,0,0)
					  logicCarrier,            //its logical volume
					  "Carrier",               //its name
					  logicWorld,            //its mother  volume
					  false,                 //no boolean operation
					  0,                     //copy number
					  checkOverlaps);        //overlaps checking
	
	
	//G4Region* carrierreg = new G4Region("CarrierReg");
	logicCarrier->SetRegion(carrierreg);
	carrierreg->AddRootLogicalVolume(logicCarrier);
	
	
	
	
	
	// Set scoring volume
	//
	
	//Pixelated CMOS
	
	fScoringVolume = logicPix;
	
	//Solid Si CMOS
	
	//fScoringVolume = logicCmos;
	
	//
	//always return the physical World
	//
	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
