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
// $Id: exampleB1.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "B1DetectorConstruction.hh"
#include "B1ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
//#include "QBBC.hh"
#include "B1PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include "G4StepLimiter.hh"
#include "G4UserLimits.hh"
#include "G4StepLimiterPhysics.hh"

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
	G4bool VisFlag=true;
	
	// Detect interactive mode (if no arguments) and define UI session
	G4UIExecutive* ui = 0;
	/*
	if ( argc == 13 ) {  //was argc==1, 7 to see geom using input parameters, 8 once added sensorchoice
		ui = new G4UIExecutive(argc, argv);
	}
	*/
	
	G4double x0Scan=0, ZValue=2, AbsorberDiam=0, TBRvalue=1,PterDiameter=6,PterThickness=5,SourceDiameter=5.25,SourceThickness=5, AbsorberThickness=1.,ProbeCaseDepth=40, ProbeCaseLateralThickness=3, ProbeCaseBackThickness=5 , HSLateralThickness=1, HSBackThickness=2;
	G4int SourceChoice=1, AbsorberMaterial=1;
	
	G4String fileName ="";

	
	for(int i=1;i<argc;i++)
		if(argv[i][0] =='-')
		{
			G4String option(argv[i]);
			G4cout<<"option: "<<i<<" "<<option<<G4endl;
			if(option.compare("-AbsD")==0)
			{
				AbsorberDiam=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-AbsT")==0)
			{
				AbsorberThickness=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-AbsMat")==0)
			{
				AbsorberMaterial=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Z")==0)
			{
				ZValue=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-TBR")==0)
			{
				TBRvalue=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Source")==0)
			{
				SourceChoice=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-X")==0)
			{
				x0Scan=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-PterD")==0)
			{
				PterDiameter=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-PterT")==0)
			{
				PterThickness=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-SourceD")==0)
			{
				SourceDiameter=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-SourceT")==0)
			{
				SourceThickness=strtod (argv[++i], NULL);;
			}else if(option.compare("-CaseDepth")==0)         //Probe Case Z_Lenght
			{
				ProbeCaseDepth=strtod (argv[++i], NULL);;
			}else if(option.compare("-CaseLT")==0)            //Lateral Probe Case Thickness
			{
				ProbeCaseLateralThickness=strtod (argv[++i], NULL);;
			}else if(option.compare("-CaseBT")==0)            //Back Probe Case Thickness
			{
				ProbeCaseBackThickness=strtod (argv[++i], NULL);;
			}else if(option.compare("-HSLT")==0)            //Lateral Probe Case Thickness
			{
				HSLateralThickness=strtod (argv[++i], NULL);;
			}else if(option.compare("-HSBT")==0)            //Back Probe Case Thickness
			{
				HSBackThickness=strtod (argv[++i], NULL);;
			}
			
		}
		else
		{
			fileName = argv[i]; //se ho trovato una macro (senza il "-" davanti) significa che NON voglio l'interattivo
			VisFlag=false;
		}
	
	if ( VisFlag ) { //Prepare for vis
		ui = new G4UIExecutive(argc, argv);
	}
	
	
	G4int SourceSelect=SourceChoice;
	//if (SourceSelect==1|| SourceSelect==2) SrSourceFlag=1; //if it is a Sr source... tell to DetCons
	
	G4String FileNamePrim="Primaries";
	G4String OutFileName="PTERmc";
	G4String FileNameCommonPart;
	
	G4String MaterialiAssorbitore[4]= {"Cu","Pb","Al","PVC"};
	
	FileNameCommonPart.append("_PDiam" + std::to_string((G4int)PterDiameter)+"_PDz" + std::to_string((G4int)PterThickness));
	
	if (AbsorberDiam>=0) FileNameCommonPart.append("_AbsDz" + std::to_string((G4int)(1000*AbsorberThickness))+"_AbsHole" + std::to_string((G4int)AbsorberDiam) +"_AbsMat" + MaterialiAssorbitore[AbsorberMaterial-1]);
	else FileNameCommonPart.append("_NoAbs");
	
	
	if (ProbeCaseDepth>0) FileNameCommonPart.append("_CaseDepth" + std::to_string((G4int)(ProbeCaseDepth))+"_CaseLT" + std::to_string((G4int)ProbeCaseLateralThickness) + "_CaseBT" + std::to_string((G4int)(ProbeCaseBackThickness))+"_HSLT" + std::to_string((G4int)HSLateralThickness)+"_HSBT" + std::to_string((G4int)HSBackThickness) );
	
	

	FileNameCommonPart.append("_X"+ std::to_string((G4int)x0Scan));
	FileNameCommonPart.append("_Z"+ std::to_string((G4int)ZValue));
	if (SourceSelect==1) FileNameCommonPart.append("_PSr");
	if (SourceSelect==2) FileNameCommonPart.append("_ExtSr");
	if (SourceSelect==3) FileNameCommonPart.append("_ExtY");
<<<<<<< HEAD
	if (SourceSelect==4) FileNameCommonPart.append("_ExtGa_Diam" + std::to_string((G4int)(10*SourceDiameter)) + "_Dz" + std::to_string((G4int)(10*SourceThickness)));
=======
	if (SourceSelect==4) FileNameCommonPart.append("_ExtGa_Diam" + std::to_string((G4int)(1000*SourceDiameter)) + "_Dz" + std::to_string((G4int)SourceThickness));
>>>>>>> ee6c8af87b18526ef1aac4034fa99d63e5ff9e2d
	if (SourceSelect==5) FileNameCommonPart.append("_Sphere511");

	FileNameCommonPart.append("");
	
	if (VisFlag) FileNameCommonPart.append("TEST"); //if it was a TEST run under vis
	
	FileNamePrim.append(FileNameCommonPart);
	OutFileName.append(FileNameCommonPart);

	/*
	if (CuDiam>=0){
		FileNameCommonPart="X"+ std::to_string((G4int)x0Scan) + "_Z" + std::to_string((G4int)(100*ZValue)) + "_CuD" + std::to_string((G4int)CuDiam) + "_Fil" + std::to_string((G4int)FilterFlag)  + "_TBR" + std::to_string((G4int)(10*TBRvalue))  ;
	}
	else	{
		FileNameCommonPart="PrimariesX" + std::to_string((G4int)x0Scan) + "_Z" + std::to_string((G4int)(100*ZValue)) + "_NoCuD"  + "_Fil" + std::to_string((G4int)FilterFlag)  + "_TBR" + std::to_string((G4int)(10*TBRvalue))  ;
	}
	
	FileNamePrim.append(FileNameCommonPart);

	
	if (fCuDiam>=0){
		FileNameCommonPart= "X"+  std::to_string((G4int)fX0Scan) + "_Z" + std::to_string((G4int)(100*fZValue)) + "_CuD" + std::to_string((G4int)fCuDiam) + "_TBR" + std::to_string((G4int)(10*fTBR))    );
	}
	else {
		fileName= fileNameBase + "X"+  std::to_string((G4int)fX0Scan) + "_Z" + std::to_string((G4int)(100*fZValue)) + "_NOCuD" + "_Fil" + std::to_string((G4int)fFilterFlag) + "_TBR" + std::to_string((G4int)(10*fTBR));
	}
	
	
	
	if (SourceSelect==1) FileNamePrim.append("_PSr");
	if (SourceSelect==2) FileNamePrim.append("_ExtSr");
	if (SourceSelect==3) FileNamePrim.append("_ExtY");
	if (SourceSelect==4) FileNamePrim.append("_ExtGa");
	*/
	/*
	if (SensorChoice==1) FileNamePrim.append("_011");
	if (SensorChoice==2) FileNamePrim.append("_115");
	if (SensorChoice==3) FileNamePrim.append("_60035");
*/
	FileNamePrim.append(+ ".dat");
	std::ofstream primFile(FileNamePrim, std::ios::out);
	
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	
	// Construct the default run manager
	//
	//#ifdef G4MULTITHREAD
	//  G4MTRunManager* runManager = new G4MTRunManager;
	//#else
	G4RunManager* runManager = new G4RunManager;
	//#endif
	
	// Set mandatory initialization classes
	// Detector construction
	runManager->SetUserInitialization(new B1DetectorConstruction(x0Scan, ZValue, AbsorberDiam, SourceSelect, AbsorberMaterial,PterDiameter,PterThickness,SourceDiameter,SourceThickness,AbsorberThickness,ProbeCaseDepth,ProbeCaseLateralThickness,ProbeCaseBackThickness,HSLateralThickness,HSBackThickness)); //DetectorConstruction needs to know if it is a SrSource to place the right geometry
	
	// Physics list
	//G4VModularPhysicsList* physicsList = new QBBC;
	//physicsList->SetVerboseLevel(1);
	
	//  runManager->SetUserInitialization(new B1PhysicsList);
	
	B1PhysicsList* physicsList=new B1PhysicsList;
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	runManager->SetUserInitialization(physicsList);
	
	// User action initialization
	//	runManager->SetUserInitialization(new B1ActionInitialization(x0Scan, ZValue, CuDiam, FilterFlag, primFile, TBRvalue,SourceSelect, SourceSelect));
	runManager->SetUserInitialization(new B1ActionInitialization(x0Scan, ZValue, AbsorberDiam,  primFile, TBRvalue, SourceSelect, AbsorberMaterial, SourceDiameter, SourceThickness, OutFileName));
	
	// Initialize visualization
	//
	G4VisManager* visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
	
	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	
	// Process macro or start UI session
	//
	
	if ( ! ui ) {
		// batch mode
		G4String command = "/control/execute ";
//		G4String fileName = argv[13];
		UImanager->ApplyCommand(command+fileName);
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;
	}
	
	delete visManager;
	delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
