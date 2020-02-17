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
#include "MyExceptionHandler.hh"


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QBBC.hh"
#include "FTFP_BERT_HP.hh"
#include "B1PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include "G4StepLimiter.hh"
#include "G4UserLimits.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4OpticalPhysics.hh"

#include "G4ScoringManager.hh"

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
	G4bool VisFlag=true;
	
	// Detect interactive mode (if no arguments) and define UI session
	G4UIExecutive* ui = 0;
	
	G4double x0Scan=0., ZValue=2., AbsorberHoleDiam=-1., TBRvalue=1.,PterDiameter=6.,PterThickness=5.,SourceDiameter=10.,SourceThickness=7., AbsorberThickness=1.,ProbeCaseDepth=-155., ProbeCaseLateralThickness=1.25, ProbeCaseBackThickness=20. , HSLateralThickness=1., HSBackThickness=2., AbsCenter=2.75;
	G4int SourceChoice=1, AbsorberMaterial=1, HousingCase=3, GaSetting=1,ApparatusMat=1,PosAbsorber=1;
	G4bool ScintFlag=0;
	
	G4String fileName ="";
	G4String FileNameLabel="";
	
	for(int i=1;i<argc;i++)
		if(argv[i][0] =='-')
		{
			G4String option(argv[i]);
			G4cout<<"option: "<<i<<" "<<option<<G4endl;
			if(option.compare("-AbsD")==0)
			{
				AbsorberHoleDiam=strtod (argv[++i], NULL);;
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
			}else if(option.compare("-HSLT")==0)            //Lateral Thickness Horseshoe
			{
				HSLateralThickness=strtod (argv[++i], NULL);;
			}else if(option.compare("-HSBT")==0)            //Back Thickness Horseshoe
			{
				HSBackThickness=strtod (argv[++i], NULL);;
			}else if(option.compare("-HSMat")==0)
			{
				HousingCase=strtod (argv[++i], NULL);;
			}else if(option.compare("-Scint")==0)
			{
				ScintFlag= argv[++i];;
			}else if(option.compare("-Label")==0)
			{
				FileNameLabel= argv[++i];;
			}else if(option.compare("-GaSet")==0)
			{
				GaSetting= strtod (argv[++i], NULL);;
			}else if(option.compare("-AppMat")==0)
			{
				ApparatusMat= strtod (argv[++i], NULL);;
			}else if(option.compare("-PosAbs")==0)
			{
				PosAbsorber= strtod (argv[++i], NULL);;
			}else if(option.compare("-ZAbs")==0)
			{
				AbsCenter= strtod (argv[++i], NULL);;
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
	G4int GaSet=GaSetting;
	
	G4String OutFileName="PTERmc";
	G4String FileNameCommonPart;
	
	G4String MaterialiAssorbitore[4]= {"Cu","Pb","Al","PVC"};
	
	// ###### PTER
	FileNameCommonPart.append("_PDiam" + std::to_string((G4int)PterDiameter)+"_PDz" + std::to_string((G4int)PterThickness));
	
	// ###### X and Z
	FileNameCommonPart.append("_X"+ std::to_string((G4int)(10*x0Scan)));
	FileNameCommonPart.append("_Z"+ std::to_string((G4int)(10*ZValue)));
	
	// ###### ABSORBER
	if (AbsorberHoleDiam<0) FileNameCommonPart.append("_NoAbs");
	else {
		if (GaSet==1)  FileNameCommonPart.append("_AbsDz" + std::to_string((G4int)(1000*AbsorberThickness))+"_AbsHoleD" + std::to_string((G4int)AbsorberHoleDiam) +"_AbsMat" + MaterialiAssorbitore[AbsorberMaterial-1]);
		if (GaSet==2 || GaSet==3) FileNameCommonPart.append("_PosAbs"+std::to_string((G4int)(PosAbsorber))+"_AbsT" + std::to_string((G4int)(100*AbsorberThickness))+"_AbsHoleD" + std::to_string((G4int)AbsorberHoleDiam) +"_AbsMat" + MaterialiAssorbitore[AbsorberMaterial-1]);
	}
	
	// ###### IF LAPAROSCOPIC CASE
	if (ProbeCaseDepth>0) FileNameCommonPart.append("_CaseDepth" + std::to_string((G4int)(ProbeCaseDepth))+"_CaseLT" + std::to_string((G4int)ProbeCaseLateralThickness) + "_CaseBT" + std::to_string((G4int)(ProbeCaseBackThickness))+"_HSLT" + std::to_string((G4int)HSLateralThickness)+"_HSBT" + std::to_string((G4int)HSBackThickness)+"_HSMat" + std::to_string(HousingCase) );
	
	// ###### SOURCES
	if (SourceSelect==1) FileNameCommonPart.append("_PSr");
	if (SourceSelect==2) FileNameCommonPart.append("_ExtSr");
	if (SourceSelect==3) FileNameCommonPart.append("_ExtY");
	if (SourceSelect==4) {
		FileNameCommonPart.append("_ExtGa");
		if (GaSet== 1) FileNameCommonPart.append("_Diam" + std::to_string((G4int)(10*SourceDiameter)) + "_Dz" + std::to_string((G4int)(10*SourceThickness)) + "_Set1");
		if (GaSet== 2 || GaSet==3) FileNameCommonPart.append("_GaSet"+std::to_string((G4int)(GaSet))+"_AluCaseT" + std::to_string((G4int)(fabs(ProbeCaseDepth))) + "_AppMat" + std::to_string((G4int)(ApparatusMat)));
	}
	if (SourceSelect==5) FileNameCommonPart.append("_Sphere511");
	if (SourceSelect==6) FileNameCommonPart.append("_FlatEle");
	if (SourceSelect==7) FileNameCommonPart.append("_FlatGamma");
	if (SourceSelect==8) FileNameCommonPart.append("_VolCu67_Diam" + std::to_string((G4int)(10*SourceDiameter)) + "_Dz" + std::to_string((G4int)(10*SourceThickness)));

	if (SourceSelect==9) FileNameCommonPart.append("_VolF18_Diam" + std::to_string((G4int)(10*SourceDiameter)) + "_Dz" + std::to_string((G4int)(10*SourceThickness)));

	if (SourceSelect<0) FileNameCommonPart.append("_Z" + std::to_string(int(-SourceSelect/100)) +"_A" + std::to_string(int(-SourceSelect%100)) + "_Diam" + std::to_string((G4int)(10*SourceDiameter)) + "_Dz" + std::to_string((G4int)(10*SourceThickness)));

	
	// ####### MISCELLANEUS
	if (ScintFlag) FileNameCommonPart.append("_Scint");
	if (FileNameLabel!="") FileNameCommonPart.append("_" + FileNameLabel);
	if (VisFlag) FileNameCommonPart.append("TEST"); //if it was a TEST run under vis
	
	OutFileName.append(FileNameCommonPart);
	
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4long seed = time(NULL);
	if (VisFlag) seed=12345; //If vis was requested same always the same seed to have reproducibility
	G4Random::setTheSeed(seed);
	// Construct the default run manager
	//
	//#ifdef G4MULTITHREAD
	//  G4MTRunManager* runManager = new G4MTRunManager;
	//#else
	
#if 0
	G4MTRunManager* runManager = new G4MTRunManager;
	//	runManager->SetNumberOfThreads( G4Threading::G4GetNumberOfCores() );
	runManager->SetNumberOfThreads( 6 );
#endif
	
	G4RunManager* runManager = new G4RunManager;
	
	// Set mandatory initialization classes
	// Detector construction
	runManager->SetUserInitialization(new B1DetectorConstruction(x0Scan, ZValue, AbsorberHoleDiam, SourceSelect, AbsorberMaterial,PterDiameter,PterThickness,SourceDiameter,SourceThickness,AbsorberThickness,ProbeCaseDepth,ProbeCaseLateralThickness,ProbeCaseBackThickness,HSLateralThickness,HSBackThickness, HousingCase, ScintFlag, GaSet, ApparatusMat, PosAbsorber, AbsCenter));
	
	// Physics list
	//G4VModularPhysicsList* physicsList = new QBBC;
	//physicsList->SetVerboseLevel(1);
	
	//  runManager->SetUserInitialization(new B1PhysicsList);
	
	//	G4VModularPhysicsList* physicsList = new QBBC;
	//	physicsList->RegisterPhysics(new G4OpticalPhysics());
	
	B1PhysicsList* physicsList=new B1PhysicsList(ScintFlag);
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	
	runManager->SetUserInitialization(physicsList);
	
	// User action initialization
	//	runManager->SetUserInitialization(new B1ActionInitialization(x0Scan, ZValue, AbsHoleDiam, FilterFlag, primFile, TBRvalue,SourceSelect, SourceSelect));
	runManager->SetUserInitialization(new B1ActionInitialization(x0Scan, ZValue, AbsorberHoleDiam, TBRvalue, SourceSelect, AbsorberMaterial, SourceDiameter, SourceThickness, OutFileName, GaSetting, ProbeCaseDepth));
	
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
