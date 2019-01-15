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
// $Id: B1RunAction.hh 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1RunAction.hh
/// \brief Definition of the B1RunAction class

#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4AccumulableManager.hh"
#include "globals.hh"
#include <vector>

#include <iostream>

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class B1RunAction : public G4UserRunAction
{
public:
	B1RunAction(G4double, G4double, G4double, G4double, G4int, G4int, G4String);
	virtual ~B1RunAction();
	
	// virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void   EndOfRunAction(const G4Run*);
	
	void AddEdep (G4double edep);
	void AddEdepSiPM (G4double EdepSiPM);
	void AddEdkin (G4double edkin);
	
	// void AddEdepPhot (G4double edepPhot);
	// void AddEdepEl (G4double edepEl);
	
	
	std::vector<G4double>& GetRunEnPter() {return RunVPterEn; }
	std::vector<G4double>& GetRunEnPterPrim() {return RunVPterPrimEn; }
	std::vector<G4double>& GetRunPartPterPrim() {return RunVPterPrimPart; }
	std::vector<G4float>& GetRunEnPterTime() {return RunVPterTime; }
	std::vector<G4double>& GetRunXPter() {return RunVPterX; }
	std::vector<G4double>& GetRunYPter() {return RunVPterY; }
	std::vector<G4double>& GetRunZPter() {return RunVPterZ; }
	std::vector<G4double>& GetRunPartPter() {return RunVPterPart; }
	
	std::vector<G4double>& GetRunEnExit() {return RunVEnExit; }
	std::vector<G4double>& GetRunXExit() {return RunVXExit; }
	std::vector<G4double>& GetRunYExit() {return RunVYExit; }
	std::vector<G4double>& GetRunZExit() {return RunVZExit; }
	std::vector<G4double>& GetRunCosXExit() {return RunVCosXExit; }
	std::vector<G4double>& GetRunCosYExit() {return RunVCosYExit; }
	std::vector<G4double>& GetRunCosZExit() {return RunVCosZExit; }
	std::vector<G4double>& GetRunPartExit() {return RunVPartExit; }
	std::vector<G4double>& GetRunParentIDExit() {return RunVParentIDExit; }
	
	std::vector<G4int>& GetRunExitProcess() {return RunExitProcess; }
	
	std::vector<G4double>& GetRunCosX() {return RunVCosX; }
	std::vector<G4double>& GetRunCosY() {return RunVCosY; }
	std::vector<G4double>& GetRunCosZ() {return RunVCosZ; }
	std::vector<G4double>& GetRunEnGen() {return RunVEnGen; }
	std::vector<G4double>& GetRunEnPart() {return RunVEnPart; }
	std::vector<G4double>& GetRunIsotopeGen() {return RunVIsotopeGen; }
	
	std::vector<G4double>& GetRunEAbsComp() {return RunVEAbsComp; }
	std::vector<G4double>& GetRunEAbsSiPMComp() {return RunVEAbsSiPMComp; }


	std::vector<G4double>& GetRunEnPre() {return RunVPrePterEn; }
	std::vector<G4double>& GetRunPart() {return RunVPrePterPart; }
	
	
	std::vector<G4double>& GetRunAnnihX() {return RunVAnnihX; }
	std::vector<G4double>& GetRunAnnihY() {return RunVAnnihY; }
	std::vector<G4double>& GetRunAnnihZ() {return RunVAnnihZ; }
	
	std::vector<G4double>& GetRunPreAbsEn() {return RunVPreAbsEn; }
	std::vector<G4double>& GetRunPartPreAbs() {return RunVPartPreAbs; }
	std::vector<G4double>& GetRunPartPostAbs() {return RunVPartPostAbs; }


	std::vector<G4double>& GetRunAnnihT() {return RunVAnnihT; }

	

	
	
	std::vector<G4double>& GetRunPixNo() {return RunVPixNo; }
//	std::vector<G4double>& GetRunPixEneDep() {return RunVPixEneDep; }
	std::vector<G4double>& GetRunPixXpos() {return RunVPixXpos; }
	std::vector<G4double>& GetRunPixYpos() {return RunVPixYpos; }
	
	G4int GetEventNumber() {return nbEventInRun;}
	
	void SetMotherIsotope(G4double miso) {fMotherIsotope=miso;}
	G4double GetMotherIsotope() {return fMotherIsotope;}
	
	void SetMotherEnergy(G4double mene) {fMotherEnergy=mene;}
	G4double GetMotherEnergy() {return fMotherEnergy;}
	
	void SetMotherTime(G4double mtime) {fMotherTime=mtime;}
	G4float GetMotherTime() {return fMotherTime;}
	
	
	/*
	 std::vector<G4double>& GetSourceX() {return RunVSourceX; }
	 std::vector<G4double>& GetSourceY() {return RunVSourceY; }
	 std::vector<G4double>& GetSourceZ() {return RunVSourceZ; }
	 std::vector<G4double>& GetSourceEne() {return RunVSourceEne; }
	 */
	
	
private:
	G4Accumulable<G4double> fEdep;
//	G4Accumulable<G4double> fEdepSiPM;
	G4Accumulable<G4double> fEdep2;
	G4Accumulable <G4double> fEdkin;
	
	G4double fX0Scan;
	G4double fZValue;
	G4double fAbsHoleDiam;
	G4int nbEventInRun;
	G4double fTBR;
	G4int fSourceSelect;
	
	G4int fMotherIsotope=-10;
	G4int fSensorChoice;

	G4double fMotherEnergy=-10;
	G4float fMotherTime=0;
	
	//G4Accumulable <G4double> fEdepPhot;
	//G4Accumulable<G4double> fEdepEl;
	
	/////////////////
	// Histogramming
	//
	void CreateHistogram();
	void WriteHistogram();
	
	std::vector<G4double> RunVPterEn;
	std::vector<G4double>	RunVPterPrimEn;
	std::vector<G4double> RunVPterPrimPart;
	std::vector<G4float>	RunVPterTime;
	std::vector<G4double> RunVPterX;
	std::vector<G4double> RunVPterY;
	std::vector<G4double> RunVPterZ;
	std::vector<G4double> RunVPterPart;
	
	std::vector<G4double> RunVPrePterEn;
	std::vector<G4double> RunVPrePterPart;
	
	
	std::vector<G4double> RunVAnnihX;
	std::vector<G4double> RunVAnnihY;
	std::vector<G4double> RunVAnnihZ;
	
	std::vector<G4double> RunVPreAbsEn;
	std::vector<G4double> RunVPartPreAbs;
	std::vector<G4double> RunVPartPostAbs;
	
	std::vector<G4double> RunVAnnihT;

	
	std::vector<G4double> RunVPixNo;
//	std::vector<G4double> RunVPixEneDep;
	std::vector<G4double> RunVPixXpos;
	std::vector<G4double> RunVPixYpos;
	
	std::vector<G4double> RunVCosX;
	std::vector<G4double> RunVCosY;
	std::vector<G4double> RunVCosZ;
	
	std::vector<G4double> RunVEnGen;
	std::vector<G4double> RunVEnPart;
	std::vector<G4double> RunVIsotopeGen;
	
	
	std::vector<G4double> RunVEnExit;
	std::vector<G4double> RunVXExit;
	std::vector<G4double> RunVYExit;
	std::vector<G4double> RunVZExit;
	std::vector<G4double> RunVCosXExit;
	std::vector<G4double> RunVCosYExit;
	std::vector<G4double> RunVCosZExit;
	
	std::vector<G4double> RunVPartExit;
	std::vector<G4double> RunVParentIDExit;
	
	std::vector<G4int> RunExitProcess;
	
	std::vector<G4double> RunVEAbsComp;
	std::vector<G4double> RunVEAbsSiPMComp;

	
	/*
	 std::vector<G4double> RunVSourceX;
	 std::vector<G4double> RunVSourceY;
	 std::vector<G4double> RunVSourceZ;
	 std::vector<G4double> RunVSourceEne;
	 */
	G4String fOutFileName;
	
};

#endif

