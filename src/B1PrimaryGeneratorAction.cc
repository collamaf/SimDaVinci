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
// $Id: B1PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "G4IonTable.hh"
#include "G4ChargedGeantino.hh"

#include "B1RunAction.hh"
#include "B1Analysis.hh"


#include "G4Event.hh"

#include <iostream>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using std::ofstream;
using std::ios;
using std::endl;


B1PrimaryGeneratorAction::B1PrimaryGeneratorAction(B1EventAction* eventAction, G4double TBR, G4int SourceSelect, G4double SourceDiameter, G4double SourceThickness, G4int GaSetting)
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0) ,
evtPrimAction(eventAction), fTBR(TBR), fSourceSelect(SourceSelect), fSourceDiameter(SourceDiameter), fSourceThickness(SourceThickness),fGaSet(GaSetting)

{
	G4int n_particle = 1;
	fParticleGun  = new G4ParticleGun(n_particle);
	G4bool fPointLike=true;
	G4bool fExtended=false;
	G4bool fSTB=false; 
	G4bool fGa=false;


	if (fSourceSelect==1) {  //pointlike Sr
		fPointLike=true;
		fExtended=false;
		fSTB=false;
	} else if (fSourceSelect==2 ) { //extended Sr
		fPointLike=false;
		fExtended=true;
		fSTB=false;
	} else if (fSourceSelect==6 ) { //extended flat Ele
		fPointLike=false;
		fExtended=false;
		fSTB=false;
		FlatEle=true;
	} else if (fSourceSelect==7 ) { //extended flat Gamma
		fPointLike=false;
		fExtended=false;
		fSTB=false;
		FlatGamma=true;
	} else if (fSourceSelect==3) { //Ext Y
		fPointLike=false;
		fExtended=false;
		fSTB=true;
	} else if (fSourceSelect==4) {  //Ext Ga
		fPointLike=false;
		fExtended=false;
		fSTB=false;
		fGa=true;
	}
	
	if (fSTB) {
	fRadiusInt=3*mm;
	fDZInt=1*mm;
	fRadiusExt=10.48*mm; //10.48 per Rosa, 6.65 per PG
	fDZExt=4.4*mm;   //4.4 per Rosa, 5.5 per PG
	} else if (fPointLike) {
		fRadiusInt=0*mm;
		fDZInt=0*mm;
		fRadiusExt=0*mm;
		fDZExt=0*mm;
	} else if (fExtended) {
		fRadiusInt=8*mm;  //8 for RM, 10.5mm PG source
		fDZInt=0*mm;
		fRadiusExt=8*mm;
		fDZExt=0*mm;
	} else if (FlatEle) {
		fRadiusInt=6*mm;
		fDZInt=0*mm;
		fRadiusExt=6*mm;
		fDZExt=0*mm;
	} else if (fGa) {
		fRadiusInt=fSourceDiameter/2.*mm;
		fDZInt=0*mm;
		fRadiusExt=fSourceDiameter/2.*mm;
		fDZExt=fSourceThickness*mm;
	}
	
	
	ofstream SourceFile;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("geantino");
	
	fParticleGun->SetParticleDefinition(particle);
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
	delete fParticleGun;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries (G4Event* anEvent)
{
//	fSourceSelect=4;
	//Stronzium
	G4int Z = 38, A = 90;
	if (fSourceSelect==3) Z=39; //If I need Y instead of Sr
	if (fSourceSelect==4 ) {
		Z=31;
		A=68;
	} else if (fSourceSelect==6) { //Flat Ele
		FlatEle=true;
	} else if (fSourceSelect==7) { //Flat Gamma
		FlatGamma=true;
	}
	G4double ionCharge   = 0.*eplus;
	G4double excitEnergy = 0.*keV;
	
	
//Retrieve position of MiddleCase that represents the center of the whole probe to use it as offset for Sphere source
	//G4double Sphere_ZOffset=G4PhysicalVolumeStore::GetInstance()->GetVolume("MiddleCase")->GetTranslation().z();
//	G4cout<<"AAAAA "<<Sphere_ZOffset<<G4endl;
	
	G4ParticleDefinition* ion
	= G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
	fParticleGun->SetParticleDefinition(ion);
	fParticleGun->SetParticleCharge(ionCharge);
	
	G4double VolA=CLHEP::pi*fDZExt*(fRadiusExt*fRadiusExt-fRadiusInt*fRadiusInt);
	G4double VolB=CLHEP::pi*fRadiusInt*fRadiusInt*fDZInt;
	G4double VolC=CLHEP::pi*fRadiusInt*fRadiusInt*(fDZExt-fDZInt);
	G4double denominatore=VolA+VolB*fTBR+VolC;
	G4double ProbA=VolA/denominatore;
	G4double ProbB=VolB*fTBR/denominatore;
	G4double ProbC=VolC/denominatore;
	
	G4double zSource=0;
	G4double zSourceOffset=1e-6*mm; //to avoid generating particles at the very boundary of source!
	
	//	if (fRadiusExt==fRadiusInt) { //se ho un solo raggio ignoro il TBR e faccio la pasticca di sorgente
	if (fSourceSelect==1||fSourceSelect==2 || fSourceSelect==6) { //se ho una delle due pasticche di Sr ignoro il TBR e faccio la pasticca di sorgente
		fRadiusMax=fRadiusInt;
		fRadiusMin=0*mm;
		zSource = -zSourceOffset;
	} else if (fSourceSelect==3) { // Extended source with TBR
		G4double random=G4UniformRand();
		if (random<=ProbA) {  //faccio il cilindretto cavo esterno al centro (VolA)
			fRadiusMax=fRadiusExt;
			fRadiusMin=fRadiusInt;
			fZ=fDZExt;
			zSource = -G4UniformRand()*fZ-zSourceOffset;
		} else if (random>ProbA && random<=ProbA+ProbB) {    //faccio il cilindretto attivo al centro (VolB) SEGNALE!!!!
			fRadiusMax=fRadiusInt;
			fRadiusMin=0*mm;
			fZ=fDZInt;
			zSource = -G4UniformRand()*fZ-zSourceOffset;
		} else if (random>ProbA+ProbB) {     //faccio il cilindretto dietro a quello attivo al centro (VolC)
			fRadiusMax=fRadiusInt;
			fRadiusMin=0*mm;
			fZ=fDZExt-fDZInt;
			zSource = -G4UniformRand()*fZ-fDZInt-zSourceOffset;
		}
	} else if (fSourceSelect==4 && fGaSet==1) {
		fRadiusMax=fRadiusInt;
		fRadiusMin=0*mm;
		fZ=fDZExt;
		zSource = -G4UniformRand()*fZ-zSourceOffset;
	} else if (fSourceSelect==4 && fGaSet==2) {
		fRadiusMax=fRadiusInt;
		fRadiusMin=0*mm;
		fZ=7*mm;
		zSource = -G4UniformRand()*fZ-zSourceOffset;
	}else if (fSourceSelect==4 && fGaSet==3) {
		fRadiusMax=fRadiusInt;
		fRadiusMin=0*mm;
		fZ=fDZExt;
		zSource = -G4UniformRand()*fZ-zSourceOffset-(7.5*mm-fDZExt);  //source's height is less than container height
	}
	

	//Inserisco sorgente 5 (sfera di fotoni)
	
	G4double Sphere_Theta=G4UniformRand()*CLHEP::pi*2.;
	G4double Sphere_Phi=acos(2*G4UniformRand()-1);
	G4double Sphere_U=cos(Sphere_Phi);
	G4double Sphere_X=sqrt(1-Sphere_U*Sphere_U)*cos(Sphere_Theta);
	G4double Sphere_Y=sqrt(1-Sphere_U*Sphere_U)*sin(Sphere_Theta);
	G4double Sphere_Z=Sphere_U;
	G4double Sphere_Radius=5*cm;
	
	fParticleGun->SetParticleEnergy(0*MeV); //SetParticleEnergy uses kinetic energy
	
	G4double rho = sqrt(fRadiusMin*fRadiusMin + G4UniformRand()*(fRadiusMax*fRadiusMax-fRadiusMin*fRadiusMin));   //fixed square problem by collamaf with internal radius!
	G4double alpha = G4UniformRand()*CLHEP::pi*2.;

	
	G4double PhotDir_cosTheta = 2*G4UniformRand() - 1., PhotDir_phi = CLHEP::pi*2.*G4UniformRand();
	G4double PhotDir_sinTheta = std::sqrt(1. - PhotDir_cosTheta*PhotDir_cosTheta);
	G4double PhotDir_ux = PhotDir_sinTheta*std::cos(PhotDir_phi),
	PhotDir_uy = PhotDir_sinTheta*std::sin(PhotDir_phi),
	PhotDir_uz = PhotDir_cosTheta;
	
	
	
	G4double Source_X=rho*cos(alpha);
	G4double Source_Y=rho*sin(alpha);
	G4double Source_Z=zSource;

	
//	G4ThreeVector position = G4ThreeVector(rho*cos(alpha), rho*sin(alpha), zSource);

//	G4ThreeVector position = G4ThreeVector(Sphere_X, Sphere_Y, Sphere_Z);
/*
	evtPrimAction->SetSourceCosX(0);
	evtPrimAction->SetSourceCosY(0);
	evtPrimAction->SetSourceCosZ(0);
*/

	
	G4double Sphere_ZOffset=0;
	
	G4double xDirection=0;
	G4double yDirection=0;
	G4double zDirection=0;
	
	if (fSourceSelect==5) {
		Source_X=Sphere_Radius*Sphere_X;
		Source_Y=Sphere_Radius*Sphere_Y;
		if (G4PhysicalVolumeStore::GetInstance()->GetVolume("MiddleCase")!=NULL)	Sphere_ZOffset=G4PhysicalVolumeStore::GetInstance()->GetVolume("MiddleCase")->GetTranslation().z();
		Source_Z=Sphere_Radius*Sphere_Z+Sphere_ZOffset;
		G4ParticleDefinition* fotone = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
		fParticleGun->SetParticleDefinition(fotone);
		fParticleGun->SetParticleEnergy(0.511*MeV);
		xDirection=PhotDir_ux;
		yDirection=PhotDir_uy;
		zDirection=PhotDir_uz;
	}
	
	if (FlatEle) {
		fParticleGun->SetParticleDefinition(	G4ParticleTable::GetParticleTable()->FindParticle("e-"));
		G4double randomEne=G4UniformRand()*3;
		fParticleGun->SetParticleEnergy(randomEne*MeV); //SetParticleEnergy uses kinetic energy
		evtPrimAction->SetSourceEne(randomEne);
	} else if (FlatGamma) {
//		fParticleGun->SetParticleDefinition(	G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
//		fParticleGun->SetParticleCharge(0);
		G4double randomEne=G4UniformRand()*1;
		evtPrimAction->SetSourceEne(randomEne);
		Source_X=Sphere_Radius*Sphere_X;
		Source_Y=Sphere_Radius*Sphere_Y;
		if (G4PhysicalVolumeStore::GetInstance()->GetVolume("MiddleCase")!=NULL)	Sphere_ZOffset=G4PhysicalVolumeStore::GetInstance()->GetVolume("MiddleCase")->GetTranslation().z();
		Source_Z=Sphere_Radius*Sphere_Z+Sphere_ZOffset;
		G4ParticleDefinition* fotone = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
		fParticleGun->SetParticleDefinition(fotone);
		fParticleGun->SetParticleEnergy(randomEne);
		xDirection=PhotDir_ux;
		yDirection=PhotDir_uy;
		zDirection=PhotDir_uz;
	}
	
	
	const	G4ThreeVector position = G4ThreeVector(Source_X, Source_Y, Source_Z);

	
	//###################################################
	// Sampling particle initial direction
	//##########################
	if (FlatEle) { //If FlatSource (for Eff) was requested, generate only towards up
		G4double phi = G4UniformRand()*CLHEP::pi*2.;
		G4double costheta = G4UniformRand();
		G4double theta = acos(costheta);
		 xDirection = sin(theta)*cos(phi);
		 yDirection = sin(theta)*sin(phi);
		 zDirection = costheta;

//		const G4ThreeVector momentumDirection = G4ThreeVector(xDirection,yDirection,zDirection);
//		fParticleGun->SetParticleMomentumDirection(momentumDirection);
	} else {
//		G4ThreeVector momentumDirection = G4ThreeVector(0,0,0);
//		fParticleGun->SetParticleMomentumDirection(momentumDirection);
	}
	
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xDirection,yDirection,zDirection));
	fParticleGun->SetParticlePosition(position);
	
	
//	evtPrimAction->SetSourceEne(fParticleGun->GetParticleEnergy());
	evtPrimAction->SetSourceX((position.x())/mm);
	evtPrimAction->SetSourceY((position.y())/mm);
	evtPrimAction->SetSourceZ((position.z())/mm);
	
	
	fParticleGun->GeneratePrimaryVertex(anEvent);
	
	if(anEvent->GetEventID()==1) {  //stampo informazioni sorgente
		G4cout<<"Dimensioni sorgente: Raggio interno = "<<fRadiusInt<<", Raggio esterno = "<<fRadiusExt<<", H = "<<fZ<<G4endl;
		if (fSourceSelect==3) { //solo se è la sorgente ExtY..
			G4cout<<"TBR richiesto= "<<fTBR<<G4endl;
			G4cout<<"VolA= "<<VolA<<", ProbA= "<<ProbA<<G4endl;
			G4cout<<"VolB= "<<VolB<<", ProbB= "<<ProbB<<G4endl;
			G4cout<<"VolC= "<<VolC<<", ProbC= "<<ProbC<<G4endl;
			G4cout<<"Volume sorgente tot= "<<VolA+VolB+VolC<<G4endl;
		}
	}
}



G4double  B1PrimaryGeneratorAction::BetaDecaySpectrum(G4double Ek, G4double EndPoint)
{
	G4double ElMassMev = 0.510998928*MeV;
	
	G4double res=0.;
	
	G4double omega= Ek /ElMassMev  +1.;
	G4double omegamax= EndPoint /ElMassMev +1.;
	if(omega>1 && omega<omegamax)
	{
		res=(omegamax-omega)*(omegamax-omega) *Ek*sqrt(omega*omega-1.);
	}
	return res;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

