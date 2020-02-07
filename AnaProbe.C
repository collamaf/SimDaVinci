#define AnaProbe_cxx
#include "AnaProbe.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#define NBINE 200
#define NBINERED 50
#define EMIN 0
#define EMAX 2300
#define NBINZ 50
#define ZMIN -10
#define ZMAX 0
#define NANAPART 3
#define NTHR 4

//
//  Root Macro to Analyse probe simulation output
//  Last modification: 2020.02.07 by collamaf
//
//
//
//

void AnaProbe::Loop()
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	TString partName[NANAPART] = {"Ele","Pos","Fot"};
	
	std::vector<double>ParticlesId={
		11, //e-
		-11, //e+
		22, //gamma
	};
	int ParticlePosition;
	
	double eThr[NTHR]= {0, 67, 150, 600};
	
	TH1F* hExitSource[NANAPART];
	TH1F* hPostAbs[NANAPART];
	TH1F* hPreProbe[NANAPART];
	TH1F* hPrePter[NANAPART];
	TH1F* hEnDep[NANAPART];

	TH1F*	hSourceZ[NTHR];
	
	for (int ii=0; ii<NTHR; ii++) {
		hSourceZ[ii]= new TH1F(Form("hSourceZThr%d",(int)eThr[ii]),Form("Source Depth for Edep Thr = %d keV; Source Depth [mm]; ",(int)eThr[ii]), NBINZ, ZMIN, ZMAX);
	}
	
	Color_t colori[4]={kBlack, kRed, kGreen, kBlue};

	TCanvas * cSumUp[NANAPART];
	TCanvas * cSourceZ=new TCanvas("cSourceZ","cSourceZ");
	TCanvas * cSourceZEabs=new TCanvas("cSourceZEabs","cSourceZEabs");

	TH2F * hSourceZEabs = new TH2F("hSourceZEabs","Source Depth vs Edep; Edep [keV]; Source Depth [mm]", NBINERED, EMIN, EMAX, NBINZ, ZMIN, ZMAX);

	for (int ii=0; ii<NANAPART; ii++) {
		hExitSource[ii]= new TH1F(Form("hExitSource%s", partName[ii].Data()),Form("Exiting Source %s; E [keV];", partName[ii].Data()), NBINE, EMIN, EMAX);
		hPostAbs[ii]= new TH1F(Form("hPostAbs%s", partName[ii].Data()),Form("Post Absorber %s; E [keV];", partName[ii].Data()), NBINE, EMIN, EMAX);
		hPreProbe[ii]= new TH1F(Form("hPreProbe%s", partName[ii].Data()),Form("Pre Probe %s; E [keV];", partName[ii].Data()), NBINE, EMIN, EMAX);
		hPrePter[ii]= new TH1F(Form("hPrePter%s", partName[ii].Data()),Form("Pre Pter %s; E [keV];", partName[ii].Data()), NBINE, EMIN, EMAX);
		hEnDep[ii]= new TH1F(Form("hEnDep%s", partName[ii].Data()),Form("Depositing Energy %s; E [keV];", partName[ii].Data()), NBINE, EMIN, EMAX);
		cSumUp[ii] = new TCanvas(Form("cSumUp%s",partName[ii].Data()),Form("cSumUp%s",partName[ii].Data()));
		
		hExitSource[ii]->SetLineColor(colori[0]);
		hPostAbs[ii]->SetLineColor(colori[1]);
		hPreProbe[ii]->SetLineColor(colori[2]);
		hPrePter[ii]->SetLineColor(colori[3]);
	}

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		if (jentry%(nentries/10)==0) cout<<"I am analyzing entry: "<<jentry<<" of "<<nentries<<" - "<< (int)(((double)jentry)/nentries*100 )<<" %"<<endl;
		
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		
		
		for (int ii=0; ii<PrePterTrackN; ii++) {
			ParticlePosition= std::find(ParticlesId.begin(), ParticlesId.end(), PrePterPart->at(ii)) - ParticlesId.begin();
			hPrePter[ParticlePosition]->Fill(PrePterEn->at(ii));
		}
		for (int ii=0; ii<PreProbeTrackN; ii++) {
			ParticlePosition= std::find(ParticlesId.begin(), ParticlesId.end(), PreProbePart->at(ii)) - ParticlesId.begin();
			hPreProbe[ParticlePosition]->Fill(PreProbeEn->at(ii));
		}
		for (int ii=0; ii<ExitEne->size(); ii++) {
			ParticlePosition= std::find(ParticlesId.begin(), ParticlesId.end(), ExitPart->at(ii)) - ParticlesId.begin();
			hExitSource[ParticlePosition]->Fill(ExitEne->at(ii));
		}
		for (int ii=0; ii<SourceEne->size(); ii++) {
			ParticlePosition= std::find(ParticlesId.begin(), ParticlesId.end(), SourcePart->at(ii)) - ParticlesId.begin();
			if (Eabs>0) hEnDep[ParticlePosition]->Fill(SourceEne->at(ii));
		}
		for (int ii=0; ii<PostAbsTrackN; ii++) {
			ParticlePosition= std::find(ParticlesId.begin(), ParticlesId.end(), PostAbsPart->at(ii)) - ParticlesId.begin();
			hPostAbs[ParticlePosition]->Fill(PostAbsEn->at(ii));
		}
		
		for (int ii = 0; ii <NTHR; ii++) {
			if (Eabs>eThr[ii])  {
				hSourceZ[ii]->Fill(SourceZ);
			}
		}

		hSourceZEabs->Fill(Eabs, SourceZ);
		
	} //fine loop sulle entries
	
	
	for (int ii=0; ii<NANAPART; ii++) {
		cSumUp[ii]->cd();
		hExitSource[ii]->Draw();
		hPostAbs[ii]->Draw("sames");
		hPreProbe[ii]->Draw("sames");
		hPrePter[ii]->Draw("sames");
		cSumUp[ii]->BuildLegend();
		cSumUp[ii]->Write();

		
		hExitSource[ii]->Write();
		hPostAbs[ii]->Write();
		hPreProbe[ii]->Write();
		hPrePter[ii]->Write();
		hEnDep[ii]->Write();

	}
	
	cSourceZ->cd();

	for (int ii=0; ii<NTHR; ii++) {
		hSourceZ[ii]->Draw("samesPLC");
		hSourceZ[ii]->Write();
	}
	cSourceZ->BuildLegend();
	
	cSourceZEabs->cd();
	hSourceZEabs->Draw("colz");
	cSourceZEabs->SetLogz();
	cSourceZEabs->Write();

	hSourceZEabs->Write();
	
}
