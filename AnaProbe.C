#define AnaProbe_cxx
#include "AnaProbe.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#define NBINE 200
#define NBINERED 50
#define EMIN 0
#define EMAX 2300

#define NBINZ 100
#define ZMIN -10
#define ZMAX 0

#define NBINX 100
#define XMIN -10
#define XMAX 10



//
//  Root Macro to Analyse probe simulation output
//  Last modification: 2020.02.19 by collamaf
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
	
	double eThr[NTHR]= {0, 67, 150, 300};
	double	zCumThr = 0.9;
	double probeDiam=6, probeDepth=3; //probe dimensions in mm
	probeDiam=strtod(inputFileName(inputFileName.Index("_PDiam")+6,1).Data(), NULL);
	probeDepth=strtod(inputFileName(inputFileName.Index("_PDz")+4,1).Data(), NULL);
	
	TH1F* hExitSource[NANAPART];
	TH1F* hPostAbs[NANAPART];
	TH1F* hPreProbe[NANAPART];
	TH1F* hPrePter[NANAPART];
	TH1F* hEnDep[NANAPART];

	TH1F*	hSourceZ[NTHR];
	TH1F*	hSourceZcum[NTHR];
	TH2F* hSourceZX[NTHR];
	
	for (int ii=0; ii<NTHR; ii++) {
		hSourceZ[ii]= new TH1F(Form("hSourceZThr%d",(int)eThr[ii]),Form("Source Depth for Edep Thr = %d keV; Source Depth [mm]; ",(int)eThr[ii]), NBINZ, ZMIN, ZMAX);
		hSourceZcum[ii]= new TH1F(Form("hSourceZcumThr%d",(int)eThr[ii]),Form("Cum. Source Depth for Edep Thr = %d keV; Source Depth [mm]; ",(int)eThr[ii]), NBINZ, ZMIN, ZMAX);
		hSourceZX[ii] = new TH2F(Form("hSourceZX%d",(int)eThr[ii]),Form("Position of primary particles giving a signal > %d keV ; X [mm]; Z [mm]",(int)eThr[ii]), NBINX, XMIN, XMAX, NBINZ, ZMIN, ZMAX);
		dirThr[ii]=outputfile->mkdir(Form("Thr%d",(int)eThr[ii]));
	}
	
	Color_t colori[5]={kBlack, kRed, kGreen, kBlue, kMagenta};

	TCanvas * cSourceZX[NTHR];
	TCanvas * cSumUp[NANAPART];
	TCanvas * cSourceZ=new TCanvas("cSourceZ","cSourceZ");
	TCanvas * cSourceZcum=new TCanvas("cSourceZcum","cSourceZcum");
	TCanvas * cSourceZEabs=new TCanvas("cSourceZEabs","cSourceZEabs");

	TH2F * hSourceZEabs = new TH2F("hSourceZEabs","Source Depth vs Edep; Edep [keV]; Source Depth [mm]", NBINERED, EMIN, EMAX, NBINZ, ZMIN, ZMAX);

	for (int ii=0; ii<NANAPART; ii++) {
		
		dirPart[ii]= outputfile->mkdir(partName[ii]);
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
		hEnDep[ii]->SetLineColor(colori[4]);
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
				hSourceZX[ii]->Fill(SourceX, SourceZ);
			}
		}

		if (Eabs>0) hSourceZEabs->Fill(Eabs, SourceZ);
		
	} //fine loop sulle entries
	
	
	for (int ii=0; ii<NANAPART; ii++) {
		cSumUp[ii]->cd();
		hExitSource[ii]->Draw();
		hPostAbs[ii]->Draw("sames");
		hPreProbe[ii]->Draw("sames");
		hPrePter[ii]->Draw("sames");
		hEnDep[ii]->Draw("sames");

		cSumUp[ii]->BuildLegend();
		dirCanvas->cd();
		cSumUp[ii]->Write();
		cSumUp[ii]->SaveAs(Form("%s/%s%s.pdf",canvDir.Data(),inputFileName.Data(),cSumUp[ii]->GetName()));
		dirPart[ii]->cd();
		
		hExitSource[ii]->Write();
		hPostAbs[ii]->Write();
		hPreProbe[ii]->Write();
		hPrePter[ii]->Write();
		hEnDep[ii]->Write();
		outputfile->cd();
	}
	

	for (int ii=0; ii<NTHR; ii++) {
		dirThr[ii]->cd();
		cSourceZ->cd();
		hSourceZ[ii]->Draw("samesPLC");
		hSourceZ[ii]->Write();
		cSourceZcum->cd();
		hSourceZcum[ii]=(TH1F*)hSourceZ[ii]->GetCumulative(kFALSE);
		hSourceZcum[ii]->SetName(Form("hSourceZcumThr%d",(int)eThr[ii]));
//		hSourceZcum[ii]->SetTitle("Source Depth of Edep signals cumulative function");
		hSourceZcum[ii]->Scale(1/hSourceZcum[ii]->GetMaximum());
		hSourceZcum[ii]->Draw("samePLC");
		double tempX= hSourceZcum[ii]->GetBinCenter(hSourceZcum[ii]->FindLastBinAbove(zCumThr));
		double tempY= hSourceZcum[ii]->GetBinContent(hSourceZcum[ii]->FindLastBinAbove(zCumThr));
		TLine* cumThrLine= new TLine(tempX, 0, tempX, tempY);
		cumThrLine->Draw("samePLC");
		hSourceZcum[ii]->Write();
		if (ii==1) cout<<" ZCUM Fraction: "<<100*zCumThr<<"% of signal comes from the last "<<-tempX<<"mm of tissue "<<endl;
	}
	dirCanvas->cd();
	cSourceZ->BuildLegend();
	cSourceZ->Write();
	cSourceZ->SaveAs(Form("%s/%s%s.pdf",canvDir.Data(),inputFileName.Data(),cSourceZ->GetName()));
	cSourceZcum->BuildLegend();
	cSourceZcum->Write();
	cSourceZcum->SaveAs(Form("%s/%s%s.pdf",canvDir.Data(),inputFileName.Data(),cSourceZcum->GetName()));
	outputfile->cd();

//	cSourceZX->Divide(2,2);
	for (int ii=0; ii<NTHR; ii++) {
		dirThr[ii]->cd();
		cSourceZX[ii]=new TCanvas(Form("cSourceZX%d", (int)eThr[ii]),Form("cSourceZX%d", (int)eThr[ii]));
//		cSourceZX->cd(ii);
		hSourceZX[ii]->Draw("colz");
		TLine* probeShadowX= new TLine(-probeDiam/2., -probeDepth, probeDiam/2., -probeDepth);
		TLine* probeShadowY1= new TLine(-probeDiam/2., -probeDepth, -probeDiam/2., 0);
		TLine* probeShadowY2= new TLine(probeDiam/2., -probeDepth, probeDiam/2., 0);
		probeShadowX->SetLineColor(kRed);
		probeShadowY1->SetLineColor(kRed);
		probeShadowY2->SetLineColor(kRed);
		probeShadowX->SetLineWidth(3);
		probeShadowY1->SetLineWidth(3);
		probeShadowY2->SetLineWidth(3);
		probeShadowX->SetLineStyle(5);
		probeShadowY1->SetLineStyle(5);
		probeShadowY2->SetLineStyle(5);
		probeShadowX->Draw("same");
		probeShadowY1->Draw("same");
		probeShadowY2->Draw("same");
		hSourceZX[ii]->Write();
		dirCanvas->cd();
		cSourceZX[ii]->Write();
		cSourceZX[ii]->SaveAs(Form("%s/%s%s.pdf",canvDir.Data(),inputFileName.Data(),cSourceZX[ii]->GetName()));

		outputfile->cd();
	}
	
	cSourceZEabs->cd();
	hSourceZEabs->Draw("colz");
	TLine* lineThr[NTHR];
	for (int ii=0; ii<NTHR; ii++) { //not sure I want to show here all the threshold lines...
		lineThr[ii]=new TLine(eThr[ii], ZMIN, eThr[ii], 0);
		lineThr[ii]->SetLineStyle(5);
		lineThr[ii]->SetLineColor(kRed);
		lineThr[ii]->SetLineWidth(3);
		if (ii==1) lineThr[ii]->Draw("same");
	}
//	cSourceZEabs->SetLogz();
	dirCanvas->cd();
	cSourceZEabs->Write();
	cSourceZEabs->SaveAs(Form("%s/%s%s.pdf",canvDir.Data(),inputFileName.Data(),cSourceZEabs->GetName()));
	outputfile->cd();
	hSourceZEabs->Write();
	
}
