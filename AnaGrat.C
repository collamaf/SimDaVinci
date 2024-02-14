{
	// ############ Macro per analizzare il grattugiamento del PTER con Ga68 - 23.05.2022 -> 12.2.2024

	const int NThr = 6;
	const int NFiles = 5;
	double eThr[NThr] = {0, 30, 55, 67, 150, 300};
	// double spessori[NFiles]={3,2.5,2,1.5,1};
	// double spessori[NFiles] = {3, 2, 1.5, 1, 0.8, 0.6, 0.2};
	double spessori[NFiles] = {3, 2, 1, 0.8, 0.5};
	const int eMaxSign = 800;
	const int eMaxBkg = 400;
	const int nBin = 200;
	const int eMaxZoom = 200;
	double NSiPmBSign[NFiles];
	double NSiPmGSign[NFiles];
	double NSiPmBBkg[NFiles];
	double NSiPmGBkg[NFiles];

	TGraph *graphSign[NThr];
	TGraph *graphBkg[NThr];
	TGraph *graphRatio[NThr];

	TGraph *graphSiPmBSign = new TGraph();
	graphSiPmBSign->SetNameTitle("graphSiPmBSign", "Direct Betas to SiPM for Signal; Spessore Pter [mm]; ");
	graphSiPmBSign->SetMarkerStyle(20);
	TGraph *graphSiPmGSign = new TGraph();
	graphSiPmGSign->SetNameTitle("graphSiPmGSign", "Direct Gammas to SiPM for Signal; Spessore Pter [mm]; ");
	graphSiPmGSign->SetMarkerStyle(20);

	TGraph *graphSiPmBBkg = new TGraph();
	graphSiPmBBkg->SetNameTitle("graphSiPmBBkg", "Direct Betas to SiPM for Bkg; Spessore Pter [mm]; ");
	graphSiPmBBkg->SetMarkerStyle(20);
	TGraph *graphSiPmGBkg = new TGraph();
	graphSiPmGBkg->SetNameTitle("graphSiPmGBkg", "Direct Gammas to SiPM for Bkg; Spessore Pter [mm]; ");
	graphSiPmGBkg->SetMarkerStyle(20);

	TFile *fFileSign[NFiles];
	TFile *fFileBkg[NFiles];

	TMultiGraph *multiGraphSign = new TMultiGraph();
	TMultiGraph *multiGraphBkg = new TMultiGraph();
	TMultiGraph *multiGraphAll = new TMultiGraph();
	TMultiGraph *multiGraphRatio = new TMultiGraph();

	TMultiGraph *multiGraphSiPMSign = new TMultiGraph();
	TMultiGraph *multiGraphSiPMBkg = new TMultiGraph();

	double dataSign[NFiles][NThr];
	double dataBkg[NFiles][NThr];

	double dataSignNorm[NFiles][NThr];
	double dataBkgNorm[NFiles][NThr];

	TFile *fOut = new TFile(Form("%s.root", "OutGratF18"), "RECREATE");
	THStack *hStackSign = new THStack("hstackSign", "hstackSign;Eabs [keV];");
	THStack *hStackBkg = new THStack("hstackBkg", "hstackBkg;Eabs [keV];");

	THStack *hStackSignZ = new THStack("hstackSignZ", "hstackSignZ;Eabs [keV];");
	THStack *hStackBkgZ = new THStack("hstackBkgZ", "hstackBkgZ;Eabs [keV];");

	THStack *hStackSignCum = new THStack("hstackSignCum", "hstackSignCum;Eabs [keV];");
	THStack *hStackBkgCum = new THStack("hstackBkgCum", "hstackBkgCum;Eabs [keV];");

	THStack *hStackSiPmBSign = new THStack("hStackSiPmBSign", "hStackSiPmBSign;E [keV];");
	THStack *hStackSiPmBBkg = new THStack("hStackSiPmBBkg", "hStackSiPmBBkg;E [keV];");

	THStack *hStackSiPmGSign = new THStack("hStackSiPmGSign", "hStackSiPmGSign;E [keV];");
	THStack *hStackSiPmGBkg = new THStack("hStackSiPmGBkg", "hStackSiPmGBkg;E [keV];");

	for (int iThr = 0; iThr < NThr; iThr++)
	{
		graphSign[iThr] = new TGraph();
		graphSign[iThr]->SetNameTitle(Form("Sign%d", (int)eThr[iThr]), Form("grSign%d", (int)eThr[iThr]));
		graphBkg[iThr] = new TGraph();
		graphBkg[iThr]->SetNameTitle(Form("Bkg%d", (int)eThr[iThr]), Form("grBkg%d", (int)eThr[iThr]));
		graphRatio[iThr] = new TGraph();
		graphRatio[iThr]->SetNameTitle(Form("Ratio%d", (int)eThr[iThr]), Form("grRatio%d", (int)eThr[iThr]));

		graphSign[iThr]->SetMarkerStyle(20);
		graphBkg[iThr]->SetMarkerStyle(22);
		graphRatio[iThr]->SetMarkerStyle(21);

		graphSign[iThr]->SetMarkerColor(iThr + 1);
		graphBkg[iThr]->SetMarkerColor(iThr + 1);
		graphRatio[iThr]->SetMarkerColor(iThr + 1);

		graphSign[iThr]->SetLineColor(iThr + 1);
		graphBkg[iThr]->SetLineColor(iThr + 1);
		graphRatio[iThr]->SetLineColor(iThr + 1);
	}

	// //File del 2022 - Grattugiamento PTER ma su Ga esteso
	// 	fFileSign[0] = TFile::Open("PTERmc_PDiam60_PDz30_X0_Z0_NoAbs_Z31_A68_Diam60_Dz35_Gratt_N1000000.root");
	// 	fFileSign[1] = TFile::Open("PTERmc_PDiam60_PDz25_X0_Z0_NoAbs_Z31_A68_Diam60_Dz35_Gratt_N1000000.root");
	// 	fFileSign[2] = TFile::Open("PTERmc_PDiam60_PDz20_X0_Z0_NoAbs_Z31_A68_Diam60_Dz35_Gratt_N1000000.root");
	// 	fFileSign[3] = TFile::Open("PTERmc_PDiam60_PDz15_X0_Z0_NoAbs_Z31_A68_Diam60_Dz35_Gratt_N1000000.root");
	// 	fFileSign[4] = TFile::Open("PTERmc_PDiam60_PDz10_X0_Z0_NoAbs_Z31_A68_Diam60_Dz35_Gratt_N1000000.root");

	// 	fFileBkg[0] = TFile::Open("PTERmc_PDiam60_PDz30_X0_Z0_NoAbs_511keV_Gratt_Light_N50000000.root");
	// 	fFileBkg[1] = TFile::Open("PTERmc_PDiam60_PDz25_X0_Z0_NoAbs_511keV_Gratt_Light_N50000000.root");
	// 	fFileBkg[2] = TFile::Open("PTERmc_PDiam60_PDz20_X0_Z0_NoAbs_511keV_Gratt_Light_N50000000.root");
	// 	fFileBkg[3] = TFile::Open("PTERmc_PDiam60_PDz15_X0_Z0_NoAbs_511keV_Gratt_Light_N50000000.root");
	// 	fFileBkg[4] = TFile::Open("PTERmc_PDiam60_PDz10_X0_Z0_NoAbs_511keV_Gratt_Light_N50000000.root");

	// File del 2024 - Grattugiamento PTER ma su Sr e Ba in lab con centratore (quindi ~300um ABS davanti)
	// fFileSign[0] = TFile::Open("PTERmc_PDiam120_PDz30_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	// fFileSign[1] = TFile::Open("PTERmc_PDiam120_PDz15_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	// fFileSign[2] = TFile::Open("PTERmc_PDiam120_PDz10_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	// fFileSign[3] = TFile::Open("PTERmc_PDiam120_PDz8_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	// fFileSign[4] = TFile::Open("PTERmc_PDiam120_PDz6_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	// fFileSign[5] = TFile::Open("PTERmc_PDiam120_PDz2_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");

	// fFileBkg[0] = TFile::Open("PTERmc_PDiam120_PDz30_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[1] = TFile::Open("PTERmc_PDiam120_PDz15_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[2] = TFile::Open("PTERmc_PDiam120_PDz10_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[3] = TFile::Open("PTERmc_PDiam120_PDz8_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[4] = TFile::Open("PTERmc_PDiam120_PDz6_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[5] = TFile::Open("PTERmc_PDiam120_PDz2_PVCDiam160_PVCLT19_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");

	// FIles per confronto con Lab: Sr-Ba
	//  fFileSign[0] = TFile::Open("PTERmc_PDiam60_PDz30_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	//  fFileSign[1] = TFile::Open("PTERmc_PDiam60_PDz20_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	//  fFileSign[2] = TFile::Open("PTERmc_PDiam60_PDz15_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	//  fFileSign[3] = TFile::Open("PTERmc_PDiam60_PDz10_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	//  fFileSign[4] = TFile::Open("PTERmc_PDiam60_PDz8_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	//  fFileSign[5] = TFile::Open("PTERmc_PDiam60_PDz6_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");
	//  fFileSign[6] = TFile::Open("PTERmc_PDiam60_PDz2_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_ExtSr_RESTART_N1000000.root");

	// fFileBkg[0] = TFile::Open("PTERmc_PDiam60_PDz30_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[1] = TFile::Open("PTERmc_PDiam60_PDz20_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[2] = TFile::Open("PTERmc_PDiam60_PDz15_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[3] = TFile::Open("PTERmc_PDiam60_PDz10_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[4] = TFile::Open("PTERmc_PDiam60_PDz8_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[5] = TFile::Open("PTERmc_PDiam60_PDz6_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");
	// fFileBkg[6] = TFile::Open("PTERmc_PDiam60_PDz2_PVCDiam120_PVCLT29_X0_Z3_AbsDz300_AbsHoleD0_AbsMatABS_RESTART_N1000000.root");

	// FIles per previsione test F18
	fFileSign[0] = TFile::Open("PTERmc_PDiam60_PDz30_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ExtGa_GaSet3_AluCaseT155_AppMat1_NEW_N1000000.root");
	fFileSign[1] = TFile::Open("PTERmc_PDiam60_PDz20_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ExtGa_GaSet3_AluCaseT155_AppMat1_NEW_N1000000.root");
	fFileSign[2] = TFile::Open("PTERmc_PDiam60_PDz10_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ExtGa_GaSet3_AluCaseT155_AppMat1_NEW_N1000000.root");
	fFileSign[3] = TFile::Open("PTERmc_PDiam60_PDz8_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ExtGa_GaSet3_AluCaseT155_AppMat1_NEW_N1000000.root");
	fFileSign[4] = TFile::Open("PTERmc_PDiam60_PDz5_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ExtGa_GaSet3_AluCaseT155_AppMat1_NEW_N1000000.root");

	fFileBkg[0] = TFile::Open("PTERmc_PDiam60_PDz30_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ContF18_Thick9_Dz80_CATAFLUORO_Light_N100000000.root");
	fFileBkg[1] = TFile::Open("PTERmc_PDiam60_PDz20_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ContF18_Thick9_Dz80_CATAFLUORO_Light_N100000000.root");
	fFileBkg[2] = TFile::Open("PTERmc_PDiam60_PDz10_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ContF18_Thick9_Dz80_CATAFLUORO_Light_N100000000.root");
	fFileBkg[3] = TFile::Open("PTERmc_PDiam60_PDz8_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ContF18_Thick9_Dz80_CATAFLUORO_Light_N100000000.root");
	fFileBkg[4] = TFile::Open("PTERmc_PDiam60_PDz5_PVCDiam120_PVCLT20_X0_Z0_NoAbs_ContF18_Thick9_Dz80_CATAFLUORO_Light_N100000000.root");

	//	dataBkg[0][0]=1;
	//	dataSign[0][0]=1;
	// TCanvas *canvEneSign = new TCanvas("canvEneSign", "canvEneSign");
	// TCanvas *canvEneBkg = new TCanvas("canvEneBkg", "canvEneBkg");
	// TCanvas *canvEneSignZoom = new TCanvas("canvEneSignZoom", "canvEneSignZoom");
	// TCanvas *canvEneBkgZoom = new TCanvas("canvEneBkgZoom", "canvEneBkgZoom");
	//	TCanvas* canvTemp=new TCanvas();

	// TH1F* SpettroFondoLontano=new TH1F("SpettroFondoLontano","SpettroFondoLontano; E dep [keV];", nBin, 0, eMax);
	// if (caso==0 || caso==2)	B1->Draw("Eabs>>SpettroFondoLontano","Eabs>00");

	// TH1* histoTemp=new TH1F("histoTemp","histoTemp", 200, 0, 1200);
	// cout<<"CIAO PRE "<<hStack->GetNhists()<<endl;
	for (int iFile = 0; iFile < NFiles; iFile++)
	{
		cout << "File= " << iFile << endl;
		fFileSign[iFile]->cd();
		// canvEneSign->cd();
		// iFile==0? B1->Draw("Eabs>>histoTemp","Eabs>0",""):B1->Draw("Eabs>>histoTemp","Eabs>0","samesPLC");
		TH1F *histoTempSign = new TH1F("histoTempSign", Form("Sign_%:.1f;Eabs [keV];", spessori[iFile]), nBin, 0, eMaxSign);
		TH1F *histoTempSignZ = new TH1F("histoTempSignZ", Form("Sign_%:.1f;Eabs [keV];", spessori[iFile]), nBin, 0, eMaxZoom);

		TH1F *histoTempSiPmBSign = new TH1F("histoTempSiPmBSign", Form("SignB_%:.1f;Eabs [keV];", spessori[iFile]), nBin / 4, 0, eMaxSign);
		TH1F *histoTempSiPmGSign = new TH1F("histoTempSiPmGSign", Form("SignG_%:.1f;Eabs [keV];", spessori[iFile]), nBin / 4, 0, eMaxSign);

		B1->Draw("Eabs>>histoTempSign", "Eabs>0", "goff");
		B1->Draw("Eabs>>histoTempSignZ", "Eabs>0&&Eabs<200", "goff");
		NSiPmBSign[iFile] = B1->Draw("SiPMEne>>histoTempSiPmBSign", "fabs(SiPMPart)==11", "goff");
		NSiPmGSign[iFile] = B1->Draw("SiPMEne>>histoTempSiPmGSign", "SiPMPart==22", "goff");
		// histoTemp = (TH1F*)gPad->GetPrimitive("htemp");
		fOut->cd();
		histoTempSign->Write();
		histoTempSignZ->Write();
		histoTempSiPmBSign->Write();
		histoTempSiPmGSign->Write();
		hStackSign->Add(histoTempSign);
		hStackSignZ->Add(histoTempSignZ);
		hStackSiPmBSign->Add(histoTempSiPmBSign);
		hStackSiPmGSign->Add(histoTempSiPmGSign);
		TH1 *histoTempSignCum = histoTempSign->GetCumulative();
		hStackSignCum->Add(histoTempSignCum);

		fFileBkg[iFile]->cd();
		// canvEneBkg->cd();
		TH1F *histoTempBkg = new TH1F("histoTempBkg", Form("Bkg_%:.1f;Eabs [keV];", spessori[iFile]), nBin, 0, eMaxBkg);
		TH1F *histoTempBkgZ = new TH1F("histoTempBkgZ", Form("Bkg_%:.1f;Eabs [keV];", spessori[iFile]), nBin, 0, eMaxZoom);

		TH1F *histoTempSiPmBBkg = new TH1F("histoTempSiPmBBkg", Form("BkgB_%:.1f;Eabs [keV];", spessori[iFile]), nBin / 4, 0, eMaxSign);
		TH1F *histoTempSiPmGBkg = new TH1F("histoTempSiPmGBkg", Form("BkgG_%:.1f;Eabs [keV];", spessori[iFile]), nBin / 4, 0, eMaxSign);

		B1->Draw("Eabs>>histoTempBkg", "Eabs>0", "goff");
		B1->Draw("Eabs>>histoTempBkgZ", "Eabs>0&&Eabs<200", "goff");
		NSiPmBBkg[iFile] = B1->Draw("SiPMEne>>histoTempSiPmBBkg", "fabs(SiPMPart)==11", "goff");
		NSiPmGBkg[iFile] = B1->Draw("SiPMEne>>histoTempSiPmGBkg", "SiPMPart==22", "goff");
		// histoTemp = (TH1F*)gPad->GetPrimitive("htemp");
		fOut->cd();
		histoTempBkg->Write();
		histoTempBkgZ->Write();
		histoTempSiPmBBkg->Write();
		histoTempSiPmGBkg->Write();
		hStackBkg->Add(histoTempBkg);
		hStackBkgZ->Add(histoTempBkgZ);
		hStackSiPmBBkg->Add(histoTempSiPmBBkg);
		hStackSiPmGBkg->Add(histoTempSiPmGBkg);

		graphSiPmBSign->AddPoint(spessori[iFile], NSiPmBSign[iFile]);
		graphSiPmGSign->AddPoint(spessori[iFile], NSiPmGSign[iFile]);

		graphSiPmBBkg->AddPoint(spessori[iFile], NSiPmBBkg[iFile]);
		graphSiPmGBkg->AddPoint(spessori[iFile], NSiPmGBkg[iFile]);

		// fFileBkg[iFile]->cd();
		// canvEneBkg->cd();
		// iFile==0? B1->Draw("Eabs>>histo(180)","Eabs>0",""):B1->Draw("Eabs>>histo(180)","Eabs>0","samesPLC");

		// fFileSign[iFile]->cd();
		// canvEneSignZoom->cd();
		// iFile==0? B1->Draw("Eabs>>histoZ(180)","Eabs>0&&Eabs<200",""):B1->Draw("Eabs>>histoZ(180)","Eabs>0&&Eabs<200","samesPLC"); //Ci concentriamo a bassa energia depositata, dove gioca la soglia
		// fFileBkg[iFile]->cd();
		// canvEneBkgZoom->cd();
		// iFile==0? B1->Draw("Eabs>>histoZ(180)","Eabs>0&&Eabs<200",""):B1->Draw("Eabs>>histoZ(180)","Eabs>0&&Eabs<200","samesPLC");

		for (int iThr = 0; iThr < NThr; iThr++)
		{
			char *whatToPlot = Form("Eabs>%f", eThr[iThr]);
			fFileSign[iFile]->cd();

			//			canvTemp->cd();
			dataSign[iFile][iThr] = B1->Draw("Eabs", whatToPlot, "goff");
			fFileBkg[iFile]->cd();
			dataBkg[iFile][iThr] = B1->Draw("Eabs", whatToPlot, "goff");
			dataSignNorm[iFile][iThr] = dataSign[iFile][iThr] / dataSign[0][iThr];
			dataBkgNorm[iFile][iThr] = dataBkg[iFile][iThr] / dataBkg[0][iThr];

			cout << "Thr=\t" << eThr[iThr] << "\t" << dataSign[iFile][iThr] << "\t" << dataSignNorm[iFile][iThr] << "\t" << dataBkg[iFile][iThr] << "\t" << dataBkgNorm[iFile][iThr] << endl;
			graphSign[iThr]->AddPoint(spessori[iFile], dataSignNorm[iFile][iThr]);
			graphBkg[iThr]->AddPoint(spessori[iFile], dataBkgNorm[iFile][iThr]);
			graphRatio[iThr]->AddPoint(spessori[iFile], dataSignNorm[iFile][iThr] / dataBkgNorm[iFile][iThr]);
		}
	}

	// cout<<"CIAO POST "<<hStack->GetNhists()<<endl;

	for (int iThr = 0; iThr < NThr; iThr++)
	{
		multiGraphSign->Add(graphSign[iThr]);
		multiGraphBkg->Add(graphBkg[iThr]);

		multiGraphAll->Add(graphSign[iThr]);
		multiGraphAll->Add(graphBkg[iThr]);

		multiGraphRatio->Add(graphRatio[iThr]);
	}

	TCanvas *canvBkg = new TCanvas("canvBkg", "canvBkg");
	multiGraphBkg->Draw("APL");
	multiGraphBkg->SetTitle("Fondo; Spessore [mm]; Perdita Segnale Bkg [%]");
	canvBkg->BuildLegend();
	canvBkg->SaveAs("Gratt_ScanBkg.pdf");

	TCanvas *canvSign = new TCanvas("canvSign", "canvSign");
	multiGraphSign->Draw("APL");
	multiGraphSign->SetTitle("Segnale; Spessore [mm]; Perdita Segnale Sign [%]");
	canvSign->BuildLegend();
	canvSign->SaveAs("Gratt_ScanSign.pdf");

	TCanvas *cAll = new TCanvas("canvAll", "canvAll");
	multiGraphAll->Draw("APL");
	multiGraphAll->SetTitle("; Spessore [mm]; Perdita Segnale [%]");
	cAll->BuildLegend();
	cAll->SaveAs("Gratt_ScanAll.pdf");

	TCanvas *cRatio = new TCanvas("canvRatio", "canvRatio");
	multiGraphRatio->Draw("APL");
	multiGraphRatio->SetTitle("Rapporto; Spessore [mm]; Rapporto");
	cRatio->BuildLegend();
	cRatio->SaveAs("Gratt_Ratio.pdf");

	fOut->cd();
	hStackSign->Write();
	hStackBkg->Write();
	hStackSignZ->Write();
	hStackBkgZ->Write();
	hStackSiPmBBkg->Write();
	hStackSiPmGBkg->Write();
	hStackSiPmBSign->Write();
	hStackSiPmGSign->Write();
	graphSiPmBSign->Write();
	graphSiPmGSign->Write();
	graphSiPmBBkg->Write();
	graphSiPmGBkg->Write();

	TCanvas *cStackSign = new TCanvas("cStackSign", "cStackSign");
	hStackSign->Draw("nostackPLC");
	cStackSign->BuildLegend();
	cStackSign->SaveAs("Gratt_Sign.pdf");

	TCanvas *cStackBkg = new TCanvas("cStackBkg", "cStackBkg");
	hStackBkg->Draw("nostackPLC");
	cStackBkg->SetLogy();
	cStackBkg->BuildLegend();
	cStackBkg->SaveAs("Gratt_Bkg.pdf");

	TCanvas *cStackSignZ = new TCanvas("cStackSignZ", "cStackSignZ");
	hStackSignZ->Draw("nostackPLC");
	cStackSignZ->BuildLegend();
	cStackSignZ->SaveAs("Gratt_SignZ.pdf");

	TCanvas *cStackBkgZ = new TCanvas("cStackBkgZ", "cStackBkgZ");
	hStackBkgZ->Draw("nostackPLC");
	cStackBkgZ->SetLogy();
	cStackBkgZ->BuildLegend();
	cStackBkgZ->SaveAs("Gratt_BkgZ.pdf");

	TCanvas *cStackSignCum = new TCanvas("cStackSignCum", "cStackSignCum");
	hStackSignCum->Draw("nostackPLC");
	cStackSignCum->BuildLegend();
	cStackSignCum->SaveAs("Gratt_Cum.pdf");

	TCanvas *canvSiPmBSign = new TCanvas("canvSiPmBSign", "canvSiPmBSign");
	hStackSiPmBSign->Draw("nostackPLC");
	// multiGraphBkg->SetTitle("Fondo; Spessore [mm]; Perdita Segnale Bkg [%]");
	canvSiPmBSign->BuildLegend();
	canvSiPmBSign->SaveAs("Gratt_SiPmBSign.pdf");

	TCanvas *canvSiPmGSign = new TCanvas("canvSiPmGSign", "canvSiPmGSign");
	hStackSiPmGSign->Draw("nostackPLC");
	canvSiPmGSign->BuildLegend();
	canvSiPmGSign->SetLogy();
	canvSiPmGSign->SaveAs("Gratt_SiPmGSign.pdf");

	TCanvas *canvSiPmBBkg = new TCanvas("canvSiPmBBkg", "canvSiPmBBkg");
	hStackSiPmBBkg->Draw("nostackPLC");
	canvSiPmBBkg->BuildLegend();
	canvSiPmBBkg->SaveAs("Gratt_SiPmBBkg.pdf");

	TCanvas *canvSiPmGBkg = new TCanvas("canvSiPmGBkg", "canvSiPmGBkg");
	hStackSiPmGBkg->Draw("nostackPLC");
	canvSiPmGBkg->BuildLegend();
	canvSiPmGBkg->SetLogy();
	canvSiPmGBkg->SaveAs("Gratt_SiPmGBkg.pdf");

	TCanvas *canvSiPmGrattBBkg = new TCanvas("canvSiPmGrattBBkg", "canvSiPmGrattBBkg");
	graphSiPmBBkg->Draw("APL");
	canvSiPmGrattBBkg->SaveAs("Gratt_SiPmGrattBBkg.pdf");

	TCanvas *canvSiPmGrattGBkg = new TCanvas("canvSiPmGrattGBkg", "canvSiPmGrattGBkg");
	graphSiPmGBkg->Draw("APL");
	canvSiPmGrattGBkg->SaveAs("Gratt_SiPmGrattGBkg.pdf");

	TCanvas *canvSiPmGrattBSign = new TCanvas("canvSiPmGrattBSign", "canvSiPmGrattBSign");
	graphSiPmBSign->Draw("APL");
	canvSiPmGrattBSign->SaveAs("Gratt_SiPmGrattBSign.pdf");

	TCanvas *canvSiPmGrattGSign = new TCanvas("canvSiPmGrattGSign", "canvSiPmGrattGSign");
	graphSiPmGSign->Draw("APL");
	canvSiPmGrattGSign->SaveAs("Gratt_SiPmGrattGSign.pdf");
}
