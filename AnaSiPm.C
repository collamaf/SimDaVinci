#define AnaSiPm_cxx
#include "AnaSiPm.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void AnaSiPm::Loop()
{
	// USAGE
	// .L AnaSiPm.C
	// c=new AnaSiPm("build/CMOSmcX0_Z10_NOCuD_Fil0_TBR10_ExtSr_60035")
	// c->Loop()
	
	// This Macro generates an histogram with the number of fired pixels per primary with an energy above a given threshold; this can be choose by the user from line comand and is set from default to 5.0 KeV.
	
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	TH1F* NPixelALL=new TH1F("NPixelALL","Number of ALL fired pixel per primary particle", 10, 0, 10);
	TH1F* NPixelGOOD=new TH1F("NPixelGOOD",Form("Number of fired pixel per primary particle (Thr=%.1lf)", fSoglia), 10, 0, 10);
	TH1F* NPixelOLD=new TH1F("NPixelOLD","Number of fired pixel per primary particle OLD METHOD", 10, 0, 10);
	
	int NPixTot=137*137;
	int debug=0;
	//vector<int> PixelAcceso(NPixTot,0.);
	//int PixelAcceso[NPixTot]=0;
	vector<int> FiredPixels;
	vector<double> PixelEn;
	vector<int> PixelIndex;
	double Threshold=fSoglia;
	cout<<"SOGLIA RICHIESTA= "<<Threshold<<endl;
	int ii=0;
	vector<int> CounterAll(nentries,0.);
	vector<int> CounterGood(nentries,0.); //above thr
	vector<int> CounterOLD(nentries,0.);
	
	std::vector<int>::iterator iteratore, iteratoreOLD;                        // Declare an iterator to a vector of int - will be used to search across vectors
	
	Long64_t nbytes = 0, nb = 0;
	//	nentries/=10; //to go quickly
	//	nentries=10;
	
	
	for (Long64_t jentry=0; jentry<nentries;jentry++) {   // INIZIO CICLO SULLE PARTICELLE PRIMARIE
		
		if (jentry%(nentries/10)==0) cout<<"Sto analizzando l'evento "<<jentry<<" di "<<nentries<<": %= "<<100*jentry/nentries<<endl;
		
		Long64_t ientry = LoadTree(jentry);                 // Set the environment to read one entrye (see def in .h) (ientry is just a return value to check if the loading operation was successful)
		if (ientry < 0) break; //if no entry could be loaded (e.g. because the tree is ended) break the loop
		nb = fChain->GetEntry(jentry); //Here I load the entry
		nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		
		
		
		//NEW METHOD (With Threshold)
		
		
		for (ii=0; ii<PixelID->size(); ii++){               // PixelID->size() coincides with the total number of pixels that have a >0 energy release in tracks originated from the jth-primary (also its secondaries...)
			if (debug) cout<<" Evento n= "<<jentry<<" ii= "<<ii <<" Pixel acceso: "<< PixelID->at(ii)<<endl;
			if (InCmosEn->at(ii)>0) {
				iteratore = find(PixelIndex.begin(), PixelIndex.end(),PixelID->at(ii)); //Search in PixelIndex vector if the current ID of the pixel is already in
				if (iteratore == PixelIndex.end()) {           //Pixel's ID not found
					PixelEn.push_back(InCmosEn->at(ii));         //Add the energy of InCmosEn->at(ii) at PixelEn vector
					PixelIndex.push_back(PixelID->at(ii));       //Add the ID of the pixel hit in PixelIndex vector
					CounterAll[jentry]++;
					//	cout<<" NUOVO PIXEL "<< PixelID->at(ii)<<" : Evento n= "<<jentry<<" ii= "<<ii <<" Pixel acceso: "<< PixelID->at(ii)<<" ha E= "<<InCmosEn->at(ii)<<endl;
					
				} else{                                                              //If the pixel's ID is already in PixelIndex
					auto indice = std::distance(PixelIndex.begin(), iteratore);        //index of pixel's ID hit in ii-th hit
					PixelEn[indice]+=InCmosEn->at(ii);                                 //Upgrade the energy of the pixel with ID given by indice
					
					
					//	cout<<"PIXEL "<< PixelID->at(ii) <<" GIA PRESENTE Stava in posizione indice= "<<indice<<" ha E= "<<InCmosEn->at(ii)<< " ora ho EnPix= "<<PixelEn[indice]<<endl;
					//	 PixelEn.at(iteratore-PixelIndex.begin())+=InCmosEn->at(ii);
				}
			}
		}
		//cout<<"Per il "<<jentry<<" primario, l'energia nel pixel "<<index<<" vale: "<<EnPerPix[index]<<endl;
		
		
		
		
		// OLD METHOD (Without Threshold)
		
		
		for (ii=0; ii<PixelID->size(); ii++) {
			iteratoreOLD = find(FiredPixels.begin(), FiredPixels.end(),PixelID->at(ii));     // Search for the element  PixelID->at(ii) in the list of "already fired" pixels, to avoid double counting
			if (iteratoreOLD == FiredPixels.end()&& InCmosEn->at(ii)>0) {                    //if it is a "new" pixel...
				FiredPixels.push_back(PixelID->at(ii));         //...add the PixelID found as last member of the vector FiredPixels
				if (debug) cout<<" Nuovo, aggiungo!"<<endl;
			} else if (debug) cout<<" C'era gia, stica!"<<endl;
			
		}
		
		CounterOLD[jentry]=FiredPixels.size();         //Put the number of pixels that have been hit by tracks originated from the jth-primary (also its secondaries...) into the CounterOLD vector (has one entry per each primary particle)
		
		
		for (ii=0; ii<CounterAll[jentry]; ii++){
			if (PixelEn[ii]>Threshold) CounterGood[jentry]++;
		}
		if (debug) cout<<"Per il "<<jentry<<" primario, CounterAll "<<CounterAll[jentry]<<" CounterGood: "<<CounterGood[jentry]<<endl;
		//PixelAcceso[PixelID->at(ii)]=1;
		
		PixelEn.assign(PixelEn.size(),0);
		PixelEn.clear();
		
		PixelIndex.assign(PixelIndex.size(),0);
		PixelIndex.clear();
		
		//		CounterAll[jentry]=std::count(PixelAcceso.begin(), PixelAcceso.end(), 1);
		
		if (debug)		cout<<"                     "<<jentry<<" numero pixel accesi: "<<CounterAll[jentry] <<endl;
		NPixelALL->Fill(CounterAll[jentry]);     //Fills the NPixel histogram with the elements of CounterAll[jentry] (the number of different pixels hit in each primary event)
		NPixelGOOD->Fill(CounterGood[jentry]);   //Fills the NPixel histogram with the elements of CounterGood[jentry] (the number of different pixels hit in each primary event with energy above the threshold)
		NPixelOLD->Fill(CounterOLD[jentry]);
		
		FiredPixels.clear();               // Clear the vector FiredPixels and put his size to 0; ready for the next iteration (new primary particle)
		
		//	PixelAcceso.assign(PixelAcceso.size(),0);
		
		
		/*
		 for (ii=0; ii<NPixTot; ii++) {
		 if (PixelAcceso[ii]==1) PixelAcceso[ii]=0;
		 }
		 */
		
		
		//for (ii=0; ii<PixelID->size(); ii++){
		//		EnPerPix[ii]=0;
		//	}
		
		
		
		
	} // FINE CICLO SULLE PARTICELLE PRIMARIE
	
	NPixelALL->Draw();
	NPixelALL->SetLineColor(kBlack);
	NPixelALL->Write();
	NPixelGOOD->SetLineColor(kRed);
	NPixelGOOD->Draw("same");
	NPixelGOOD->Write();
	//	NPixelOLD->Draw("same");
	NPixelOLD->Write();
}
