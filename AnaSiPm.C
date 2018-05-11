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
	
	// This Macro generates an histogram
	
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	TH1F* NPixel=new TH1F("NPixel","Number of fired pixel per primary particle", 10, 0, 10);
	
	int NPixTot=137*137;
	int debug=0;
	vector<int> PixelAcceso(NPixTot,0.);
	vector<int> FiredPixels;
	//	int PixelAcceso[NPixTot]=0;
	int ii=0;
	vector<int> Counter(nentries,0.);
	
	std::vector<int>::iterator it;                        // Declare an iterator to a vector of int
	
	Long64_t nbytes = 0, nb = 0;
	//	nentries/=10; //to go quickly
	//	nentries=10;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {   // INIZIO CICLO SULLE PARTICELLE PRIMARIE
		
		if (jentry%(nentries/10)==0) cout<<"Sto analizzando l'evento "<<jentry<<" di "<<nentries<<": %= "<<100*jentry/nentries<<endl;
		
		Long64_t ientry = LoadTree(jentry);                 //ientry is an intermediate variable and LoadTree(jentry) sets is value to that of jentry
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		
		
		for (ii=0; ii<PixelID->size(); ii++){               // PixelID->size() coincides with the total number of pixels that have been hit by jth-primary
			if (debug) cout<<" Evento n= "<<jentry<<" ii= "<<ii <<" Pixel acceso: "<< PixelID->at(ii)<<endl;
			it = find(FiredPixels.begin(), FiredPixels.end(),PixelID->at(ii));     // Search for the element  PixelID->at(ii) in the intervall [.begin , .end]
			if (it == FiredPixels.end()) {
				FiredPixels.push_back(PixelID->at(ii));         //Add the PixelID found as last member of the vector FiredPixels
				if (debug) cout<<" Nuovo, aggiungo!"<<endl;
			} else if (debug) cout<<" C'era gia, stica!"<<endl;
			
			//			PixelAcceso[PixelID->at(ii)]=1;
			
			
		}
		
		//		Counter[jentry]=std::count(PixelAcceso.begin(), PixelAcceso.end(), 1);
		Counter[jentry]=FiredPixels.size();        //Put the number of pixels that have been hit by jth-primary into the vector Counter
#if 0
		if (debug) cout<<" Evento n= "<<jentry<<" elenco pixel accesi: "<<endl;
		
		for (ii=0; ii<NPixTot; ii++) {
			if (PixelAcceso[ii]!=0) {
				if (debug) cout<< ii<<endl;
				Counter[jentry]++;
			}
		}
#endif
		if (debug)		cout<<"                     "<<jentry<<" numero pixel accesi: "<<Counter[jentry] <<endl;
		NPixel->Fill(Counter[jentry]);     //Fills the NPixel histogram with the elements of Counter[jentry] (the number of pixels hit)
		
	
		//		PixelAcceso.assign(PixelAcceso.size(),0);
		FiredPixels.clear();               // Clear the vector FiredPixels and put his size to 0; ready for the next iteration of the primary particle
		/*
		 for (ii=0; ii<NPixTot; ii++) {
		 if (PixelAcceso[ii]==1) PixelAcceso[ii]=0;
		 }
		 */
		
		
		
		
		
		
		
		
		
		
		
		
	} // FINE CICLO SULLE PARTICELLE PRIMARIE
	
	NPixel->Draw();
	NPixel->Write();
	
}
