//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  7 12:17:14 2020 by ROOT version 6.18/05
// from TTree B1/physics
// found on file: PTERmc_PDiam6_PDz3_X0_Z0_NoAbs_ExtF_Diam100_Dz70_Set1.root
//////////////////////////////////////////////////////////

#ifndef AnaProbe_h
#define AnaProbe_h

#define NANAPART 3
#define NTHR 4

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class AnaProbe {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Eabs;
   vector<double>  *EabsComp;
   Double_t        InPterTrackN;
   vector<double>  *InPterPart;
   vector<double>  *InPterEn;
   vector<double>  *InPterPrimEn;
   vector<double>  *InPterPrimPart;
   vector<float>   *InPterTime;
   vector<double>  *InPterX;
   vector<double>  *InPterY;
   vector<double>  *InPterZ;
   Double_t        PrePterTrackN;
   vector<double>  *PrePterPart;
   vector<double>  *PrePterEn;
   Double_t        PreProbeTrackN;
   vector<double>  *PreProbePart;
   vector<double>  *PreProbeEn;
   Double_t        PostAbsTrackN;
   vector<double>  *PostAbsPart;
   vector<double>  *PostAbsEn;
   vector<double>  *ExitEne;
   vector<double>  *ExitPart;
   Double_t        SourceX;
   Double_t        SourceY;
   Double_t        SourceZ;
   vector<double>  *SourceCosX;
   vector<double>  *SourceCosY;
   vector<double>  *SourceCosZ;
   vector<double>  *SourceEne;
   vector<double>  *SourcePart;
   vector<double>  *SourceIsotope;
   Int_t           Npmt;
   Int_t           EnterPterFlag;
   vector<double>  *AnnihilationX;
   vector<double>  *AnnihilationY;
   vector<double>  *AnnihilationZ;
   vector<double>  *AnnihilationTime;
   Int_t           Nev;

   // List of branches
   TBranch        *b_Eabs;   //!
   TBranch        *b_EabsComp;   //!
   TBranch        *b_InPterTrackN;   //!
   TBranch        *b_InPterPart;   //!
   TBranch        *b_InPterEn;   //!
   TBranch        *b_InPterPrimEn;   //!
   TBranch        *b_InPterPrimPart;   //!
   TBranch        *b_InPterTime;   //!
   TBranch        *b_InPterX;   //!
   TBranch        *b_InPterY;   //!
   TBranch        *b_InPterZ;   //!
   TBranch        *b_PrePterTrackN;   //!
   TBranch        *b_PrePterPart;   //!
   TBranch        *b_PrePterEn;   //!
   TBranch        *b_PreProbeTrackN;   //!
   TBranch        *b_PreProbePart;   //!
   TBranch        *b_PreProbeEn;   //!
   TBranch        *b_PostAbsTrackN;   //!
   TBranch        *b_PostAbsPart;   //!
   TBranch        *b_PostAbsEn;   //!
   TBranch        *b_ExitEne;   //!
   TBranch        *b_ExitPart;   //!
   TBranch        *b_SourceX;   //!
   TBranch        *b_SourceY;   //!
   TBranch        *b_SourceZ;   //!
   TBranch        *b_SourceCosX;   //!
   TBranch        *b_SourceCosY;   //!
   TBranch        *b_SourceCosZ;   //!
   TBranch        *b_SourceEne;   //!
   TBranch        *b_SourcePart;   //!
   TBranch        *b_SourceIsotope;   //!
   TBranch        *b_Npmt;   //!
   TBranch        *b_EnterPterFlag;   //!
   TBranch        *b_AnnihilationX;   //!
   TBranch        *b_AnnihilationY;   //!
   TBranch        *b_AnnihilationZ;   //!
   TBranch        *b_AnnihilationTime;   //!
   TBranch        *b_Nev;   //!

   AnaProbe(TString filename);
	TString inputFileName;
   virtual ~AnaProbe();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
		TFile* outputfile;
	TString canvDir="CanvDir";
	TDirectory* dirCanvas;
//	TDirectory* dirPart=new TDirectory();
	TDirectory* dirPart[NANAPART];
	TDirectory* dirThr[NTHR];


};

#endif

#ifdef AnaProbe_cxx
AnaProbe::AnaProbe(TString filename) : fChain(0)
{
	inputFileName=filename;
	TFile *file = new TFile(Form("%s.root",filename.Data()));
	TTree *tree = (TTree*)gDirectory->Get("B1");
	if (file!=NULL) outputfile=new TFile(Form("%s_Ana.root",filename.Data()),"RECREATE");
	Init(tree);
	
	dirCanvas = outputfile->mkdir("Canvas");
//	dirPart[0] = outputfile->mkdir("Ele");
//	dirPart[1] = outputfile->mkdir("Pos");
//	dirPart[2] = outputfile->mkdir("Fot");

	
	
//	fVectorNoise->Dump();
}

AnaProbe::~AnaProbe()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaProbe::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaProbe::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnaProbe::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EabsComp = 0;
   InPterPart = 0;
   InPterEn = 0;
   InPterPrimEn = 0;
   InPterPrimPart = 0;
   InPterTime = 0;
   InPterX = 0;
   InPterY = 0;
   InPterZ = 0;
   PrePterPart = 0;
   PrePterEn = 0;
   PreProbePart = 0;
   PreProbeEn = 0;
   PostAbsPart = 0;
   PostAbsEn = 0;
   ExitEne = 0;
   ExitPart = 0;
   SourceCosX = 0;
   SourceCosY = 0;
   SourceCosZ = 0;
   SourceEne = 0;
   SourcePart = 0;
   SourceIsotope = 0;
   AnnihilationX = 0;
   AnnihilationY = 0;
   AnnihilationZ = 0;
   AnnihilationTime = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Eabs", &Eabs, &b_Eabs);
   fChain->SetBranchAddress("EabsComp", &EabsComp, &b_EabsComp);
   fChain->SetBranchAddress("InPterTrackN", &InPterTrackN, &b_InPterTrackN);
   fChain->SetBranchAddress("InPterPart", &InPterPart, &b_InPterPart);
   fChain->SetBranchAddress("InPterEn", &InPterEn, &b_InPterEn);
   fChain->SetBranchAddress("InPterPrimEn", &InPterPrimEn, &b_InPterPrimEn);
   fChain->SetBranchAddress("InPterPrimPart", &InPterPrimPart, &b_InPterPrimPart);
   fChain->SetBranchAddress("InPterTime", &InPterTime, &b_InPterTime);
   fChain->SetBranchAddress("InPterX", &InPterX, &b_InPterX);
   fChain->SetBranchAddress("InPterY", &InPterY, &b_InPterY);
   fChain->SetBranchAddress("InPterZ", &InPterZ, &b_InPterZ);
   fChain->SetBranchAddress("PrePterTrackN", &PrePterTrackN, &b_PrePterTrackN);
   fChain->SetBranchAddress("PrePterPart", &PrePterPart, &b_PrePterPart);
   fChain->SetBranchAddress("PrePterEn", &PrePterEn, &b_PrePterEn);
   fChain->SetBranchAddress("PreProbeTrackN", &PreProbeTrackN, &b_PreProbeTrackN);
   fChain->SetBranchAddress("PreProbePart", &PreProbePart, &b_PreProbePart);
   fChain->SetBranchAddress("PreProbeEn", &PreProbeEn, &b_PreProbeEn);
   fChain->SetBranchAddress("PostAbsTrackN", &PostAbsTrackN, &b_PostAbsTrackN);
   fChain->SetBranchAddress("PostAbsPart", &PostAbsPart, &b_PostAbsPart);
   fChain->SetBranchAddress("PostAbsEn", &PostAbsEn, &b_PostAbsEn);
   fChain->SetBranchAddress("ExitEne", &ExitEne, &b_ExitEne);
   fChain->SetBranchAddress("ExitPart", &ExitPart, &b_ExitPart);
   fChain->SetBranchAddress("SourceX", &SourceX, &b_SourceX);
   fChain->SetBranchAddress("SourceY", &SourceY, &b_SourceY);
   fChain->SetBranchAddress("SourceZ", &SourceZ, &b_SourceZ);
   fChain->SetBranchAddress("SourceCosX", &SourceCosX, &b_SourceCosX);
   fChain->SetBranchAddress("SourceCosY", &SourceCosY, &b_SourceCosY);
   fChain->SetBranchAddress("SourceCosZ", &SourceCosZ, &b_SourceCosZ);
   fChain->SetBranchAddress("SourceEne", &SourceEne, &b_SourceEne);
   fChain->SetBranchAddress("SourcePart", &SourcePart, &b_SourcePart);
   fChain->SetBranchAddress("SourceIsotope", &SourceIsotope, &b_SourceIsotope);
   fChain->SetBranchAddress("Npmt", &Npmt, &b_Npmt);
   fChain->SetBranchAddress("EnterPterFlag", &EnterPterFlag, &b_EnterPterFlag);
   fChain->SetBranchAddress("AnnihilationX", &AnnihilationX, &b_AnnihilationX);
   fChain->SetBranchAddress("AnnihilationY", &AnnihilationY, &b_AnnihilationY);
   fChain->SetBranchAddress("AnnihilationZ", &AnnihilationZ, &b_AnnihilationZ);
   fChain->SetBranchAddress("AnnihilationTime", &AnnihilationTime, &b_AnnihilationTime);
   fChain->SetBranchAddress("Nev", &Nev, &b_Nev);
   Notify();
}

Bool_t AnaProbe::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaProbe::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaProbe::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnaProbe_cxx
