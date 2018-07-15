# Full simulation of Pter beta- probe in Geant4

## HOW TO RUN:
```
cd build
cmake -DGeant4_DIR=$G4INSTALL ../
make
./exampleB1
./exampleB1 {1-CuDiam (<0->no Cu)} {2-ZOffs} {3-FilterFlag} {4-TBR} {5-SourceChoice} {6-x0Scan} {7-Absorber Material} {8-PterDiameter} {9-PterThickness} {10-SourceDiameter} {11-SourceThickness} {12-AbsorberThickness} ../run1.mac
e.g.:

./exampleB1 -AbsD -1 -SourceD 6 -SourceT 3 -X 0 -PterD 6 -PterT 4 -Z 1 -TBR 0 -Source 4 ../run1.mac

 
 
 -AbsD= Absorber Diameter, -AbsT= Absorber Thickness , -AbsMat= Absorber's Material , -Z distance from origin to FrontsShield
 
 -GaSet = 1 for Old Configuration, 2 for 3D-Printed Apparatus 
 
 
 
```
{all distances/sizes in mm}
Source Choice:
1 - Pointlike Sr
2 - Extended Sr
3 - ExtY
4 - ExtGa



## GEOMETRY
- Extended Sr Source ending at Z=0
- Cu collimator on top of Sr source (toggleble)
- Pter Detector starting at Z offset (Z distance is from source surface to Front-Shield )
- Front-Shield in contact with Pter (towards source)
- Dummy volume for scoring purposes between source (or if presente CU Collimator) and World to score what exits the primary generation


## PHYSICS
Process				Type		SubType
RadioActiveDecay	6			210
eIoni				2			2
eBrem				2			3

CPU TIMES NEEDED FOR 1e5 PRIMARIES:


## OUTPUT:
The usual PrimariesX{}_Z{}_CuD{}_Fil{}_TBR{}{_Sr}.dat is created to keep track of the progress
A root file named MCsondaGEANT_Z{XX}.root is created, reporting the Z offset value, in which on an event (i.e. a primary particle) by event basis it is stored:
### SOURCE vector (one entry per primary particle, only for first 100k primaries if more primaries are requested):
- AllX: X coordinate of primary particle [mm];
- AllY: Y coordinate of primary particle [mm];
- AllZ: Z coordinate of primary particle [mm];
- AllCosX[2]: X directive cosine of produced electron;
- AllCosY[2]: Y directive cosine of produced electron;
- AllCosZ[2]: Z directive cosine of produced electron;
- AllEne[2]: kinetic energy of produced electron [keV];
- AllIsotope[2]: isotope of primary particle (0=Sr, 1=Y);
- ExitX: X coordinate of primary particle exiting the source volume [mm];
- ExitY: Y coordinate of primary particle exiting the source volume [mm];
- ExitZ: Z coordinate of primary particle exiting the source volume [mm];
- ExitCosX[2+]: X directive cosine of primary particle exiting the source volume;
- ExitCosY[2+]: Y directive cosine of primary particle exiting the source volume;
- ExitCosZ[2+]: Z directive cosine of primary particle exiting the source volume;
- ExitEne[2+]: kinetic energy of primary particle exiting the source volume [keV];
- ExitPart[2+]: kind of primary particle (11=e-, -11=e+, 22=gamma, 13=mu-...) exiting the source volume;
- ExitParentID[2+]: partent-id of particle exiting the source
- ExitProcess[2+]: process that created the particles that exits the source (see table above)
- ExitTrackN: number of different tracks exiting the source per event

### B1 vector (one entry per primary particle that gives a >0 energy deposition):
- Eabs: energy absorbed in Pter [keV];
- EAbsComp[2]: vector containing energy absorbed in Pter [keV] due to Sr (comp 1) and to Y (comp 2)
- PrePterTrackN: number of tracks entering Pter per primary (from front Alluminum) (it's the length of the following vector);
- PrePterPart[PrePterTrackN]: kind of particle of each track entering Pter(from front resin);
- PrePterEn[PrePterTrackN]: kinetic energy of particle of each tracks entering Pter (from front resin) [keV];
- InPterTrackN: number of hits inside Pter (it's the length of the following vector);
- InPterPart[InPterTrackN]: kind of particle of hit inside Pter;
- InPterEn[InPterTrackN]: energy deposit of single hit of particle inside Pter;
- InPterEnPrim[InPterTrackN]: energy of the primary particle that origined the hit of particle inside Pter [keV];
- InPterTime[InPterTrackN]: time of interaction of hit inside Pter [ns] (To be really undersood);
- InPterX[InPterTrackN]: X position of hit inside Pter [mm];
- InPterY[InPterTrackN]: Y position of hit inside Pter [mm];
- InPterZ[InPterTrackN]: Z position of hit inside Pter [mm];
- PixelID[InPterTrackN]: number of pixel in which the hit occurred (from 1 to NpixMax);
- PixXPos[InPterTrackN]: x position of the pixel in which the hit occurred [mm];
- PixYPos[InPterTrackN]: y position of the pixel in which the hit occurred [mm];
- SourceX: X coordinate of primary particle (isotope) giving a signal in Pter [mm];
- SourceY: Y coordinate of primary particle (isotope) giving a signal in Pter [mm];
- SourceZ: Z coordinate of primary particle (isotope) giving a signal in Pter [mm];
- SourceCosX[2]: X directive cosine of decay electron(s) giving a signal in Pter;
- SourceCosY[2]: Y directive cosine of  decay electron(s) giving a signal in Pter;
- SourceCosZ[2]: Z directive cosine of decay electron(s) giving a signal in Pter;
- SourceEne[2]: kinetic energy of  decay electron(s)  giving a signal in Pter [keV];
- SourceIsotope: isotope of primary particle (0=Sr, 1=Y) giving a signal in Pter;
- Nev: storing number of events generated


```
Source->Draw("ExitEne","ExitPart==11&&ExitProcess==6")

file=$(ls -t Primaries_*.dat | head -n1); tail -f $file


```
ELE
POS
FOT


## TO MERGE IN CASE OF MT

TString nomefile="PTERmc_PDiam6_PDz2_NoAbs_CaseDepth40_CaseLT1_CaseBT5_HSLT4_HSBT3_HSMat3_X0_Z2_Sphere511";
TChain * chain = new TChain("B1")
chain->Add(Form("%s_t*.root",nomefile.Data()))
TChain * chain2 = new TChain("Source")
chain2->Add(Form("%s_t*.root",nomefile.Data()))
TFile *file = TFile::Open(Form("%s.root",nomefile.Data()),"RECREATE");
chain->CloneTree(-1,"fast");
chain2->CloneTree(-1,"fast");
file->Write();


## TO ANALYZE HOUSING BACK MATERIAL

_file2->cd()
B1->Draw("EabsComp[2]","EabsComp[2]>0")
B1->Draw("EabsComp[0]","EabsComp[0]>0")
_file1->cd()
B1->Draw("EabsComp[2]","EabsComp[2]>0")
B1->Draw("EabsComp[0]","EabsComp[0]>0")
_file0->cd()
B1->Draw("EabsComp[2]","EabsComp[2]>0")
B1->Draw("EabsComp[0]","EabsComp[0]>0")

_file2->cd()
B1->Draw("PrePterEn","PrePterPart==22")
B1->Draw("PrePterEn","PrePterPart==11")
_file1->cd()
B1->Draw("PrePterEn","PrePterPart==22")
B1->Draw("PrePterEn","PrePterPart==11")
_file0->cd()
B1->Draw("PrePterEn","PrePterPart==22")
B1->Draw("PrePterEn","PrePterPart==11")




to see energy spectrum of electrons created by Sr/Y that exit the source

Per disegnare contributi Sr e Y:
```
B1->Draw("Eabs")
B1->SetLineColor(kBlue)
B1->Draw("InPterEnSr","","same")
B1->SetLineColor(kRed)
B1->Draw("InPterEnY","","same")
````

## CHANGELOG
2017.12.1 by collamaf
- Try to fix problem of primary particles double counting by putting a check in Stepping Action. Seems to reduce by about ~9% the number of exiting particles (NoCudZ2 test)

2017.12.4 by collamaf
- Added reset of pass counter in stacking action. Seems to increase back of about another 5% the number of exiting particles. now we have about 96% of the beginning. Verified this should be the right approach.

2017.12.12 by collamaf
- Corrected z_resin:  no need to add 1mm of copper since distance is always from source top

2018.01.17 by collamaf
- Deep reorganization of DetConstr, much clearer now
- Code extended to both sensors smoothly (new argument to be passed by terminal, default sensor 1).
- Fixed ExtY source problem

2018.01.30 by collamaf
- Quite deep reorganization of output root file, cleaned a little in ordering and name (and updated readme)
- Changed "ExtY" to "ExtY" in output file naming

2018.04.23 by collamaf
- Now resin is always present in front of Pter, if flag not selected made of air (useful for scoring)
- Added double crossing check also for particles entering Pter

2018.04.26 by collamaf
- Added InPterEnPrim to bring primary particle info to Riduzione and DataAnalysis

2018.05.7 by collamaf
- first implementation of storage of time of interection. Still not clear which time to save...
- Introduced possibility to simulate bare SiPm (assumed to be a particular version of Pter detector): SensorChoice=3


2018.05.31 by collamaf
- Now arguments are taken with labels, not necessary to give them all! If a macro is provided no visualization is init.. cool! (removed useless FilterFlag)

2018.06.05 by collamaf
- Added structure to classify energy release due to e-/e+/gamma

2018.06.11 by MorettiR
- Added Probe's Casing in Geometry with parameters from line comand.

2018.06.11 by collamaf
- Fixed spotted error in fSourceSelect==4 condition in PrimAct
- Added fSourceSelect==5 to generate a sphere (now with R=10cm) around origin of isotropic 511keV gammas to simulate the "far bacgkround" in a Ga68 environment
- Simplification of materials for the probe laparo housing: 3 materials: "Ext" for external housing, "Metal" for metal covering inside and "Inner" for inner part (maybe air)

2018.06.13 by MorettiR
- Probe Casing modified with boolean solids and possibility to choose the presence/absence of the casing, from line comand, simply putting respectively the CaseDepth > 0 or <=0 ( this is true only in the case -GaSet 1 ). 

2018.06.13 by collamaf
- Now when requesting Sphere Source the sphere is centered on the center of the probe case to gain in isotropy
- Changed "entering PTER" condition: now is "from Not PTER to PTER" to consider also tracks from behind (before was "from frontshield")
- Added in output file name also info about the probe case

2018.06.14 by collamaf
- Substituted "Midle" with "Middle"
- Fix sphere source offset when no case is requested

2018.06.20 by collamaf
- Added new argument from command line to choose material for inner probe case: "HSMat": 1 (default) is Pb+plastic, 2 is Pb+Air, 3 is Air+Air
- Added possibility to run in Multi Thread

2018.06.21 by collamaf
- Added scintillation for PTER (with ScintFlag from terminal, default NO). Actually addedd all Optical Physics but disabled for now Cerenkov directly in PhysList
- Added SiPm volume behind pter to score optical photons entering it
- Added scoring of "Npmt" in root file and "_Scint" to filename
- Now the progress status is printed on screen (since on Primaries file sometimes does not work) and every 10% of evts
- Changed condition to score primary decay product: added requeste StepN==0 to avoid double counting particles interacting via optical processes. Should not have any other undesired effect
- Excluded OptFot from InPter scoring to avoid huge size of root file
- Changed X and Z in filename to be in 0.1mm

2018.06.22 by collamaf
- Corrected error in FrontShield Thickness: was 5um is now 15

2018.06.27 by MorettiR
- Boolean Geometry Implemented for SiPM + TopCase solids

2018.07.06 by MorettiR
- Added G4ScoringManager in exampleB1.cc line 263 (Problem to be fixed).
- Addition of "GaSetting" flag to choiche experimental setup in Ga-Source Case.
- Addition of the 3D-Printed expreimental setup.
- In -GaSet 2  -CaseDepth>0 case, every distances are setted and can't be modified by line comand exept the ProbeCase Depth.
- In -GaSet 2 case exists only 2 configurations: if CaseDepth>0 we have the one with ProbeCase while if CaseDepth<0 we have the classic beta minus probe (that is an aluminum case behind the Pter).
- In -GaSet 2 case if CaseDepth<0 his absolute value correspondes to the aluminum case's length.
- In -GaSet 1 case we have the usual simulation (the one before this upgrade).
- To must be add, 3D-Printed expreimental setup's material.

2018.07.15 by MorettiR
- Added Absorber in -GaSet 2.
- Added flag PosAbsorber to choice absorber position in -GaSet2; if PosAbsorber == 1 it will be placed in the hole just near the source, that had 2mm depth, while if PosAbsorber == 2 it will be placed  in the hole that had 6mm depth.
- The Absorber thickness could be choosen by AbsT flag while the diameter is fixed in both -PosAbsorber cases.
- The flag AbsD must be >=0.


## TO DO's

- flag per non piazzare la struttura di supporto se uno dei valori passati da terminale a riguardo è negativo
- Sistemare l'overlap fra il volume SiPm e il tappo presente dietro il PTER quando c'è il case
- capire perche quando si accende la scintillazione poi guardando il verbose dell'evento al posto di RadioactiveDecay compare sempre "Scintillation" (ma lo scoring della sorgente lo risonosce uguale..)

