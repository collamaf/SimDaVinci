# Full simulation of Pter beta- probe in Geant4

## HOW TO RUN:
```
cd build
cmake -DGeant4_DIR=$G4INSTALL ../
make
./exampleB1 {flags (see below)} 
```
### Some Use cases:
- Gallium68 Campaign Gemelli measurements (late 2018) with PVC
```
./exampleb1 -GaSet 3 -Source 4 -AbsD 0 -AbsT 5.5 -AbsMat 4  -SourceT 6.4 -Z 5.501  -Vis 1

```
- Gallium68 Campaign Gemelli measurements (late 2018) No PVC
```
./exampleb1 -GaSet 3 -Source 4 -AbsD -1 -AbsT 5.5 -AbsMat 4  -SourceT 6.4 -Z 0  -Vis 1

```

- Electron Efficiency of drop-in NL-probe (Post Gallium68 Campaign Gemelli measurements (late 2018))
```
./exampleb1  -Source 6   -Z 0 -CaseDepth 50 -NPrim 10000000
```

- Gamma Efficiency of drop-in NL-probe (Post Gallium68 Campaign Gemelli measurements (late 2018))
```
./exampleb1  -Source 7   -Z 0 -CaseDepth 50 -NPrim 100000000
```

- Electron Efficiency of "standard" (pen-like) probe (Post Gallium68 Campaign Gemelli measurements (april 2019))
```
./exampleb1  -Source 6  -Z 0 -NPrim 10000000
```

- Gamma Efficiency of  "standard" (pen-like) probe (Post Gallium68 Campaign Gemelli measurements (april 2019))
```
./exampleb1  -Source 7  -Z 0 -NPrim 100000000
```
- Monte Carlo tests with F18 with Robotic Probe (GEN 2020)
```
./exampleb1  -Source 9   -Z 0 -CaseDepth 50 -SourceD 30 -Vis 1
```

- Monte Carlo tests with variable isotope with Standard Probe (FEB 2020)
```
./exampleb1  -Source -3990    -Z 0 -SourceD 20 -SourceT 10  -Vis 1
```

- Gamma Efficiency of drop-in NL-probe in F18 application via external generated file (from GAMOS) (OCT 2020)
```
./exampleb1 -Source 10 -SourceExtFile "sourceFileModP2_123MBq.txt" -Z -30 -CaseDepth 50 -CaseLT 1 -HSBT 2 -NPrim 4523120
./exampleb1 -Source 10 -SourceExtFile "sourceFileModP2_123MBq.txt" -Z -30 -CaseDepth 50 -NPrim 4523120

```

- Gamma Efficiency of drop-in NL-probe in F18 application on classical F18 signal sample (OCT 2020)
```
./exampleb1  -Source -0918 -SourceD 20 -SourceT 10 -Z 0 -CaseDepth 50 -Vis 1
```

- New Container for 18F measurments for gyneco applications (DEC 2020): Needs GaSet3, SourceT gives the height in mm of the fluorine column besides GaSet, SourceD is such that 4cm means a 2cm thick ring around the container. The bottom thickness of the liquid is hard coded and fixed to 1cm for now

```
./exampleb1 -GaSet 3  -Source 11 -AbsD -10 -AbsT 5.5 -AbsMat 4  -Z 5.5 -SourceT 80 -SourceD 40 -Vis 1
./exampleb1 -GaSet 3  -Source 11 -Z 0 -SourceT 80 -SourceD 40 -Vis 1
```

### Misure per capire sonda morelli in vista di sperabile in-vivo (GEN2021): sonda open con scotch davanti su sorgente Sr estesa

- Base
```
./exampleb1  -Source 2  -PterD 6  -PterT 3 -SecondShield 1 -Z 0 -NPrim 1000000 
```
- Spessori Alluminio
```
./exampleb1  -Source 2 -AbsT 0.1 -AbsMat 3 -AbsD 0 -PterD 6  -PterT 3 -Z 0.101 -SecondShield 1  -NPrim 1000000
./exampleb1  -Source 2 -AbsT 0.2 -AbsMat 3 -AbsD 0 -PterD 6  -PterT 3 -Z 0.201 -SecondShield 1  -NPrim 1000000
./exampleb1  -Source 2 -AbsT 0.3 -AbsMat 3 -AbsD 0 -PterD 6  -PterT 3 -Z 0.301 -SecondShield 1  -NPrim 1000000
./exampleb1  -Source 2 -AbsT 0.4 -AbsMat 3 -AbsD 0 -PterD 6  -PterT 3 -Z 0.401 -SecondShield 1  -NPrim 1000000
./exampleb1  -Source 2 -AbsT 0.5 -AbsMat 3 -AbsD 0 -PterD 6  -PterT 3 -Z 0.501 -SecondShield 1  -NPrim 1000000
```
- Spessori Rame

```
./exampleb1  -Source 2 -AbsT 0.3 -AbsMat 1 -AbsD 0 -PterD 6  -PterT 3 -Z 0.301 -SecondShield 1 -NPrim 1000000
./exampleb1  -Source 2 -AbsT 0.6 -AbsMat 1 -AbsD 0 -PterD 6  -PterT 3 -Z 0.601 -SecondShield 1 -NPrim 1000000
./exampleb1  -Source 2 -AbsT 0.9 -AbsMat 1 -AbsD 0 -PterD 6  -PterT 3 -Z 0.901 -SecondShield 1 -NPrim 1000000
```
- Spessori ABS
```
./exampleb1  -Source 2 -AbsT 1.46 -AbsMat 5 -AbsD 0 -PterD 6  -PterT 3 -Z 1.47 -NPrim 1000000 -SecondShield 1 -Label HiDens
./exampleb1  -Source 2 -AbsT 1.67 -AbsMat 5 -AbsD 0 -PterD 6  -PterT 3 -Z 1.68 -NPrim 1000000 -SecondShield 1 -Label HiDens
./exampleb1  -Source 2 -AbsT 1.85 -AbsMat 5 -AbsD 0 -PterD 6  -PterT 3 -Z 1.86 -NPrim 1000000 -SecondShield 1 -Label HiDens
```
- Spessori aggiuntivi bianchi
```
./exampleb1  -Source 2 -AbsT 1.27 -AbsMat 5 -AbsD 0 -PterD 6  -PterT 3 -Z 1.28 -NPrim 1000000 -SecondShield 1 -Label HiDens
./exampleb1  -Source 2 -AbsT 3.13 -AbsMat 5 -AbsD 0 -PterD 6  -PterT 3 -Z 3.14 -NPrim 1000000 -SecondShield 1 -Label HiDens
./exampleb1  -Source 2 -AbsT 3.30 -AbsMat 5 -AbsD 0 -PterD 6  -PterT 3 -Z 3.31 -NPrim 1000000 -SecondShield 1 -Label HiDens
./exampleb1  -Source 2 -AbsT 4.97 -AbsMat 5 -AbsD 0 -PterD 6  -PterT 3 -Z 4.97 -NPrim 1000000 -SecondShield 1 -Label HiDens
```
- Spessori neri (scoperto che sono probabilmente PVC ma con densità 1.4)
```
./exampleb1  -Source 2 -AbsT 2.01 -AbsMat 6 -AbsD 0 -PterD 6  -PterT 3 -Z 2.02 -NPrim 1000000 -SecondShield 1
./exampleb1  -Source 2 -AbsT 2.35 -AbsMat 6 -AbsD 0 -PterD 6  -PterT 3 -Z 2.36 -NPrim 1000000 -SecondShield 1
./exampleb1  -Source 2 -AbsT 3.38 -AbsMat 6 -AbsD 0 -PterD 6  -PterT 3 -Z 3.39 -NPrim 1000000 -SecondShield 1
```

- To obtain efficiency curve:
 With Eabs
```
B1->Draw("SourceEne>>num(100)","Eabs>68","")
B1->Draw("SourceEne>>denom(100)","EnterPterFlag==1","")
num->Sumw2()
num->Divide(denom)
num->Draw("E")
```
 With NPMT
```
B1->Draw("SourceEne>>num(200)","Npmt>55","")
B1->Draw("SourceEne>>denom(200)","EnterPterFlag==1","")
num->Sumw2()
num->Divide(denom)
num->Draw("E")

```



## GEOMETRY
- All sources ending at Z=0 (including dummy)
- DummyExitSorg: after source (with possible air gap in GaSet 3)
- DummyExitAbs: after absorber (if any)
- DummyEnterProbe: before FrontShield


## FLAGS

### GaSet 
To choose the experimental setup; 
- 1: "bare probe", no "catafalco" [default]
- 2: "catafalco" in PVC for liquid source measurements 

### Z
Distance from probe frontshield to source surface [mm, def: 2]

### AbsT
Absorber's thickness [mm, def: 1]

### AbsD
Diameter of the absorber's hole in its center. If <0 no absorber is placed. [mm, def: -1] 

### AbsMat
Absorber's material:
- 1: Cu [def]
- 2: Pb
- 3: Alu
- 4: PVC
- 5: ABS
- 6: Strange black plastic with middle density (1.3 g/cm3)
- 7: Nucleomed 3d printed material (1.21 g/cm3)

### Source
Source type: 
- 1: Pointlike Sr [def]
- 2: Extended Sr
- 3: ExtY
- 4: ExtGa
- 5: 511KeV gamma sphere
- 6: FlatEle
- 7: FlatGamma sphere
- 8: ExtCu67 (vol for MC studies, not exp. meas.)
- 9: ExtF18 (vol for MC studies, not exp. meas.)
- 10: ExternaFile (typically no lab source)
- 11: 18F container around GaSet3
- <0: Generic Ion, in forma ZZAA (eg 3990 is Y)

### X 
Horizontal offset of the probe [mm, def: 0]

### SourceD
Source Diameter for variable size souces [mm, def: 10]

### SourceT 
Source Thickness for variable size souces [mm, def: 7]

### PterD
PTER Diameter [mm, def: 6]

### PterT 
PTER Thickness [mm, def: 3]

### CaseLT
Lateral thickness of robotic probe case [mm, def: 1]

### CaseBT
Bottom thickness of robotic probe case [mm, def: 20]

### HSLT
Lateral thickness of horseshoe structure inside robotic probe case [mm, def: 1]

### HBLT
Bottom thickness of horseshoe structure inside robotic probe case [mm, def: 2]

### CaseDepth 
Length of the probe case [mm, def: -155]:
- >0: robotic probe
- ==0: "bare probe" (just PTER+Delrin+PVC)
- <0: open surgery probe

### HSMat
Material of robotic probe case. From inside out
- 1: PVC+Pb+PVC [default]
- 2: AIR+Pb+AIR
- 3: AIR+AIR+AIR

### AppMat
Material of the "catafalco":
- 1: AIR [default]
- 2: Pb (to kill almost everything, maybe..)


## PHYSICS
Process				Type		SubType
RadioActiveDecay	6			210
eIoni				2			2
eBrem				2			3

CPU TIMES NEEDED FOR 1e5 PRIMARIES:

## SCORING:
- Primary paricles generated (PrimAction);
- Decay products (StackingAction);
- What exits source (!->DummyExitSorg) (StepAction);
- What exit Absorber (if any) (Abs->DummyExitAbs) (StepAction);
- What enters Probe (DummyEnterProbe->FrontShield) (StepAction);
- What enters Pter (!->Pter) (StepAction)

## OUTPUT:
A root file named MCsondaGEANT_{Relevant sim info}.root is created, in which on an event (i.e. a primary particle) by event basis it is stored:

### SOURCE vector (one entry per primary particle, only for first 100k primaries if more primaries are requested):
- 0: AllX: X coordinate of primary particle [mm];
- 1: AllY: Y coordinate of primary particle [mm];
- 2: AllZ: Z coordinate of primary particle [mm];
- 3: AllCosX[1+]: X directive cosine of produced particle;
- 4: AllCosY[1+]: Y directive cosine of produced particle;
- 5: AllCosZ[1+]: Z directive cosine of produced particle;
- 6: AllEne[1+]: kinetic energy of produced particle [keV];
- 7: AllPart[1+]: PID of produced particle;
- 8: AllIsotope[1+]: isotope of primary particle (0=Sr, 1=Y ??, ParentID - 1);
- 9: ExitX: X coordinate of primary particle exiting the source volume (-> DummyExitSorg) [mm];
- 10: ExitY: Y coordinate of primary particle exiting the source volume [mm];
- 11: ExitZ: Z coordinate of primary particle exiting the source volume [mm];
- 12: ExitCosX[2+]: X directive cosine of primary particle exiting the source volume;
- 13: ExitCosY[2+]: Y directive cosine of primary particle exiting the source volume;
- 14: ExitCosZ[2+]: Z directive cosine of primary particle exiting the source volume;
- 15: ExitEne[2+]: kinetic energy of primary particle exiting the source volume [keV];
- 16: ExitPart[2+]: kind of primary particle (11=e-, -11=e+, 22=gamma, 13=mu-...) exiting the source volume;
- 17: ExitParentID[2+]: partent-id of particle exiting the source
- 18: ExitProcess[2+]: process that created the particles that exits the source (see table above)
- 19: ExitTrackN: number of different tracks exiting the source per event
- 20: AnnihilationX: X coord of annihilation point of beta+ source particle [mm]
- 21: AnnihilationY: Y coord of annihilation point of beta+ source particle [mm]
- 22: AnnihilationZ: Z coord of annihilation point of beta+ source particle [mm]

### B1 vector (one entry per primary particle that gives a >0 energy deposition):
- 0: Eabs: energy absorbed in Pter [keV];
- 1: EAbsComp[3]: vector containing energy absorbed in Pter [keV] due to Electrons (EAbsComp[0]), Positrons (EAbsComp[1]) and photons (EAbsComp[2])
- 2: InPterTrackN: number of hits inside Pter (it's the length of the following vector);
- 3: InPterPart[InPterTrackN]: kind of particle of hit inside Pter;
- 4: InPterEn[InPterTrackN]: energy deposit of single hit of particle inside Pter;
- 5: InPterPrimEn[InPterTrackN]: energy of the primary particle (decay product, or fresh primary e/gamma) that origined the hit of particle inside Pter [keV];
- 6: InPterPrimPart[InPterTrackN]: PID of the primary particle (decay product, or fresh primary e/gamma) that origined the hit of particle inside Pter [keV];
- 7: InPterTime[InPterTrackN]: time of interaction of hit inside Pter [ns] (To be really undersood);
- 8: InPterX[InPterTrackN]: X position of hit inside Pter [mm];
- 9: InPterY[InPterTrackN]: Y position of hit inside Pter [mm];
- 10: InPterZ[InPterTrackN]: Z position of hit inside Pter [mm];
- 11: PrePterTrackN: number of tracks entering Pter (from wherever) per primary (it's the length of the following vector);
- 12: PrePterPart[PrePterTrackN]: kind of particle of each track entering Pter (from wherever);
- 13: PrePterEn[PrePterTrackN]: kinetic energy of particle of each tracks entering Pter (from wherever) [keV];
- 14: PrePterX[PrePterTrackN]: X position of particle of each tracks entering Pter (from wherever) [mm];
- 15: PrePterY[PrePterTrackN]: Y position of particle of each tracks entering Pter (from wherever) [mm];
- 16: PrePterZ[PrePterTrackN]: Z position of particle of each tracks entering Pter (from wherever) [mm];
- 17: PreProbeTrackN: number of tracks entering Probe (DummyEnterProbe->FrontShield) per primary (it's the length of the following vector);
- 18: PreProbePart[PreProbeTrackN]: kind of particle of each track entering Probe (DummyEnterProbe->FrontShield);
- 19: PreProbeEn[PreProbeTrackN]: kinetic energy of particle of each tracks entering Probe (DummyEnterProbe->FrontShield) [keV];
- 20: PostAbsTrackN: number of tracks exiting Absorber (Abs->DummyExitAbs) (if any, otherwise 0) per primary (it's the length of the following vector);
- 21: PostAbsPart[PrePterTrackN]: kind of particle of each track entering  Absorber (Abs->DummyExitAbs) (if any, otherwise empty);
- 22: PostAbsEn[PrePterTrackN]: kinetic energy of particle of each tracks entering  Absorber (Abs->DummyExitAbs) (if any, otherwise empty)  [keV];
- 23:  ExitEne[2+]: kinetic energy of primary particle exiting the source volume (-> DummyExitSorg) [keV];
- 24: SourceX: X coordinate of primary particle (isotope) giving a signal in Pter [mm];
- 25: SourceY: Y coordinate of primary particle (isotope) giving a signal in Pter [mm];
- 26: SourceZ: Z coordinate of primary particle (isotope) giving a signal in Pter [mm];
- 27: SourceCosX[2]: X directive cosine of decay electron(s) giving a signal in Pter;
- 28: SourceCosY[2]: Y directive cosine of  decay electron(s) giving a signal in Pter;
- 29: SourceCosZ[2]: Z directive cosine of decay electron(s) giving a signal in Pter;
- 30: SourceEne[2]: kinetic energy of  decay product(s)  giving a signal in Pter [keV];
- 31: SourcePart[2]: PID of  decay product(s)  giving a signal in Pter [keV];
- 32: SourceIsotope: isotope of primary particle (0=Sr, 1=Y) giving a signal in Pter;
- 33: Npmt: number of optical photons entering SiPm per event;
- 34: EnterPterFlag: bool flag to save if the primary particle had a secondary track entering PTER;
- 35: AnnihilationX: X coord of annihilation point of beta+ source particle [mm]
- 36: AnnihilationY: Y coord of annihilation point of beta+ source particle [mm]
- 37: AnnihilationZ: Z coord of annihilation point of beta+ source particle [mm]
- 38: AnnihilationTime: time of annihilation point of beta+ source particle [ns]
- 39: Nev: storing number of events generated


```


file=$(ls -t Primaries_*.dat | head -n1); tail -f $file


```
ELE
POS
FOT



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
- Added flag PosAbsorber to choose absorber position in -GaSet2; if -PosAbs == 1 it will be placed in the hole just near the source, that has 2mm depth, while if -PosAbs == 2 it will be placed  in the hole that has 6mm depth.
- The Absorber thickness could be choosen by AbsT flag while the diameter is fixed in both PosAbsorber cases.
- The flag AbsD must be >=0.


2018.07.17 by MorettiR
- Added opportunity to select placement position of the absorber; the PosAbsorber now is used only to fix absorber diameter, in fact this depends in which hole I'll place the absorber.
- The -Z flag is used to choice the absorber's placement in -GaSet2 case.
- Example of a line comad for the optimal choice of the absorber type: ./exampleb1 -GaSet 2 -CaseDepth -40 -Source 4 -AbsD 0 -AbsT 1 -AbsMat 2 -AppMat 1 -Z 7 -PosAbs 2 ../run1.mac
- Stepping Action fixed when fCudiam>=0.

2018.07.17 by MorettiR
- Fixed bug in geometry relative to Dummy placement.
- Addition of another solid Dummy2 for stepping action in -GaSet 2 case.
- Example of a line comad for the optimal choice of the absorber type: ./exampleb1 -GaSet 2 -CaseDepth -40 -Source 4 -AbsD 0 -AbsT 1 -AbsMat 2 -AppMat 1 -Z 0.5 -PosAbs 1 -PterD 5 -PterT 3 ../run1.mac

2018.07.19 by MorettiR
- File name fixed.

2018.07.20 by Collamaf + MorettiR
- Fixed scoring macro
- Fixed source thickness in PrimGenAct
- Added keep event for interesting annihilation events

2018.07.20 by MorettiR
- Addition of annihilation coordinates of positrons in the .root file

2018.07.25 by MorettiR
- Addition of a leaf in B1 branch about energy of particles arriving on Absorber from source  called "PreDummy2En"
- Addition of a leaf in B1 branch about energy of particles arriving on Dummy2 from source  called "InDummy2En"
- Change in stepping action's condition relative to positrons annihilation region; now I'm looking annihilations in all over the world.

2018.07.26 by MorettiR
- New leafs added in B1 and Source Branches 

2018.09.19 by MorettiR
- The position of the Probe in GaSet2 has been modified; now the top of the probe is placed at 5.5mm from the source (before at 8mm).
This because the configuration used at Gemelli's hospital in which the probe was placed on the absrober's surface (We used absorber of 5.5mm thickness). In the configuration without the absorber the probe was placed at the same distance of 5.5mm from the source.
- Fixed the position of the probe, now is take in consideration the dummy2's thickness.
- SourceDiameter now fixed at 10.*mm and SourceThickness at 7.*mm

2018.10.2 by collamaf
- Fixed error in StackingAction: we were writing info on produced particles only if PDGcode>0, thus cutting all anti particles

2018.11.9 by MorettiR
- README file uploaded with flags specifications.
- Fixing of line 1639 and 1640 of DetectorConstruction.cc .

2018.11.11 by MorettiR
- Now the position of the probe's head ( the frontshield ) varies automatically with the placement of the absorber in GaSet 2.  
  N.B. in this way you have always to give AbsT and Z even if you do not want to place the absorber, for more details take a look to the flags section in this ReadMe file. 
- Name of the file with no abs in GaSet 2 has been modified. Now there is inside the name the distance of the probe from the source.


2018.11.14 by MorettiR
- A new leaf called "AnnihilationTime" has been added to B1 vector. This gives you the time elapsed since the positron's creation and his annihilation. 
- Now we can count the number of hits in the SiPM thanks to the new leaf "EabsSiPM" in B1 Vector.

2018.11.23 by MorettiR
- Added new configuration GaSet3; in this one the GaContainer is obtained by PVC and not printed by 3D printer. 
- File name modified to take into account GaSet3 configuration.
- New stepping action condition in source's branch (to be cheked). Take a look to line 144 an line 76 in stepping action to see changes. 

2018.12.03 by MorettiR
- Table added in GaSet1 configuration.

2018.12.05 by MorettiR
- Added new flag for absorber placement -ZAbs. This corresponds to the center of the absorber on z axes.
- Now in GaSet 2 and GaSet3, the flag -Z corresponds to the distance of the frontshield from the source on z axes (-Z must be always >= -AbsT when the absorber is placed to avoid superposition of volumes).
- Fixed superposistion of absorber in GaSet1
- Added study of Catafalco's material in GaSet3

2018.12.07 by MorettiR (morning)
- Fixed GaSet3 dimension of catafalco's Gacontainer.

2018.12.07 by MorettiR (afternoon)
- Fixed enumeration in B1RunAction (take a look also to B1EventAction).

2018.12.08 by MorettiR
- Fixed line 99 of B1RunAction.

2018.12.11 by MorettiR
- Fixed error in geometry relative to the source.
- Fixed error in PrimaryGeneratorAction relativo to GaSet 3.

2018.12.12 by MorettiR
- Fixed error in geometry relative to the dimension of the dummy2 in GaSet3.
- Fixed error in geometry relative to the dimension of the absorber in GaSet3.

2018.12.12 by MorettiR (afternoon)
- Now you can choose di thickness of the source in GaSet3.

2018.12.12 by collamaf
- Added flat sources (ele-6, gamma-7) for efficiency studies

2019.01.09 by collamaf
- First attempt of deep reorganization and cleaning of DetConst: seems to work fine at least up to GaSet 1. Please note that in this case the absorber is always ontop of source (no need of AbsorberPos) and the Z distance of the probe is always from source surface

2019.01.10 by collamaf
- Dummy reorganisation seems to work at least for Gaset 1 and 3. Need to fix DummyExitSource for GaSet 3, but first need to fix ExtGa z position in this case to start from the bottom of the hole

2019.01.14 by collamaf
- Finally got rid of all old "dummies". New 3 dummy scheme seems to work for all 3 gaset!
- Removed all references to FilePrim

2019.01.16 by collamaf
- Deep Cleaning of scoring, root file, PrimAction etc
- Now primary particle position in GaSet 2/3 is based on "Source" volume in DetConst
- Still need to clean OutFIleName creation in main

2019.01.17 by collamaf
- Cleaned OutFileName generation in main
- Fixed InPterPrimPart (was same as SourcePart)
- Updated readme about scoring

2019.04.01 by collamaf
- Updated readme about Efficiency (also for stanard pen like probe)

2019.10.02 by collamaf
- Added sourcechoice 8: Cu67 "ex-vivo" simple sample
- Added scoring of Vertex XYZ of particles exiting the source

2020.02.07 by collamaf
- Add extended F18 source (like Cu), source 9
- Fix StackAct condition to intercept decay product with new RadioactiveDecayBase process (should however be retro compatible)
- Now ScintFlag is passed to PhysList also (not only DetConst), since the mere process somehow overshadowed the RadioactiveDecay one preventing the StackAction scoring
- Add "AnaProbe" analysis macro
- Add ExitPart scoring also in B1 tree

2020.02.17 by collamaf
- Add new <0 SourceSelect, to choose generic radioactive ion: -SourceSelect -ZZAA (has the same source geometry than SourceSelect 8-9: simple active water voulme (manageble via SourceD SourceT) 

2020.02.25 by collamaf
- Fix physics list: opt4 and proper cuts (maybe excessive too..)

2020.10.08 by collamaf
- Add possibility to use external HepMC generator for 18F PET/CT simulations performed with GAMOS. SourceSelect=10, no source volume is built. Scoring of this kind of primary particle has to be fixed.

2020.10.16 by collamaf
- Add scoring of X Y Z position for PrePter

2020.12.23 by collamaf
- Add new source (11) for F18 container, to be used only around GaSet3. Seems to work and be retrocompatible with everything, but is a little too much accrocchiated. Probably time to perform a deep bottom-up reorganization of the whole simulation
- Restore MT capability and add automatic merging of root files

2021.01.19 by collamaf
- Go back to manual merging of root files in MT since the G4 way does not work on vectors

2021.01.19 by collamaf
- Add possibility to add a "SecondShield", i.e. 10um of scotch in front. Command line argument to activate it

2021.02.01 by collamaf
- Add absorber material 6 for black plastic, maybe PVC but with slightly higher density (1.4 vs 1.3))
- Fix missing water around GaContainer in 18F studies (source 11)

2021.02.08 by collamaf
- Massive cleaning of DetConst and README

2021.02.09 by collamaf
- Don't place dummyExitSorg with external file source

2021.02.09 by collamaf
- General cleaning of REAMDE
- Fix some wrong pos and sizes (dummy volumes..)
- Move to boolean solid for external layer of 18F (source 11)
- Add scoring of source particle region
- Add NUCLEOMED plastic material for absorbers

## TO DO's

- capire perche quando si accende la scintillazione poi guardando il verbose dell'evento al posto di RadioactiveDecay compare sempre "Scintillation" (ma lo scoring della sorgente lo risonosce uguale..)
