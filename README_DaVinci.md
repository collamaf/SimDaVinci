# Full simulation of Pter beta- probe in Geant4

## HOW TO RUN:
```
cd build
cmake -DGeant4_DIR=$G4INSTALL ../
make
./exampleB1
./exampleB1 {CuDiam (<0->no Cu)} {ZOffs} {FilterFlag} {TBR} {SourceChoice} {x0Scan} {SensorChoice} {PterDiameter} {PterThickness} {SourceDiameter} {SourceThickness} {AbsorberThickness} ../run1.mac
e.g.:
./exampleB1 -5 2 0 1 2 0 1 ../run1.mac
./exampleB1 -5 0.01 0 1 2 0 3 ../run1.mac
```
{all distances/sizes in mm}
Source Choice:
1 - Pointlike Sr
2 - Extended Sr
3 - ExtY
Sensor Choice:
1 - MT9V011
2 - MT9V115
3 - bare 60035 SiPm



## GEOMETRY
- Extended Sr Source ending at Z=0
- Cu collimator on top of Sr source (toggleble)
- Pter Detector starting at Z offset (Z distance is from source surface and possible resin (in case of sensor 1 or up to the sensor in case of sensor 2, even if with filter)
- Sensor resin filter in contact with Pter (towards source)
- Dummy volume for scoring purposes between source (or if presente CU Collimator) and World to score what exits the primary generation
- "infinite" carrier volume behind Pter to simulate mechanical support


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
```
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



## TO DO's

- Manage other sources (gamma)


