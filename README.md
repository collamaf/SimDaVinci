# Full simulation of CMOS beta- probe in Geant4

## HOW TO RUN:
```
cd build
cmake -DGeant4_DIR=$G4INSTALL ../
make
./exampleB1
./exampleB1 {CuDiam (<0->no Cu)} {ZOffs} {FilterFlag} {TBR} {SourceChoice} {x0Scan} {SensorChoice} ../run1.mac
e.g.:
./exampleB1 5 2 0 1 2 0 1 ../run1.mac
```
{all distances/sizes in mm}
Source Choice:
1 - Pointlike Sr
2 - Extended Sr
3 - DOTA
Sensor Choice:
1 - MT9V011
2 - MT9V115

## GEOMETRY
- Extended Sr Source ending at Z=0
- Cu collimator on top of Sr source (toggleble)
- CMOS Detector starting at Z offset (Z distance is from source surface and possible resin (in case of sensor 1 or up to the sensor in case of sensor 2, even if with filter)
- Sensor resin filter in contact with CMOS (towards source)
- Dummy volume for scoring purposes between source (or if presente CU Collimator) and World to score what exits the primary generation
- "infinite" carrier volume behind CMOS to simulate mechanical support


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
- XAll: X coordinate of primary particle [mm];
- YAll: Y coordinate of primary particle [mm];
- ZAll: Z coordinate of primary particle [mm];
- AllCosX[2]: X directive cosine of produced electron;
- AllCosY[2]: Y directive cosine of produced electron;
- AllCosZ[2]: Z directive cosine of produced electron;
- AllEne[2]: kinetic energy of produced electron;
- AllIsotope[2]: isotope of primary particle (0=Sr, 1=Y);
- ExitX: X coordinate of primary particle exiting the source volume [mm];
- ExitY: Y coordinate of primary particle exiting the source volume [mm];
- ExitZ: Z coordinate of primary particle exiting the source volume [mm];
- CosXExit: X directive cosine of primary particle exiting the source volume;
- CosYExit: Y directive cosine of primary particle exiting the source volume;
- CosZExit: Z directive cosine of primary particle exiting the source volume;
- EneExit: kinetic energy of primary particle exiting the source volume;
- PartExit: kind of primary particle (11=e-, -11=e+, 22=gamma, 13=mu-...) exiting the source volume;

### B1 vector (one entry per primary particle that gives a >0 energy deposition):
- Eabs: energy absorbed by Cmos [keV];
- PreCmosTrackN: number of tracks entering Cmos per primary (from front shield if present otherwise from world) (it's the length of the following vector);
- PreCmosPart[PreCmosTrackN]: kind of particle of each track entering Cmos (from front shield if present otherwise from world);
- PreCmosEn[PreCmosTrackN]: kinetic energy of particle of each tracks entering Cmos (from front shield if present otherwise from world) [keV];
- InCmosTrackN: number of hits inside Cmos (it's the length of the following vector);
- InCmosPart[InCmosTrackN]: kind of particle of hit inside Cmos;
- InCmosEn[InCmosTrackN]: energy deposit of single hit of particle of hit inside Cmos;
- InCmosX[InCmosTrackN]: X position of hit inside Cmos [mm];
- InCmosY[InCmosTrackN]: Y position of hit inside Cmos [mm];
- InCmosZ[InCmosTrackN]: Z position of hit inside Cmos [mm];
- InCmosSr[InCmosTrackN]: energy released by sons of Sr [keV];
- InCmosEnY[]InCmosTrackN: energy released by sons of Y [keV];
- PixelID[InCmosTrackN]: number of pixel in which the hit occurred (from 1 to NpixMax);
- EDepInPixel[InCmosTrackN]: energy deposited in single pixel [keV]; *** same as InCmosEn
- PixXPos[InCmosTrackN]: x position of the pixel in which the hit occurred [mm];
- PixYPos[InCmosTrackN]: y position of the pixel in which the hit occurred [mm];
- SourceX: X coordinate of primary particle (isotope) giving a signal in Cmos [mm];
- SourceY: Y coordinate of primary particle (isotope) giving a signal in Cmos [mm];
- SourceZ: Z coordinate of primary particle (isotope) giving a signal in Cmos [mm];
- SourceCosX: X directive cosine of decay electron(s) giving a signal in Cmos;
- SourceCosY: Y directive cosine of  decay electron(s) giving a signal in Cmos;
- SourceCosZ: Z directive cosine of decay electron(s) giving a signal in Cmos;
- SourceEne: kinetic energy of  decay electron(s)  giving a signal in Cmos [keV];
- SourceIsotope: isotope of primary particle (0=Sr, 1=Y) giving a signal in Cmos;
- Nev: storing number of events generated
- OptPhotNo: number of optical photons passing from Cmos to SiPm per event;


```
Source->Draw("ExitEne","ExitPart==11&&ExitProcess==6")
```
to see energy spectrum of electrons created by Sr/Y that exit the source

Per disegnare contributi Sr e Y:
```
B1->Draw("Eabs")
B1->SetLineColor(kBlue)
B1->Draw("InCmosEnSr","","same")
B1->SetLineColor(kRed)
B1->Draw("InCmosEnY","","same")
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
- Fixed DOTA source problem






