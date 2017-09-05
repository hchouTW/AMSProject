# Theoretical Model


## Propagation Model

### Lorentz Force

### Multiple-Scattering

### Ionisation Energy Loss

### Bremsstrahlung Energy Loss


## Track Fitting Model









# Program Testing

## Propagation Testing

### Lorentz Force

Particle Beam:
```
Amount(1M)
Electron & Proton & Deuterium & Helium4
Position(0, 0, 60cm)
Direction(0, 0, -1)

Momentum(GeV) (15 sets)
0.1  0.3  0.5  0.8  1  3  5  8  10  30  50  80  100  300  500
```

Magnetic Field:
```
Uniform, 1 kGaus in +X
```

Measurement Plane: (Silicon Detector)
```
Density(0.083 mol/cm^3 = 2.331099 g/cm^3) 
Size(x y z := 100cm 100cm 300micrometre)
Central Location(X Y Z := 0cm 0cm (58cm 50cm 0cm))
```

Variables: (Particle Information)
```
Position (cm)
Momentum (GeV)
```

#### Example
```c++
double tagz;
double cx, cy, cz;
double mx, my, mz;

PhySt part(PartType::Proton);
part.set_state(cx, cy, cz, mx, my, mz);

PropMgnt::PropToZ(tagz, part);
part.print();
```

### Multiple Scattering & Energy Loss

Particle Beam:
```
Amount(1M)
Electron & Proton & Deuterium & Helium4
Position(0, 0, 60cm)
Direction(0, 0, -1)

Momentum(GeV)  (15 sets)
0.1  0.3  0.5  0.8  1  3  5  8  10  30  50  80  100  300  500
```

Magnetic Field:
```
None
```

Material Box: (Carbon)
```
Density(0.08 mol/cm^3 = 0.9608624 g/cm^3) 
Size(x y z := 100cm 100cm 7cm)
Central Location(X Y Z := 0cm 0cm 54cm)
```

Measurement Plane: (Silicon Detector)
```
Density(0.083 mol/cm^3 = 2.331099 g/cm^3) 
Size(x y z := 100cm 100cm 300micrometre)
Central Location(X Y Z := 0cm 0cm (58cm 50cm 0cm))
```

Variables: (Particle Information)
```
Position (cm)
Momentum (GeV)
```

#### Example
```
codes
```

## Track Fitting Testing

Particle Beam:
```
Amount(1M)
Electron & Proton & Deuterium & Helium4
Position(0, 0, 110cm)
Direction(0, 0, -1)

Momentum(GeV) (21 sets)
0.1  0.3  0.5  0.8  1  3  5  8  10  30  50  80  100  300  500  800  1000  3000  5000  8000  10000
```

Magnetic Field:
```
Uniform, 1 kGaus in +X
```

Material Box: (Carbon)
```
Density(0.08 mol/cm^3 = 0.9608624 g/cm^3) 
Size(x y z := 100cm 100cm 8cm)
Central Location(X Y Z := 0cm 0cm (65cm -65cm))
```

Material Box: (Aluminium)
```
Density(0.1 mol/cm^3 = 2.698153868 g/cm^3) 
Size(x y z := 100cm 100cm 0.05cm)
Central Location(X Y Z := 0cm 0cm (26.9cm 23.1cm 1.9cm -1.9cm -23.1cm -26.9cm))
```

Measurement Plane: (Silicon Detector)
```
Density(0.083 mol/cm^3 = 2.331099 g/cm^3) 
Size(x y z := 100cm 100cm 300micrometre)
Central Location(X Y Z := 0cm 0cm (100cm 70cm 60cm 27cm 23cm 2cm -2cm -23cm -27cm -60cm -70cm -100cm))
```

Variables: (Particle Information)
```
Position (cm)
Momentum (GeV)
```

#### Example
```
codes
```

## Special Environment Testing

# Usage
## Library Requirement
```c++
#include <CCPLibs/CCPLibs.h>
#include <ROOTLibs/ROOTLibs.h>
#include <TRACKLibs/TRACKLibs.h>

// This is for special case
// #define __HAS_TESTPROP__
// #define __HAS_TESTFIT__
// #define __HAS_AMS_OFFICE_LIBS__
```
## Particle Info
```c++
PartType type = PartType::Proton;
PhySt part(type);

// coordinate and momentum setting
double coo[3];
double mom[3];
part.set_state(coo[0], coo[1], coo[2], mom[0], mom[1], mom[2]);

// show particle information
part.print();
```
## Magnetic Field
Build Magnetic Field Data Base
```c++
MagGeoBoxCreator creator;
```
Read Magnetic Field Data Base
```c++
MagGeoBoxReader reader;
```
Get Magentic Field
```c++
SVecD<3> coo;
MagFld&& mfld = MagMgnt::Get(coo);
mfld.print();
```
## Material Field
Build Material Field Data Base
```c++
MatGeoBoxCreator creator;
```
Read Material Field Data Base
```c++
MatGeoBoxReader reader;
```
Get Material Field
```c++
SVecD<3> coo;
MatFld&& mfld = MatMgnt::Get(coo);
mfld.print();
```
## Propagation
Normal Method
```c++
bool mscatOpt = true;
bool elossOpt = true;
MatArg marg(mscatOpt, elossOpt);
marg.set_mscat(0., 0.);
marg.set_eloss(0., 0.);

PhySt part;

// Propagate by step length
double tagetS;
PropMgnt::Prop(tagetS, part, marg);

// Propagate to coordinate Z
double tagetZ;
PropMgnt::PropToZ(tagetZ, part, marg);

part.print();
```
Toy Monte Carlo
```c++
bool mscatOpt = true;
bool elossOpt = true;
MatArg marg(mscatOpt, elossOpt);

PhySt part;

// Propagate by step length
double tagetS;
PropMgnt::PropWithMC(tagetS, part, marg);

// Propagate to coordinate Z
double tagetZ;
PropMgnt::PropToZWithMC(tagetZ, part, marg);

part.print();
```
## Hit Info
```c++
HitSt hit;
```
## Track Fitting
```c++
PhyTr track;
```
