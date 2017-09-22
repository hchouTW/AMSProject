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
0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1  2  3  4  5  8  10  30  50  80  100  300  500
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

Momentum(GeV)  (22 sets)
0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1  2  3  4  5  8  10  30  50  80  100  300  500
```

Magnetic Field:
```
None
```

Material Box: (Carbon)
```
Density(0.01 mol/cm^3 = 0.1201078 g/cm^3) 
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

Momentum(GeV) (28 sets)
0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1  2  3  4  5  8  10  30  50  80  100  300  500  800  1000  3000  5000  8000  10000
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
int xn, yn, zn;
double xmin, xmax;
double ymin, ymax;
double zmin, zmax;
std::string output_path;
MagGeoBoxCreator creator(xn, xmin, xmax, yn, ymin, ymax, zn, zmin, zmax, output_path);

double mag[3];
creator.save_and_close(mag[0], mag[1], mag[2]);
```
Read Magnetic Field Data Base
```c++
std::string input_path;
MagGeoBoxReader reader(input_path);

if (reader.exist()) {
    SVecD<3> coo;
    MagFld&& mfld = reader.get(coo);
    mfld.print();
}
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
int xn, yn, zn;
double xmin, xmax;
double ymin, ymax;
double zmin, zmax;
std::string output_path;
MatGeoBoxCreator creator(xn, xmin, xmax, yn, ymin, ymax, zn, zmin, zmax, output_path);

bool   elm[MatProperty::NUM_ELM];
double den[MatProperty::NUM_ELM];
creator.save_and_close(elm, den);
```
Read Material Field Data Base
```c++
std::string input_path;
MatGeoBoxReader reader(input_path);

// get material field by point
if (reader.exist()) {
    SVecD<3> coo;
    MatFld&& mfld = reader.get(coo);
    mfld.print();
}

// get material field by path
if (reader.exist()) {
    bool is_std_step;
    SVecD<3> vcoo, wcoo;
    MatFld&& mfld = reader.get(coo, wcoo, is_std_step);
    mfld.print();
}
```
Get Material Field
```c++
SVecD<3> coo;
MatFld&& mfld = MatMgnt::Get(coo);
mfld.print();
```
Get Material Physical Field
```c++
PhySt part;
double step_len;
MatPhyFld&& mpfld = MatPhy::Get(step_len, part);
// TODO : add function
// mpfld.print();
```
## Propagation
Normal Method
```c++
// Particle
bool mscatOpt = true;
bool elossOpt = true;
PhyArg::SetOpt(mscatOpt, elossOpt);

double coo[3];
double dir[3];
double rig;
PhySt part(PartType::Proton);
part.set_state_with_cos(coo[0], coo[1], coo[2], dir[0], dir[1], dir[2]);
part.set_rig(rig);

// Material Field
MatFld mfld;

// Jacobian for fitting
PhyJb* phyJb = nullptr;

// Propagate by step length
double tagetS;
PropMgnt::Prop(tagetS, part, &mfld, phyJb);

part.print();

// pdf parameters

// angle (multi-gaussian)
part.vst().mscatu();

// position (multi-gaussian)
part.vst().mscatcu();
part.vst().mscatcl();

// ion-energy loss (modified-Moyal eq.)
part.vst().eloss_ion_kpa();
part.vst().eloss_ion_mos();

// Bremsstrahlung-energy loss (gamma) (Now, Turn off)
part.vst().eloss_brm_men();

// 1st way
part.arg().set_mscat();
part.arg().set_eloss();

part.symbk(); // symmetry break

// 2nd way
part.symbk(true); // symmetry break

part.print();
```
Toy Monte Carlo
```c++
// Particle
bool mscatOpt = true;
bool elossOpt = true;

double coo[3];
double dir[3];
double rig;
PhySt part(PartType::Proton, mscatOpt, elossOpt);
part.set_state_with_cos(coo[0], coo[1], coo[2], dir[0], dir[1], dir[2]);
part.set_rig(rig);

// Propagate to coordinate Z
double tagetZ;
PropMgnt::PropToZWithMC(tagetZ, part);

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
