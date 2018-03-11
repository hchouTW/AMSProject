# Theoretical Model


## Propagation Model

### Lorentz Force

### Multiple-Scattering

### Ionisation Energy Loss

### Bremsstrahlung Energy Loss

## Track Fitting Model

# Usage
## Library Requirement
```c++
#include <TRACKSys.h>
using namespace TrackSys;
```
## Particle Info
```c++
PhySt part(PartType::Proton);

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
long long n[3];
double    min[3];
double    max[3];
std::string opath;
MagGeoBoxCreator creator(n, min, max, opath);

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
long long n[3];
double    min[3];
double    max[3];
std::string opath;
MatGeoBoxCreator creator(n, min, max, opath);

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
// position (multi-gaussian)
// ion-energy loss (modified-Moyal eq.)
// Bremsstrahlung-energy loss (gamma) (Now, Turn off)
part.arg().nrl();
part.arg().ela();

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
PhyFit trFit;
```
