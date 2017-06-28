# Theoretical Model


## Propagation Model

### Lorentz Force

### Multiple-Scattering

### Ionisation Energy Loss

### Bremsstrahlung Energy Loss


## Track Fitting Model



# Propagation Testing


## Lorentz Force

Magnetic Field: <br />
Uniform, 1 kGaus <br />
Particle Beam: <br />
(Electron Proton), Position(0, 0, 55cm), Direction(0, 0, -1), Momentum(0.3GeV 0.5GeV 1GeV 3GeV 10GeV 100GeV) <br />
Measurement Plane: <br />
Silicon Detector, Size(x y z := 100cm 100cm 300micrometre), Central Location(X Y Z := 0cm 0cm (50cm 0cm)) <br />

## Multiple-Scattering

### Direction Testing

Magnetic Field: <br />
None <br />
Particle Beam: <br />
(Electron Proton), Position(0, 0, 110cm), Direction(0, 0, -1), Momentum(0.3GeV 0.5GeV 1GeV 3GeV 10GeV 100GeV) <br />
Measurement Plane: <br />
Silicon Detector, Size(x y z := 100cm 100cm 300micrometre), Central Location(X Y Z := 0cm 0cm (106cm 0cm)) <br />
Material Box: <br />
Carbon, Size(x y z := 100cm 100cm 5cm), Central Location(X Y Z := 0cm 0cm 102.5cm) <br />

### Position & Direction Testing

Magnetic Field: <br />
None <br />
Particle Beam: <br />
(Electron Proton), Position(0, 0, 15cm), Direction(0, 0, -1), Momentum(0.3GeV 0.5GeV 1GeV 3GeV 10GeV 100GeV) <br />
Measurement Plane: <br />
Silicon Detector, Size(x y z := 100cm 100cm 300micrometre), Central Location(X Y Z := 0cm 0cm (12cm 0cm)) <br />
Material Box: <br />
Carbon, Size(x y z := 100cm 100cm 10cm), Central Location(X Y Z := 0cm 0cm 6cm) <br />

## Ionisation Energy Loss

Magnetic Field:<br /> 
None <br />
Particle Beam: <br />
(Proton Helium Carbon), Position(0, 0, 15cm), Direction(0, 0, -1), Momentum(0.3GeV 0.5GeV 1GeV 3GeV 10GeV 100GeV) <br />
Measurement Plane: <br />
Silicon Detector, Size(x y z := 100cm 100cm 300micrometre), Central Location(X Y Z := 0cm 0cm (12cm 0cm)) <br />
Material Box: <br />
Carbon, Size(x y z := 100cm 100cm 10cm), Central Location(X Y Z := 0cm 0cm 6cm) <br />

## Bremsstrahlung Energy Loss

Magnetic Field:<br /> 
None <br />
Particle Beam: <br />
Electron, Position(0, 0, 15cm), Direction(0, 0, -1), Momentum(0.3GeV 0.5GeV 1GeV 3GeV 10GeV 100GeV) <br />
Measurement Plane: <br />
Silicon Detector, Size(x y z := 100cm 100cm 300micrometre), Central Location(X Y Z := 0cm 0cm (12cm 0cm)) <br />
Material Box: <br />
Carbon, Size(x y z := 100cm 100cm 10cm), Central Location(X Y Z := 0cm 0cm 6cm) <br />


# Track Fitting Testing

Magnetic Field: <br />
(None) (Uniform, 1 kGaus) (Gaussian, 1kGaus 0cm 40cm) <br />
Particle Beam: <br />
(Electron Proton Helium Carbon), Position(0, 0, 120cm), Direction(0, 0, -1), Momentum(0.3GeV 0.5GeV 1GeV 3GeV 10GeV 100GeV 300GeV 800GeV 1000GeV 2000GeV 4000GeV 8000GeV) <br />
Measurement Plane: <br />
Silicon Detector, Size(x y z := 100cm 100cm 300micrometre), Central Location(X Y Z := 0cm 0cm (110cm 100cm 32cm 28cm 2cm -2cm =28cm -32cm -100cm -110cm)) <br />
Material Box: <br />
(None) (Carbon, Size(x y z := 100cm 100cm 8cm), Central Location(X Y Z := 0cm 0cm (60cm -60cm))) <br />


## Special Environment Testing

### Position Measurement
Resultion Models: <br />
(Single Gaussian, ) (Doublet Gaussion, ) <br />



















# Monte Carlo

## Silicon Detector (x * y * z := 150cm * 150cm * 300 micrometre)
### Central Position
#### Layer    : Position (x y z)

Layer 01 : (0 0  160)      (AMS)  
Layer 02 : (0 0  155)      (Test on 2 extenal layers)  
Layer 03 : (0 0  55)       (AMS)  
Layer 04 : (0 0  30)       (AMS)  
Layer 05 : (0 0  25)       (AMS)  
Layer 06 : (0 0   2.5)     (AMS)  
Layer 07 : (0 0  -2.5)     (AMS)  
Layer 08 : (0 0 -25)       (AMS)  
Layer 09 : (0 0 -30)       (AMS)  
Layer 10 : (0 0 -135)      (AMS)  
Layer 11 : (0 0 -140)      (Test on 2 extenal layers)  

## Aluminium (x * y * z := 150cm * 150cm * 4cm)
### Central Position
#### Layer   : Position (x y z)

Layer01 : (0 0  157.5)  
Layer02 : (0 0   70.0)  
Layer03 : (0 0   27.5)  
Layer04 : (0 0    0.0)  
Layer05 : (0 0  -27.5)  
Layer06 : (0 0  -70.0)  
Layer07 : (0 0 -137.5)  

## Carbon (x * y * z := 150cm * 150cm * 10cm)
### Central Position
#### Layer   : Position (x y z)

Layer01 : (0 0  62)  
Layer02 : (0 0 -62)  

## Magnetic Field B (unit : 1 T = 10000 Gauss) (Only have X-Component)
#### Uniform Model

Bx(unit : T) = 0.10  

#### Gaussian Model

Bx(unit : T) = 0.15 * Exp(-0.5 * (z/50cm)^2 )

## Trigger
### Hits are in (Silicon L04 or L05) and (Silicon L08 or L09)

## Particle Generator
### Initial particles position and direction

Position (x y z) = ((-150 to 150) (-150 to 150) 200) cm  
Direction (theta phi) = (0~-PI Uniform)  

#### Electron
Momentum from 0.25 GeV to 2000 GeV

#### Proton
Momentum from 0.25 GeV to 4000 GeV  

#### Helium
Momentum from 1.0 GeV to 6000 GeV  

#### Carbon
Momentum from 3.0 GeV to 12000 GeV  


# Template
## Test on different effect
### Tow Silicon Layers 
### Material (Al or C or None)   (1cm 3cm 5cm 10cm 20cm 30cm)  
### Magnetic Field (Uniform or Without)  
### (Ionisation enregy loss and/or Multiple-Scattering)






