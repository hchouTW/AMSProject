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
Electron & Proton & Helium
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
Position (mm)
Momentum (MeV)
```

#### Example
```
codes
```

### Multiple Scattering & Energy Loss

Particle Beam:
```
Amount(1M)
Electron & Proton & Helium
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
Position (mm)
Momentum (MeV)
```

#### Example
```
codes
```

## Track Fitting Testing

Particle Beam:
```
Amount(1M)
Electron & Proton & Helium
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
Position (mm)
Momentum (MeV)
```

#### Example
```
codes
```

## Special Environment Testing
