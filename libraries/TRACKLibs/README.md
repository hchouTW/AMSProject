# Monte Carlo

## Silicon Detector (x * y * z := 150cm * 150cm * 300 micrometre)
### Central Position
#### Layer    : Position (x y z)
Layer 01 : (0 0  160)      (AMS)  <\br>
Layer 02 : (0 0  155)      (Test on 2 extenal layers)  <\br>
Layer 03 : (0 0  55)       (AMS)  <\br>
Layer 04 : (0 0  30)       (AMS)  <\br>
Layer 05 : (0 0  25)       (AMS)  <\br>
Layer 06 : (0 0   2.5)     (AMS)  <\br>
Layer 07 : (0 0  -2.5)     (AMS)  <\br>
Layer 08 : (0 0 -25)       (AMS)  <\br>
Layer 09 : (0 0 -30)       (AMS)  <\br>
Layer 10 : (0 0 -135)      (AMS)  <\br>
Layer 11 : (0 0 -140)      (Test on 2 extenal layers)  <\br>

## Aluminium (x * y * z := 150cm * 150cm * 4cm)
### Central Position
#### Layer   : Position (x y z)

Layer01 : (0 0  157.5)  <\br>
Layer02 : (0 0   70.0)  <\br>
Layer03 : (0 0   57.5)  <\br>
Layer04 : (0 0   27.5)  <\br>
Layer05 : (0 0    0.0)  <\br>
Layer06 : (0 0  -27.5)  <\br>
Layer07 : (0 0  -70.0)  <\br>
Layer08 : (0 0 -137.5)  <\br>

## Carbon (x * y * z := 150cm * 150cm * 10cm)
### Central Position
#### Layer   : Position (x y z)
Layer01 : (0 0  62)  <\br>
Layer02 : (0 0 -62)  <\br>

## Magnetic Field B (unit : 1 T = 10000 Gauss) (Only have X-Component)
### Uniform Model

Bx = 0.10T  <\br>

### Gaussian Model

Bx = 0.15 * Exp(-0.5 * (z/50cm)^2 )  <\br>

## Trigger
### Hits are in (Silicon L04 or L05) and (Silicon L08 or L09)
