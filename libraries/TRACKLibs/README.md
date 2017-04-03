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
Layer03 : (0 0   57.5)  
Layer04 : (0 0   27.5)  
Layer05 : (0 0    0.0)  
Layer06 : (0 0  -27.5)  
Layer07 : (0 0  -70.0)  
Layer08 : (0 0 -137.5)  

## Carbon (x * y * z := 150cm * 150cm * 10cm)
### Central Position
#### Layer   : Position (x y z)
Layer01 : (0 0  62)  
Layer02 : (0 0 -62)  

## Magnetic Field B (unit : 1 T = 10000 Gauss) (Only have X-Component)
### Uniform Model

Bx = 0.10T  <\br>

### Gaussian Model

Bx = 0.15 * Exp(-0.5 * (z/50cm)^2 )  

## Trigger
### Hits are in (Silicon L04 or L05) and (Silicon L08 or L09)
