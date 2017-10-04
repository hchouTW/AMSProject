## Certificate Authority (CA)

### CA from CERN
Go To https://ca.cern.ch/ca
New Grid User certificate
Get Grid User certificate (myCertificate.p12)

### CA from AS
https://canew.twgrid.org/request/new/request-new.php

### Copy CA to Host (hep068.phy.ncu.edu.tw)
```
mv myCertificate.p12 to ~/.globus/
```

### Create
create usercert.pem and userkey.pem by openssl
```
openssl pkcs12 -in xxxx.p12 -clcerts -nokeys -out usercert.pem
openssl pkcs12 -in xxxx.p12 -nocerts -out userkey.pem
```

## VOMS (from AS)
https://voms.grid.sinica.edu.tw:8443/

### Init VOMS
```
cd ~/.globus
voms-proxy-init --voms ams02.cern.ch --out ams02
voms-proxy-init --voms ams02.cern.ch
voms-proxy-info --all
```

## Storage
### EOS (CERN)
/eos/ams/Data/AMS02/2014/ISS.B950/pass6
number of files = 132115
size = 1.6P 

### EOS (AS)
xrd tw-eos02.grid.sinica.edu.tw ls /eos/ams/amsdatadisk/Data/2014/ISS.B950/pass6

### DPM (NCU)
xrdfs root://grid71.phy.ncu.edu.tw:1094/
dpns-ls /dpm/phy.ncu.edu.tw/home/ams02/ams02datadisk
dpm-qryconf
dpm-listspaces

### TEST (Copy files from AS to NCU)
xrdcp root://tw-eos03.grid.sinica.edu.tw//eos/ams/amsdatadisk/Data/2014/ISS.B950/pass6/1411871829.00000001.root root://grid71.phy.ncu.edu.tw:1094//ams02/ams02datadisk/1411871829.00000001.root
