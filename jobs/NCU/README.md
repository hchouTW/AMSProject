# Submit Jobs System In NCU
## My Workflow

### Make Job and Run
```bash
submit RUN job.ncu.testing.conf
```

### Config File
```bash
cat job.ncu.testing.conf
```
```bash
[PROJECT]
PROJPATH     =  prod
PROJVERSION  =  vdev
PROJTITLE    =  TESTING
PROJBIN      =  YiProdNtuple
PROJFLIST    =  flist.dpm.iss.B950.pass6
EVENTTYPE    =  ISS
FILEPEREXE   =  1
JOBREGION    =  PART
EXESATID     =  1
EXEENDID     =  10

[QUEUE]
QUEUE      =  ams
STORAGE    =  DPM
EXEPERJOB  =  3
CONFIRM    =  NONE
```

## Certificate Authority (CA)

### CA from CERN
Go To https://ca.cern.ch/ca <br/>
-> New Grid User certificate <br/>
-> Get Grid User certificate (myCertificate.p12) <br/>

### CA from AS
https://canew.twgrid.org/request/new/request-new.php

### Copy CA to Host (hep068.phy.ncu.edu.tw)
```bash
cp myCertificate.p12 to ~/.globus/
```

### Create
create usercert.pem and userkey.pem by openssl
```bash
openssl pkcs12 -in myCertificate.p12 -clcerts -nokeys -out usercert.pem
openssl pkcs12 -in myCertificate.p12 -nocerts -out userkey.pem
chmod 400 userkey.pem
```

## VOMS (from AS)
https://voms.grid.sinica.edu.tw:8443/

### Init VOMS
In ~/.bashrc file
```bash
export X509_USER_PROXY=~/.ams02
```
Init voms
```bash
voms-proxy-init -voms ams02.cern.ch -valid 24:0 -hours 24 -out ~/.ams02
voms-proxy-info -all -file ~/.ams02
```

## PBS (Portable Batch System)
Commands
```bash
qsub -q ams -N JobName -j oe -o job.out script.sh
qstat
qselect
qdel
```

## Storage
### EOS (CERN)
```bash
eos ls /eos/ams/Data/AMS02/2014/ISS.B950/pass6
```
ISS Data in CERN <br/>
number of files = 132115 <br/>
size = 1.6P <br/>

### EOS (AS)
In ~/.bashrc file
```bash
export EOS_MGM_URL=root://tw-eos03.grid.sinica.edu.tw
export EOS_HOME=/eos/ams
```
Example:
```bash
xrd tw-eos03.grid.sinica.edu.tw ls /eos/ams/amsdatadisk/Data/2014/ISS.B950/pass6
```

### DPM (Disk Pool Manager)
In ~/.bashrc file
```bash
export DPM_HOST=grid71.phy.ncu.edu.tw
export DPNS_HOST=grid71.phy.ncu.edu.tw
```
Example:
```bash
xrdfs grid71.phy.ncu.edu.tw
dpns-ls /dpm/phy.ncu.edu.tw/home/ams02
```
DPM Commands
```bash
dpm-qryconf
dpm-listspaces
dpns-ls
rfcp
rfrm
```

### Tranfer Data from AS to NCU (Xrd tools)
Example:
```bash
xrdcp root://tw-eos03.grid.sinica.edu.tw//eos/ams/amsdatadisk/Data/2014/ISS.B950/pass6/1411871829.00000001.root root://grid71.phy.ncu.edu.tw:1094//ams02/ams02datadisk/1411871829.00000001.root
```
