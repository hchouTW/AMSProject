# Submit Jobs System In CERN
## Build Config File
```
[PROJECT]
PROJPATH     =  prod
PROJVERSION  =  vdev
PROJTITLE    =  TESTING
PROJBIN      =  YiProdNtuple
PROJFLIST    =  flist.cern.iss.B950.pass6
EVENTTYPE    =  ISS
FILEPEREXE   =  1
JOBREGION    =  PART
EXESATID     =  1
EXEENDID     =  10

[QUEUE]
QUEUE      =  ams1nd
STORAGE    =  EOS
EXEPERJOB  =  3
CONFIRM    =  NONE
```

## Make Job and Run
```bash
cat job.cern.testing.conf
submit RUN job.cern.testing.conf
```
