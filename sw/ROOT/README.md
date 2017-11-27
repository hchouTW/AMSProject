#AMS Offline Libs Install
```bash
amsenv_root5icc         # AFS Offline
unset AMSLIB            #
cvs co AMS              # AMS02Offline
export AMSWD=$PWD/AMSWD # AMS02Offline/AMSWD
cd AMS/install          # AMS/install
make -j8 gbatch lib     # using 8-process
cvs update              # AMS/ updata Offline
```
