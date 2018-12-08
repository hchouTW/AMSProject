#!/bin/bash
export CPLUS_INCLUDE_PATH=${AMSProjLibs}/AMSLibs/EcalHadron:${CPLUS_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${AMSProjLibs}/AMSLibs/JFbeta:${CPLUS_INCLUDE_PATH}
source ${AMSProjLibs}/AMSLibs/TRDVertex/bin/trdvtx.sh
