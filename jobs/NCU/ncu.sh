#!/bin/bash

# VOMS
export X509_USER_PROXY=~/.ams02
export VMOS_WEB=https://voms.grid.sinica.edu.tw:8443
alias voms_info='voms-proxy-info -all -file ~/.ams02'
alias voms_init='voms-proxy-init -voms ams02.cern.ch -valid 168:00 -hours 168 -out ~/.ams02'
alias voms_auto_init='cat ~/.globus/passwd | voms-proxy-init -rfc -voms ams02.cern.ch -cert ~/.globus/usercert.pem -key ~/.globus/userkey.pem -hours 168 -vomslife 168:0 -pwstdin'

# NCU DPM
export DPM_HOST=grid71.phy.ncu.edu.tw
export DPNS_HOST=grid71.phy.ncu.edu.tw
export DPM_HOME=/dpm/phy.ncu.edu.tw/home

# AS EOS
export EOS_MGM_URL=root://tw-eos01.grid.sinica.edu.tw
export EOS_HOME=/eos/ams

# Libs
source ${AMSProjLibs}/external/castor/bin/castor.sh

# Submit Script
alias mkjob='sh ${AMSProjJobs}/NCU/mkjob.sh'
alias submit='sh ${AMSProjJobs}/NCU/submit.sh'
alias readflist='sh ${AMSProjJobs}/readflist.sh'
alias xrdcp_cern2ncu='sh ${AMSProjJobs}/NCU/xrdcp_cern2ncu.sh'
