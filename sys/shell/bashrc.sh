#!/bin/bash

# Backspace and Delete
stty erase '^?'

# OPENSTACK
export OS_PASSWORD='$<'

# Git
export GIT_EDITOR=vim

# Remote
alias ssh='ssh -o "StrictHostKeyChecking no" -Y'

# Sys
alias dstat='dstat -cdlmnpsy'
alias sysinfo='cat /proc/cpuinfo && cat /proc/meminfo && free -h'

# Bash
alias df="df -h"
alias ls="ls --color=auto"
alias ll="ls --color=auto -al -h --group-directories-first"
alias rm="sh ~/AMSProject/sys/shell/rmsoft.sh"
alias rmfc="/bin/rm"

# ROOT
alias root='root -l'

# AMSProject
#export AMSProj=~/AMSProject # define in ~/.bashrc
export AMSProjLibs=${AMSProj}/libs
export AMSProjProd=${AMSProj}/prod
export AMSProjSubj=${AMSProj}/subj
export AMSProjJobs=${AMSProj}/jobs

# AMSCore
#export AMSCore=~/AMSCore # define in ~/.bashrc
export AMSCoreProd=${AMSCore}/prod
export AMSCoreSubj=${AMSCore}/subj

# AMS Software
export AMSMKfile=${AMSProj}/sw/ROOT/makefile.env
alias amsenv_root5icc="source ${AMSProj}/sw/ROOT/setup_amsenv_root5icc.sh"
alias amsenv_root5gcc="source ${AMSProj}/sw/ROOT/setup_amsenv_root5gcc.sh"
alias amsenv=amsenv_root5gcc

# Read ROOT file list
alias readflist="sh ${AMSProj}/sys/shell/readflist.sh"

# Job Config
source ${AMSProj}/sys/shell/ini_parser.sh
if [[ $HOSTNAME == *"lxplus"* ]]; then
    alias mkjob='sh ${AMSProj}/jobs/CERN/mkjob.sh'
    alias submit='sh ${AMSProj}/jobs/CERN/submit.sh'
else
    alias mkjob='sh ${AMSProj}/jobs/NCU/mkjob.sh'
    alias submit='sh ${AMSProj}/jobs/NCU/submit.sh'
    export X509_USER_PROXY=/ams_home/hchou/ams02
    alias voms_info='voms-proxy-info --all'
    alias voms_init='voms-proxy-init --voms ams02.cern.ch --hours 168 --out ~/ams02'
fi

if [[ $HOSTNAME == *"lxplus"* ]]; then
    export CASTOR=/castor/cern.ch/user/h/hchou
    export EOS=/eos/ams/user/h/hchou
    export AFSWORK=/afs/cern.ch/work/h/hchou
    export ubackup=/afs/cern.ch/ubackup/h/hchou
else
    export VMOS_WEB=https://voms.grid.sinica.edu.tw:8443
    export EOS_MGM_URL=root://tw-eos03.grid.sinica.edu.tw
    export EOS_HOME=/eos/ams
    export DPM_HOST=grid71.phy.ncu.edu.tw
    export DPNS_HOST=grid71.phy.ncu.edu.tw
    export DPM_HOME=/dpm/phy.ncu.edu.tw/home/ams02
fi
