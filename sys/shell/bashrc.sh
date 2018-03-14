#!/bin/bash

# Backspace and Delete
stty erase '^?'

# OPENSTACK
export OS_PASSWORD='$<'

# Git
export GIT_EDITOR=vim

# Remote
alias ssh='ssh -o "StrictHostKeyChecking no" -Y'
alias rsync='rsync --ignore-existing --progress --human-readable'
alias rsync_cern="~/AMSProject/sys/shell/rsync_cern.sh"

# Sys
alias dstat='dstat -cdlmnpsy'
alias sysinfo='cat /proc/cpuinfo && cat /proc/meminfo && free -h'

# Bash
alias df="df -h"
alias ls="ls --color=auto"
alias ll="ls --color=auto -al -h --group-directories-first"
alias rm="sh ~/AMSProject/sys/shell/rmsoft.sh"
alias rmfc="/bin/rm"

# AMSProject
if [[ "${AMSProj}" == "" ]]; then export AMSProj=~/AMSProject; fi
export AMSProjLibs=${AMSProj}/libs
export AMSProjProd=${AMSProj}/prod
export AMSProjSubj=${AMSProj}/subj
export AMSProjJobs=${AMSProj}/jobs

# AMSCore
if [[ "${AMSCore}" == "" ]]; then export AMSCore=~/AMSCore; fi
export AMSCoreProd=${AMSCore}/prod
export AMSCoreSubj=${AMSCore}/subj

# AMSJobs
if [[ "${AMSJobs}" == "" ]]; then export AMSJobs=~/AMSJobs; fi

# AMSData
if [[ "${AMSData}" == "" ]]; then export AMSData=~/AMSData; fi

# TEX
source ${AMSProj}/sw/tex/tex.sh

# Eigen3 Libs
source ${AMSProjLibs}/external/eigen3/bin/eigen3.sh 

# Google Libs
source ${AMSProjLibs}/external/google/bin/google.sh

# AMS Software
source ${AMSProjLibs}/AMSLibs/bin/amslibs.sh
alias amsenv_root5icc="source ${AMSProj}/sw/ROOT/setup_amsenv_root5icc.sh"
alias amsenv_root5gcc="source ${AMSProj}/sw/ROOT/setup_amsenv_root5gcc.sh"
alias amsenv=amsenv_root5gcc

# TRACKSys Software
source ${AMSProjLibs}/TRACKSys/bin/TRACKSys.sh

# Makefile
export AMSMKfile=${AMSProj}/sw/ROOT/makefile.env

# ROOT
alias root='root -l'

# Castor Libs
source ${AMSProjLibs}/external/castor/bin/castor.sh

# Local Job Command
source ${AMSProjJobs}/LOCAL/local.sh

## Job Config
source ${AMSProj}/sys/shell/ini_parser.sh
if [[ "$HOSTNAME" == *"cern"* ]]; then
    #alias amsenv=amsenv_root5icc
    source ${AMSProjJobs}/CERN/cern.sh
    alias scphep='scp -r -oProxyCommand="ssh -W %h:%p cms01.phy.ncu.edu.tw" hep068.phy.ncu.edu.tw:'
else
    #alias amsenv=amsenv_root5gcc
    source ${AMSProjJobs}/NCU/ncu.sh
fi

function showcmd {
    # ANSI escape codes
    GREEN='\033[0;32m'
    BLUE='\033[0;34m'
    NC='\033[0m'
    echo -e ""
    echo -e "${GREEN}AMS Software Command:${NC}"
    echo -e "  ${BLUE}amsenv${NC}"
    echo -e "  ${BLUE}root${NC}"
    echo -e "${GREEN}VOM Command: [NCU]${NC}"
    echo -e "  ${BLUE}voms_info${NC}"
    echo -e "  ${BLUE}voms_auto_init${NC}"
    echo -e "${GREEN}EOS Command:${NC}"
    echo -e "  ${BLUE}eos${NC}"
    echo -e "${GREEN}DPM Command:${NC}"
    echo -e "${GREEN}[NCU_URL  root://grid71.phy.ncu.edu.tw]${NC}"
    echo -e "${GREEN}[NCU_HOST /dpm/phy.ncu.edu.tw/home]${NC}"
    echo -e "  ${BLUE}dpm-listspaces${NC}"
    echo -e "  ${BLUE}dpm-qryconf${NC}"
    echo -e "  ${BLUE}dpns-rm${NC}"
    echo -e "  ${BLUE}dpns-ls${NC}"
    echo -e "  ${BLUE}dpns-mkdir${NC}"
    echo -e "${GREEN}Local Jobs Command:${NC}"
    echo -e "  ${BLUE}ljsearch${NC}"
    echo -e "  ${BLUE}ljkill${NC}"
    echo -e "  ${BLUE}ljcheck${NC}"
    echo -e "${GREEN}PBS Jobs Command: [NCU]${NC}"
    echo -e "  ${BLUE}qsub qstat qselect qdel${NC}"
    echo -e "${GREEN}LSF Jobs Command: [CERN]${NC}"
    echo -e "  ${BLUE}bsub bjobs bpeek${NC}"
    echo -e "${GREEN}SelfJobs Command:${NC}"
    echo -e "  ${BLUE}mkjob${NC}"
    echo -e "  ${BLUE}submit${NC}"
    echo -e "${GREEN}Copy Command: [From CERN to NCU]${NC}"
    echo -e "  ${BLUE}xrdcp_cern${NC} in out [key]"
    echo -e "  ${BLUE}rsync_cern${NC} pw in out"
    echo -e ""
}
