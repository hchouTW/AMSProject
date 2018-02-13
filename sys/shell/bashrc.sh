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

# ROOT
alias root='root -l'

# TEX
alias texmake="sh /home/hchou/AMSProject/sw/tex/texmake.sh"
alias texclean="sh /home/hchou/AMSProject/sw/tex/texclean.sh"

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
if [[ "$HOSTNAME" == *"lxplus"* ]]; then
    alias mkjob='sh ${AMSProj}/jobs/CERN/mkjob.sh'
    alias submit='sh ${AMSProj}/jobs/CERN/submit.sh'
else
    alias mkjob='sh ${AMSProj}/jobs/NCU/mkjob.sh'
    alias submit='sh ${AMSProj}/jobs/NCU/submit.sh'
    alias voms_info='voms-proxy-info -all -file ~/.ams02'
    alias voms_init='voms-proxy-init -voms ams02.cern.ch -valid 168:00 -hours 168 -out ~/.ams02'
    alias voms_auto_init='cat ~/.globus/passwd | voms-proxy-init -rfc -voms ams02.cern.ch -cert ~/.globus/usercert.pem -key ~/.globus/userkey.pem -hours 168 -vomslife 168:0 -pwstdin'
    export X509_USER_PROXY=~/.ams02
fi

if [[ "${HOSTNAME}" == *"lxplus"* ]]; then
    export CASTOR=/castor/cern.ch/user/h/hchou
    export EOS_HOME=/eos/ams/user/h/hchou
    export AFSWORK=/afs/cern.ch/work/h/hchou
    export ubackup=/afs/cern.ch/ubackup/h/hchou
else
    # system
    export EOS_MGM_URL=root://tw-eos01.grid.sinica.edu.tw
    export EOS_HOME=/eos/ams
    export DPM_HOST=grid71.phy.ncu.edu.tw
    export DPNS_HOST=grid71.phy.ncu.edu.tw
    # user
    export VMOS_WEB=https://voms.grid.sinica.edu.tw:8443
    export DPM_HOME=/dpm/phy.ncu.edu.tw/home
    export EOS_AMS_HOST=eosams.cern.ch

    alias xrdcp_cern='sh ${AMSProj}/sys/shell/xrdcp_cern.sh'
fi

# Local Job Command
function ljsearch {
    if [ $# == 1 ]; then
        ps -U $USER -o pid -o s -o time -o command | grep ${1} | grep -v grep
    else
        ps -U $USER -o pid -o s -o time -o command
    fi
}

function ljkill {
    if [ $# == 1 ]; then
        ps -U $USER -o pid -o s -o time -o command | grep ${1} | grep -v grep | awk '{print $1}' | xargs kill
    else
        echo -e "Error: No Keyword."
    fi
}

function ljcheck {
    if [ $# == 1 ]; then
        date_beg=`date`
        while true
        do
            clear
            jobs_num=`ljsearch ${1} | wc -l`
            echo -e "DATE BEGIN ${date_beg} NOW `date`"
            echo -e "NJOBS ${jobs_num}\n"
            ljsearch ${1}
            if (( ${jobs_num} == 0 )); then
                echo -e "SUCCESS."
                break
            else
                sleep 30
            fi
        done
    else
        echo -e "Error: No Keyword."
    fi
}

function showcom {
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
