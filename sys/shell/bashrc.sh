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
alias ll="ls -al -h --group-directories-first"
alias rm="sh ~/AMSProject/sys/shell/rmsoft.sh"
alias rmfc="/bin/rm"

# ROOT
alias root='root -l'

# AMSProject
#export AMSProj=~/AMSProject # define in ~/.bashrc
export AMSProjLibs=${AMSProj}/libs
export AMSProjProd=${AMSProj}/prod
export AMSProjSubj=${AMSProj}/subj

# AMSCore
#export AMSCore=~/AMSCore # define in ~/.bashrc
export AMSCoreProd=${AMSCore}/prod
export AMSCoreSubj=${AMSCore}/subj

# AMS Software
export AMSMKfile=${AMSProj}/sw/ROOT/makefile.env
alias amsenv_root5icc="source ${AMSProj}/sw/ROOT/setup_amsenv_root5icc.sh"
alias amsenv_root5gcc="source ${AMSProj}/sw/ROOT/setup_amsenv_root5gcc.sh"


# LXPLUS
if [[ $HOSTNAME == *"lxplus"* ]]; then
    export CASTOR=/castor/cern.ch/user/h/hchou
    export EOS=/eos/ams/user/h/hchou
    export AFSWORK=/afs/cern.ch/work/h/hchou
    export ubackup=/afs/cern.ch/ubackup/h/hchou
fi
