# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
stty erase '^?' #backspace and delete

# OPENSTACK
export OS_PASSWORD='$<'

# Git
export GIT_EDITOR=vim

# Remote
alias ssh='ssh -o "StrictHostKeyChecking no" -Y'

# Sys
alias dstat='dstat -cdlmnpsy'
alias sysinfo='cat /proc/cpuinfo && cat /proc/meminfo && free'

# Bash
alias df='df -h'
alias ll='ls -al -h --group-directories-first'
alias rm='sh $HOME/AMSProject/sys/shell/rmsoft.sh'
alias rmfc='/bin/rm'

# ROOT
alias root='root -l'

# AMSProject
export AMSProj=~/AMSProject
export AMSProjLibs=${AMSProj}/libs
export AMSProjProd=${AMSProj}/prod
export AMSProjSubj=${AMSProj}/subj

export AMSMKfile=${AMSProj}/env/env.mk

# Core
export AMSCore=/data3/hchou/AMSProject/core
export AMSCoreProd=${AMSCore}/prod
export AMSCoreSubj=${AMSCore}/subj

# AMS Software
alias amsenv_root5icc="source ${AMSProj}/env/setup_amsenv_root5icc.sh"
alias amsenv_root5gcc="source ${AMSProj}/env/setup_amsenv_root5gcc.sh"
