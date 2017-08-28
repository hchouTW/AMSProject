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

# Trash
mkdir -p $HOME/.trash
export trash=$HOME/.trash

# Bash
alias df='df -h'
alias ll='ls -al -h --group-directories-first'
alias rm='sh $HOME/AMSProject/sw/shell/rmsoft.sh'
alias rmfc='/bin/rm'

# ROOT
alias root='root -l'

# AMSProject
export AMSProj="$HOME/AMSProject"
export AMSProjLibs="${AMSProj}/libs"
export AMSProjProd="${AMSProj}/prod"
export AMSProjSubj="${AMSProj}/subj"

export AMSMKfile="${AMSProj}/env/env.mk"
