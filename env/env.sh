#!bin/bash

# User specific aliases and functions
stty erase '^?' #backspace and delete

# OPENSTACK
export OS_PASSWORD='$<'

# Git
export GIT_EDITOR=vim

# ROOT
alias root='root -l'

# Variables
export AMSProj="/afs/cern.ch/user/h/hchou/private/AMSProject"
export AMSProjLibs="${AMSProj}/libs"
export AMSProjProd="${AMSProj}/prod"
export AMSProjSubj="${AMSProj}/subj"

export AMSMKfile="${AMSProj}/env/env.mk"
