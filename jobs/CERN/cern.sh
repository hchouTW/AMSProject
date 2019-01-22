#!/bin/bash

# System
export CASTOR=/castor/cern.ch/user/h/hchou
export EOS_HOME=/eos/ams/user/h/hchou
export AFSWORK=/afs/cern.ch/work/h/hchou
export ubackup=/afs/cern.ch/ubackup/h/hchou

# Submit command
alias readflist='sh ${AMSProjJobs}/CERN/readflist.sh'
alias mkjob='sh ${AMSProjJobs}/CERN/mkjob.sh'
alias submit='sh ${AMSProjJobs}/CERN/submit.sh'
alias mkjob_condor='sh ${AMSProjJobs}/CERN/mkjob_condor.sh'
alias submit_condor='sh ${AMSProjJobs}/CERN/submit_condor.sh'
