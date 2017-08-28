#!/bin/bash
#shopt -s -o nounset

ulimit -c 0
ulimit -d 700000
ulimit -s 32000

unset DEBUGFLAG 
export GDBMAN=1
export GDBINFO=1
export VERBOSE=1
export LANGUAGE=US
export LANG=C
export LC_ALL=en_US

###### SLC Version
OSRelease=`grep SLC /etc/redhat-release | cut -d' ' -f6`
OSVersion=${OSRelease%.*}

#### CMS %% GCC Compiler
#GCCVAR=slc6_amd64_gcc493
#GCCTAG=4.9.3
#GCCDIR=/cvmfs/cms.cern.ch/${GCCVAR}/external/gcc/${GCCTAG}

#export PATH=${GCCDIR}/bin:${PATH}
#export LD_LIBRARY_PATH=${GCCDIR}/lib:${LD_LIBRARY_PATH}
#export LD_LIBRARY_PATH=${GCCDIR}/lib64:${LD_LIBRARY_PATH}

#### CMS %% ROOT Environment
#CMSVersion=CMSSW_7_6_3
#CMSTAG=7.6.3
#CMSSW=/cvmfs/cms.cern.ch/${GCCVAR}/cms/cmssw/${CMSVersion}/external/${GCCVAR}
#export PATH=${CMSSW}/bin:${PATH}
#export LD_LIBRARY_PATH=${CMSSW}/lib:${LD_LIBRARY_PATH}

#### CERN %% GCC Compiler
GCCTAG=4.9.1
source /cvmfs/sft.cern.ch/lcg/external/gcc/${GCCTAG}/x86_64-slc6/setup.sh

#### CERN AFS %% INTEL Compiler
ICCTAG=17
ICCDIR=/afs/cern.ch/sw/IntelSoftware/linux
source ${ICCDIR}/${ICCTAG}-all-setup.sh intel64 &> /dev/null
if [ $ICCTAG -gt 15 ]; then
    source ${ICCDIR}/x86_64/xe20${ICCTAG}/compilers_and_libraries/linux/bin/compilervars.sh intel64
else
    source ${ICCDIR}/x86_64/xe20${ICCTAG}/composerxe/bin/compilervars.sh intel64 
fi

#### CERN CVMFS %% INTEL Compiler
#ICCTAG=17
#ICCDIR=/cvmfs/projects.cern.ch/intelsw/psxe/linux
#source ${ICCDIR}/${ICCTAG}-all-setup.sh intel64 &> /dev/null
#source ${ICCDIR}/x86_64/20${ICCTAG}/compilers_and_libraries/linux/bin/compilervars.sh intel64 

#export INTEL_LICENSE_FILE=${Offline}/intel/licenses

#### AMS %% ROOT Environment
AMSSW=root-v5-34-9-icc64-slc6
export Offline=/cvmfs/ams.cern.ch/Offline
export ROOTSYS=${Offline}/root/Linux/${AMSSW}
export PATH=${ROOTSYS}/bin:${PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

#### AMS %% Software Environment
AMSVersion=vdev
ROOTARCH=linuxx8664icc5.34
export AMSSRC=${Offline}/${AMSVersion}
export AMSLIB=${AMSSRC}/lib/${ROOTARCH}

export AMSWD=${AMSSRC}
export AMSDataDir=${Offline}/AMSDataDir
export AMSDataDirRW=${AMSDataDir}

echo -e "setting environment for GCC-${GCCTAG} ICC-${ICCTAG} AMSSW-${AMSVersion}-${ROOTARCH}"
