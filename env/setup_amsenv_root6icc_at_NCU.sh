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
GCCVAR=slc6_amd64_gcc493
GCCTAG=4.9.3
GCCDIR=/cvmfs/cms.cern.ch/${GCCVAR}/external/gcc/${GCCTAG}

export PATH=${GCCDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${GCCDIR}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${GCCDIR}/lib64:${LD_LIBRARY_PATH}

#### CMS %% ROOT Environment
CMSVersion=CMSSW_7_6_3
CMSTAG=7.6.3
CMSSW=/cvmfs/cms.cern.ch/${GCCVAR}/cms/cmssw/${CMSVersion}/external/${GCCVAR}
export PATH=${CMSSW}/bin:${PATH}
export LD_LIBRARY_PATH=${CMSSW}/lib:${LD_LIBRARY_PATH}


#### CERN %% GCC Compiler
#source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6/setup.sh
#source /cvmfs/sft.cern.ch/lcg/external/gcc/6.2/x86_64-slc6/setup.sh

#### AMS %% INTEL Compiler
ICCTAG=2017
ICCDIR=/cvmfs/projects.cern.ch/intelsw/psxe/linux/x86_64/2017/compilers_and_libraries/linux
source ${ICCDIR}/bin/compilervars.sh intel64 
export PATH=${ICCDIR}/bin/intel64:${PATH}
export LD_LIBRARY_PATH=${ICCDIR}/lib/intel64:${LD_LIBRARY_PATH}

#### AMS %% ROOT Environment
AMSSW=root6-04-08-icc16
export Offline=/cvmfs/ams.cern.ch/Offline
export ROOTSYS=${Offline}/root/Linux/${AMSSW}
export PATH=${ROOTSYS}/bin:${PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

export INTEL_LICENSE_FILE=${Offline}/intel/licenses
#export INTEL_LICENSE_FILE=28518@lxlicen01.cern.ch,28518@lxlicen02.cern.ch,28518@lxlicen03.cern.ch

#### AMS %% Software Environment
AMSVersion=vdev
ROOTARCH=linuxx8664icc6.04
export AMSSRC=${Offline}/${AMSVersion}
export AMSLIB=${AMSSRC}/lib/${ROOTARCH}

export AMSWD=${AMSSRC}
export AMSDataDir=${Offline}/AMSDataDir
export AMSDataDirRW=${AMSDataDir}

echo -e "setting environment for GCC-${GCCTAG} CMSSW-${CMSTAG} ICC-${ICCTAG} AMSSW-${AMSVersion}-${ROOTARCH}"
