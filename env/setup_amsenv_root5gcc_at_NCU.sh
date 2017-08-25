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

#### SLC Version
OSRelease=`grep SLC /etc/redhat-release | cut -d' ' -f6`
OSVersion=${OSRelease%.*}

#### CMS %% GCC Compiler
Compiler=slc6_amd64_gcc491
CompilerTag=4.9.1-cms

export PATH=/cvmfs/cms.cern.ch/${Compiler}/external/gcc/${CompilerTag}/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/${Compiler}/external/gcc/${CompilerTag}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/${Compiler}/external/gcc/${CompilerTag}/lib64:${LD_LIBRARY_PATH}

#### CMS %% ROOT Environment
CMSSW=CMSSW_7_3_6
export PATH=/cvmfs/cms.cern.ch/${Compiler}/cms/cmssw/${CMSSW}/external/${Compiler}/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/${Compiler}/cms/cmssw/${CMSSW}/external/${Compiler}/lib:${LD_LIBRARY_PATH}

#### AMS %% ROOT Environment
AMSSW=root-v5-34-9-gcc64-slc6
export Offline=/cvmfs/ams.cern.ch/Offline
export ROOTSYS=${Offline}/root/Linux/${AMSSW}
export PATH=${ROOTSYS}/bin:${PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

#### AMS %% Software Environment
AMSVersion=vdev
AMSVersion_TXT=${AMSVersion}
if [[ "${AMSVersion_TXT}" != "vdev" ]];then
	AMSVersion_TXT="B${AMSVersion}_patches"
fi

ROOTARCH=linuxx8664gcc5.34
export AMSSRC=${Offline}/${AMSVersion_TXT}
export AMSLIB=${AMSSRC}/lib/${ROOTARCH}

export AMSWD=${AMSSRC}
export AMSDataDir=${Offline}/AMSDataDir
export AMSDataDirRW=${AMSDataDir}

echo -e "setting environment for CMS-${Compiler} ${CMSSW} AMS-${AMSSW} AMSSW-${AMSVersion_TXT}-${ROOTARCH}"
