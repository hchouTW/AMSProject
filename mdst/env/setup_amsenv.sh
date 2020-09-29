#!/bin/bash
#shopt -s -o nounset

export LANGUAGE=US
export LANG=C
export LC_ALL=en_US
export GDBMAN=1
export GDBINFO=1
export VERBOSE=1

#### CERN CVMFS %% GCC Compiler
GCCTag=4.9.3
#GCCTag=5.3.0
#GCCTag=6.4.0
#GCCTag=8.2.0
GCCDir=/cvmfs/sft.cern.ch/lcg/contrib/gcc/${GCCTag}/x86_64-slc6
#GCCDir=/cvmfs/sft.cern.ch/lcg/contrib/gcc/${GCCTag}/x86_64-centos7
source ${GCCDir}/setup.sh

#### CERN CVMFS && CLANG Compiler
#CLANGTag=7.0.0
#CLANGDir=/cvmfs/sft.cern.ch/lcg/contrib/clang/${CLANGTag}/x86_64-centos7
#source ${CLANGDir}/setup.sh

#### CERN CVMFS %% INTEL Compiler
ICCTag=19
ICCDir=/cvmfs/projects.cern.ch/intelsw/psxe/linux
ICCLux=${ICCDir}/x86_64/20${ICCTag}/compilers_and_libraries/linux
source ${ICCDir}/${ICCTag}-all-setup.sh intel64 &> /dev/null
source ${ICCLux}/bin/compilervars.sh intel64

export INTELBIN=${ICCLux}/bin/intel64
export IFORTBIN=${ICCLux}/bin/intel64
export INTELLIB=${ICCLux}/lib/intel64

#### CERN CVMFS %% ROOT Environment
#ROOTDir=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-centos7-gcc48-opt
#
#source ${ROOTDir}/bin/thisroot.sh
#export LIBRARY_PATH=${ROOTDir}/lib:${LIBRARY_PATH}
#export LD_LIBRARY_PATH=${ROOTDir}/lib:${LD_LIBRARY_PATH}
#export CPLUS_INCLUDE_PATH=${ROOTDir}/include:${CPLUS_INCLUDE_PATH}

##### CERN CVMFS %% AMS Environment
export Offline=/cvmfs/ams.cern.ch/Offline
export AMSDataDir=${Offline}/AMSDataDir

#### CERN CVMFS %% AMS ROOT Environment
ROOTVersion=v5-34-9-icc64
ROOTDir=${Offline}/root/Linux/root-${ROOTVersion}
export ROOTSYS=${ROOTDir}

export PATH=${ROOTDir}/bin:${PATH}
export LIBRARY_PATH=${ROOTDir}/lib:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${ROOTDir}/lib:${LD_LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${ROOTDir}/include:${CPLUS_INCLUDE_PATH}

##### CERN CVMFS %% AMS Libs Environment
AMSArch=linuxx8664icc5.34

AMSVersion=vdev
export AMSOffline=${Offline}/${AMSVersion}

#AMSVersion=V.20190813
#export AMSOffline=/eos/home-h/hchou/AMSOffLibs/${AMSVersion}
#export AMSOffline=/afs/cern.ch/work/h/hchou/public/AMSOffLibs/${AMSVersion}

export LIBRARY_PATH=${AMSOffline}/lib/${AMSArch}:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${AMSOffline}/lib/${AMSArch}:${LD_LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${AMSOffline}/include:${CPLUS_INCLUDE_PATH}

##### CERN CVMFS %% AMS CERNLib Environment
export CERNLIB=${Offline}/AMSsoft/linux_slc6_icc64/2005/lib
export LIBRARY_PATH=${CERNLIB}:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${CERNLIB}:${LD_LIBRARY_PATH}

#### AFS %% AMS Geant4 Environment
G4Dir=/afs/cern.ch/ams/Offline/vdev/g4103

export CVSEDITOR=vim
export CVSROOT=${AFSOffline}/CVS

export AMSICC=1
export AMSP=1
export G4AMS=1
export PGTRACK=1
export ECALBDT=1

#### AFS %% AMS Event Display
export AMSGeoDir=${AMSOffline}/display/
alias amsedc="${AMSOffline}/exe/${AMSArch}/amsedcPG"
alias amsedd="${AMSOffline}/exe/${AMSArch}/amseddPG"

#### EOS %% Google Libs
#GoogleDir=/eos/user/h/hchou/ExternalLibs/LIBS/external/google
GoogleDir=/eos/ams/user/h/hchou/ExternalLibs/LIBS/external/google
#GoogleDir=/afs/cern.ch/work/h/hchou/public/ExternalLibs/LIBS/external/google
export LIBRARY_PATH=${GoogleDir}/lib64:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${GoogleDir}/lib64:${LD_LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${GoogleDir}/include:${CPLUS_INCLUDE_PATH}

#### AFS %% Track Sys
export TrackSysDir=/afs/cern.ch/user/h/hchou/TrackSys
export CPLUS_INCLUDE_PATH=${TrackSysDir}/inc:${CPLUS_INCLUDE_PATH}

#### Self defination libs
export LIBRARY_PATH=${PWD}/lib:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH}

#### COMPILER
export COMPILER=ICC

echo -e "setting environment for GCC-${GCCTag} CLANG-${CLANGTag} ICC-${ICCTag} ROOT-${ROOTVersion} AMSSW-${AMSVersion}-${AMSArch}"
echo -e ""
which g++
which icpc
which root
echo -e ""
