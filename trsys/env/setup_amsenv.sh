#!/bin/bash
#shopt -s -o nounset

export LANGUAGE=US
export LANG=C
export LC_ALL=en_US
export GDBMAN=1
export GDBINFO=1
export VERBOSE=1

#### CERN CVMFS %% GCC Compiler
GCCTag=6.4.0
#GCCTag=8.2.0
GCCDir=/cvmfs/sft.cern.ch/lcg/contrib/gcc/${GCCTag}/x86_64-centos7
source ${GCCDir}/setup.sh

#### CERN CVMFS && CLANG Compiler
#CLANGTag=7.0.0
#CLANGDir=/cvmfs/sft.cern.ch/lcg/contrib/clang/${CLANGTag}/x86_64-centos7
#source ${CLANGDir}/setup.sh

#### CERN CVMFS %% INTEL Compiler
ICCTag=18
ICCDir=/cvmfs/projects.cern.ch/intelsw/psxe/linux
ICCLux=${ICCDir}/x86_64/20${ICCTag}/compilers_and_libraries/linux
source ${ICCDir}/${ICCTag}-all-setup.sh intel64 &> /dev/null
source ${ICCLux}/bin/compilervars.sh intel64

export INTELBIN=${ICCLux}/bin/intel64
export IFORTBIN=${ICCLux}/bin/intel64
export INTELLIB=${ICCLux}/lib/intel64

#### CERN CVMFS %% ROOT Environment
#ROOTVersion=6.16.00
#ROOTDir=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/${ROOTVersion}/x86_64-centos7-gcc48-opt
#
#source ${ROOTDir}/bin/thisroot.sh
#export LIBRARY_PATH=${ROOTDir}/lib:${LIBRARY_PATH}
#export LD_LIBRARY_PATH=${ROOTDir}/lib:${LD_LIBRARY_PATH}
#export CPLUS_INCLUDE_PATH=${ROOTDir}/include:${CPLUS_INCLUDE_PATH}

#### CERN CVMFS %% AMS Environment
export Offline=/cvmfs/ams.cern.ch/Offline

#### CERN CVMFS %% AMS ROOT Environment
ROOTVersion=v5-34-9-icc64-slc6
ROOTDir=${Offline}/root/Linux/root-${ROOTVersion}

export PATH=${ROOTDir}/bin:${PATH}
export LIBRARY_PATH=${ROOTDir}/lib:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${ROOTDir}/lib:${LD_LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${ROOTDir}/include:${CPLUS_INCLUDE_PATH}

#### EOS %% Google Libs
GoogleDir=/eos/home-h/hchou/ExternalLibs/LIBS/external/google
export LIBRARY_PATH=${GoogleDir}/lib64:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${GoogleDir}/lib64:${LD_LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${GoogleDir}/include:${CPLUS_INCLUDE_PATH}

#### EOS %% Eigen Libs
EigenDir=/eos/home-h/hchou/ExternalLibs/LIBS/external/eigen3
export CPLUS_INCLUDE_PATH=${EigenDir}/include:${CPLUS_INCLUDE_PATH}

#### Self defination libs
export LIBRARY_PATH=${PWD}/lib:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${PWD}/lib:${CPLUS_INCLUDE_PATH}

export CPLUS_INCLUDE_PATH=~/AMSProject/libs/CPPLibs:${CPLUS_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=~/AMSProject/libs/ROOTLibs:${CPLUS_INCLUDE_PATH}

#export CPLUS_INCLUDE_PATH=~/AMSProject/libs/TRACKSys/include:${CPLUS_INCLUDE_PATH}
#export CPLUS_INCLUDE_PATH=~/TrackSys/inc:${CPLUS_INCLUDE_PATH}

TrSysDir=/afs/cern.ch/user/h/hchou/TrackSys
export CPLUS_INCLUDE_PATH=${TrSysDir}/inc:${CPLUS_INCLUDE_PATH}
export LIBRARY_PATH=${TrSysDir}/lib:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${TrSysDir}/lib:${LD_LIBRARY_PATH}

#### COMPILER
export COMPILER=ICC

echo -e "setting environment for GCC-${GCCTag} CLANG-${CLANGTag} ICC-${ICCTag} ROOT-${ROOTVersion}"
