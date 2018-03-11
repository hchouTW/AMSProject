#!/bin/bash
################################################################################
##  1st Step: source CERN_TRACKSys.sh                                         ##
##  2nd Step: add OpenMP flag to makefile                                     ##
##  3rd Step: add GOOGLE_FLAGS to makefile                                    ##
##  4th Step: add GOOGLE_PARA to makefile (if needed)                         ##
##  5th Step: Complier with -std=c++11 or -std=c++14                          ##
##  6th Step: In C++ codes.                                                   ##
##            #include <TRACKSys.h>                                           ##
##            using namespace TrackSys;                                       ##
################################################################################

ExternalLibs=/eos/ams/user/h/hchou/ExternalLibs
export LIBRARY_PATH=${ExternalLibs}/LIBS/external/google/lib64:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${ExternalLibs}/LIBS/external/google/lib64:${LD_LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${ExternalLibs}/LIBS/external/eigen3/include:${CPLUS_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${ExternalLibs}/LIBS/external/google/include:${CPLUS_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${ExternalLibs}/LIBS/TRACKSys/include:${CPLUS_INCLUDE_PATH}

export GOOGLE_PARA='-DGLOG_NO_ABBREVIATED_SEVERITIES'
export GOOGLE_FLAGS='-lgflags -lglog -lceres'

export TRACKSys_MagBox=${ExternalLibs}/DB/magnetic/AMS02Mag.bin
export TRACKSys_MatBox=${ExternalLibs}/DB/material
