#COMPILER = icc
COMPILER = gcc
DEFS = -D_PGTRACK_ -D__ROOTSHAREDLIBRARY__
PARA = -D_X8664__ -D_AMSPARALLEL__
#IOPT = -axsse4.2,ssse3 # Only for Intel Core
IOPT =
CXX  = $(COMPILER) $(IOPT) -std=c++11 $(PARA) -O2 -D_DEBUG -Wall
FLGS = -I$(ROOTSYS)/include -isystem$(ROOTSYS)/include -I$(AMSSRC)/include -isystem$(AMSSRC)/include -I$(AMSProjLibs) -isystem$(AMSProjLibs)

OPTS = -g
STAT = -static-intel
LD   = $(shell ${ROOTSYS}/bin/root-config --ld) $(IOPT) $(PARA)
LDS  = $(shell ${ROOTSYS}/bin/root-config --ld) $(IOPT) $(PARA) $(STAT)

LIBA = $(AMSLIB)/libntuple_slc6_PG.a
LIBD = $(AMSLIB)/libntuple_slc6_PG_dynamic.so

CERNDIR = /afs/cern.ch/exp/ams/Offline/CERN/2005
LIBP  = $(shell $(ROOTSYS)/bin/root-config --libs)
LIBP += -lifcore -L$(CERNDIR)/lib -lmathlib
LIBP += -lMinuit -lMinuit2 -lRFIO -lTMVA -lXMLIO -lMLP -lTreePlayer -lMathCore
LIBP += -llzma -Llib -lcrypto -L$(XRDLIB)/lib64 -lXrdClient -lXrdUtils
LIBP += -L$(ROOTSYS)/lib -lRoot -lfreetype -pthread -lpcre $(LISX) $(LIBF)
LIBP += -lQtCore -L$(QTDIR_STATIC)/lib
LIBP += $(LIBD)
