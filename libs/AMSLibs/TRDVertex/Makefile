CXX=$(shell root-config --cxx)
MARCH=$(shell root-config --arch)
ROOTLIB=$(shell root-config --libs)

ifeq "$(CXX)" "icc"
CFLIB=-lifcore
CF=ifort
else
CFLIB=-lgfortran
CF=gfortran
endif

WFLAGS= -openmp
# -std=c++11

#Offline=/afs/cern.ch/work/w/wxu/public/AMS/
AMSSRC=$(Offline)/vdev
#B900_patches/

#AMSSRC_LOCAL=/afs/cern.ch/work/w/wxu/public/AMS/vdev/include
#CPPFLAGS=-D_PGTRACK_ -I. -Ilib -I$(AMSSRC_LOCAL) -I$(ROOTSYS)/include





##AMSSRC=/afs/cern.ch/work/z/zweng/20150729/AMS
AMSLIB=$(AMSSRC)/lib/linuxx8664icc5.34/
#AMSLIB=$(AMSSRC)/lib/linuxx8664gcc5.34/
#Offline=/afs/cern.ch/work/z/zwengams/public/AMS/
#AMSSRC=$(Offline)/vdev


CPPFLAGS=-D_PGTRACK_ -D__ROOTSHAREDLIBRARY__ -I$(AMSSRC)/include -I$(ROOTSYS)/include
CXXFLAGS=   -g $(WFLAGS) $(CPPFLAGS) 
#CXXFLAGS= -debug -O0  $(WFLAGS) $(CPPFLAGS) 


NAGDIR=/afs/cern.ch/exp/ams/Offline/CERN/NagLib
CERNDIRLIB=/afs/cern.ch/exp/ams/Offline/CERN/2005/lib/
#CERNDIRLIB=/afs/cern.ch/exp/ams/Offline/CERN/2005.gcc64.44/lib/

all: libTRDVertex.a 

#main_new


SRCS= Cell.C Segment.C HitState.C TRD2DHit.C Line.C Vertex2D.C Vertex3D.C TRDVertex.C rootdict.C


MYOBJS := $(SRCS:$(SRC)%.C=$(BIN)%.o)

$(BIN)%.o: $(SRC)%.C
	@echo Compiling  $< ...	
	$(CXX) -c $(CXXFLAGS)  $< -o $@

#main_new: main.C $(MYOBJS) rootdict.C
#	$(CXX)  $^ -o main_new  $(CXXFLAGS)  $(CFLIB)  $(ROOTLIB) -lMinuit -lNetx -lProof -lProofPlayer -lTMVA -lXMLIO -lRFIO -lTreePlayer -lMLP  -lFoam -lMathCore   -L$(AMSLIB) -L$(NAGDIR#) -lntuple_slc6_PG $(ROOTLIB) -lMinuit -lNetx -lProof -lProofPlayer -lTMVA -lXMLIO -lRFIO -lTreePlayer -lMLP  -lFoam -lMathCore  -lnag64  -L$(CERNDIRLIB) -lmathlib
#

libTRDVertex.a: $(MYOBJS) 
	@echo Compilling $< ......
	ar rcs $@ $^

rootdict.C: HitState.h SimpleHitState.h linkdef.h
	$(ROOTSYS)/bin/rootcint -f rootdict.C -c -p $(CPPFLAGS) $^

clean:
	rm -f libTRDVertex.a rootdict*.C *.o


