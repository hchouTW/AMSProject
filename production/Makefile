DEFS = -D_PGTRACK_ -D__ROOTSHAREDLIBRARY__
PARA = -D_X8664__ -D_AMSPARALLEL__
#IOPT = -axsse4.2,ssse3 # Only for Intel Core
IOPT =
CXX  = icpc $(IOPT) -std=c++11 $(PARA) -O2 -D_DEBUG -Wall
FLGS = -I$(ROOTSYS)/include -isystem$(ROOTSYS)/include -I$(AMSSRC)/include -isystem$(AMSSRC)/include

OPTS = -g
STAT = -static-intel
LD   = $(shell ${ROOTSYS}/bin/root-config --ld) $(IOPT) $(PARA)
LDS  = $(shell ${ROOTSYS}/bin/root-config --ld) $(IOPT) $(PARA) $(STAT)

LIBA = $(AMSLIB)/libntuple_slc6_PG.a
LIBD = $(AMSLIB)/libntuple_slc6_PG_dynamic.so

CERNDIR = /afs/cern.ch/exp/ams/Offline/CERN/2005
LIBP  = $(shell ${ROOTSYS}/bin/root-config --libs)
LIBP += -lifcore -L$(CERNDIR)/lib -lmathlib
LIBP += -lMinuit -lMinuit2 -lRFIO -lTMVA -lXMLIO -lMLP -lTreePlayer -lMathCore -lMathMore
LIBP += -llzma -Llib -lcrypto -L$(XRDLIB)/lib64 -lXrdClient -lXrdUtils
LIBP += -L$(ROOTSYS)/lib -lRoot -lfreetype -pthread -lpcre $(LISX) $(LIBF)
LIBP += -lQtCore -L${QTDIR_STATIC}/lib
LIBP += $(LIBD)


TagetDir = /afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/vdev

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/lib
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/obj
	@cp -fa $(PWD)/src/ClassDef.C $(TagetDir)/src/ClassDef.C 2> /dev/null
	@cp -fa $(PWD)/src/ClassDef.h $(TagetDir)/src/ClassDef.h 2> /dev/null
	@cp -fa $(PWD)/src/LinkDef.h $(TagetDir)/src/LinkDef.h 2> /dev/null
	@cp -fa $(PWD)/src/YiProdNtuple.h $(TagetDir)/src/YiProdNtuple.h 2> /dev/null
	@cp -fa $(PWD)/src/YiProdNtuple.tcc $(TagetDir)/src/YiProdNtuple.tcc 2> /dev/null
	@cp -fa $(PWD)/src/YiProdNtuple.C $(TagetDir)/src/YiProdNtuple.C 2> /dev/null
	make $(TagetDir)/obj/dict.o
	make $(TagetDir)/lib/ClassDef.so
	make $(TagetDir)/YiProdNtuple

$(TagetDir)/src/dict.C : $(TagetDir)/src/ClassDef.h $(TagetDir)/src/ClassDef.C $(TagetDir)/src/LinkDef.h
	rootcling -f $@ -c -p $(FLGS) $+

$(TagetDir)/obj/dict.o : $(TagetDir)/src/dict.C
	$(CXX) -c -o $@ $(FLGS) $<

$(TagetDir)/lib/ClassDef.so : $(TagetDir)/src/dict.C
	$(CXX) -fPIC -shared $^ -o $@ $(FLGS)
	@cp -fa $(TagetDir)/src/dict_rdict.pcm $(TagetDir)/lib/dict_rdict.pcm

$(TagetDir)/obj/YiProdNtuple.o : $(TagetDir)/src/YiProdNtuple.C $(TagetDir)/src/YiProdNtuple.tcc $(TagetDir)/src/YiProdNtuple.h
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD -include $(TagetDir)/obj/YiProdNtuple.d

$(TagetDir)/YiProdNtuple : $(TagetDir)/obj/YiProdNtuple.o $(TagetDir)/lib/ClassDef.so
	$(LD) -o $@ $+ $(LIBP)

clean :
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lib/*
	@rm -fv $(TagetDir)/YiProdNtuple
