UABI = -D_GLIBCXX_USE_CXX11_ABI=0 -DGLOG_NO_ABBREVIATED_SEVERITIES -DGOOGLE_STRIP_LOG=0 # _GLIBCXX_USE_CXX11_ABI std::__c++11::string and std::string difference
PARA = -D_X8664__ -D_DEBUG -D_PGTRACK_ -D__ROOTSHAREDLIBRARY__

CXX =
LDP =
LDS =

ifeq ($(COMPILER), ICC)
    CXX = icpc -std=c++14 -O3 -Wall $(UABI) $(PARA)
    LDP = $(shell root-config --ld) -std=c++14 $(PARA)
    LDS = $(shell root-config --ld) -std=c++14 $(PARA)
endif

ifeq ($(COMPILER), GCC)
    CXX = g++ -std=c++14 -O3 -Wall $(UABI) $(PARA)
    LDP = $(shell root-config --ld) -std=c++14 $(PARA)
    LDS = $(shell root-config --ld) -std=c++14 $(PARA)
endif

ifeq ($(CXX),)
    $(error Compiler not set correctly)
endif

LIBP = $(shell root-config --libs) -fopenmp -lEG -lntuple_slc6_PG_dynamic -lmathlib -lgflags -lglog -lceres
LIBS = $(shell root-config --libs) -fopenmp -lEG -lntuple_slc6_PG_dynamic -lmathlib -lgflags -lglog -lceres

ProjDir=${PWD}

CoreVersion=vdev
#CoreVersion=20Jul08
CoreDir=${AMSCore}/mdst/${CoreVersion}

all :
	make ${ProjDir}/lib/libClassDef.so
	make ${ProjDir}/exe/testing
	make ${ProjDir}/exe/mdst
	@rm -rf ${CoreDir}
	@mkdir -p ${CoreDir}
	@mkdir -p ${CoreDir}/env
	@mkdir -p ${CoreDir}/lst
	@mkdir -p ${CoreDir}/lib
	@mkdir -p ${CoreDir}/inc
	@cp -fa ${ProjDir}/env/setup_amsenv.sh ${CoreDir}/env/amsenv.sh
	@cp -fa ${ProjDir}/lst/* ${CoreDir}/lst
	@cp -fa ${ProjDir}/lib/* ${CoreDir}/lib
	@cp -fa ${ProjDir}/inc/* ${CoreDir}/inc
	@cp -fa ${ProjDir}/exe/* ${CoreDir}
	@echo -e "CORE: ${CoreDir}"

${ProjDir}/obj/testing.o : ${ProjDir}/inc/testing.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/testing.d

${ProjDir}/exe/testing : ${ProjDir}/obj/testing.o
	$(LDP) -o $@ $+ $(shell root-config --libs)

${ProjDir}/obj/dict.C : ${ProjDir}/inc/ClassDef.h ${ProjDir}/inc/LinkDef.h
	rootcint -f $@ -c -p $+
	#rootcling -f $@ -c -p $+

${ProjDir}/lib/libClassDef.so : ${ProjDir}/obj/dict.C
	$(CXX) -fPIC -shared $^ -o $@

${ProjDir}/obj/mdst.o : ${ProjDir}/inc/mdst.C ${ProjDir}/inc/selector.C ${ProjDir}/inc/selector.h
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/mdst.d

${ProjDir}/exe/mdst : ${ProjDir}/obj/mdst.o
	$(LDP) -o $@ $+ $(LIBP) -L${ProjDir}/lib -lClassDef

clean :
	@rm -fv ${ProjDir}/obj/*
	@rm -fv ${ProjDir}/exe/*
	@rm -fv ${ProjDir}/log/*
