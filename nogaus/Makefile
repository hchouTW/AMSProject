UABI = -D_GLIBCXX_USE_CXX11_ABI=0 -DGLOG_NO_ABBREVIATED_SEVERITIES -DGOOGLE_STRIP_LOG=0 # _GLIBCXX_USE_CXX11_ABI std::__c++11::string and std::string difference
PARA = -D_X8664__ -D_DEBUG 

CXX =
LDP =
LDS =

ifeq ($(COMPILER), ICC)
    CXX = icpc -std=c++14 -O2 -Wall $(UABI) $(PARA)
    LDP = $(shell root-config --ld) -std=c++14 $(PARA)
    LDS = $(shell root-config --ld) -std=c++14 $(PARA)
endif

ifeq ($(COMPILER), GCC)
    CXX = g++ -std=c++14 -O2 -Wall $(UABI) $(PARA)
    LDP = $(shell root-config --ld) -std=c++14 $(PARA)
    LDS = $(shell root-config --ld) -std=c++14 $(PARA)
endif

ifeq ($(CXX),)
    $(error Compiler not set correctly)
endif

LIBP = $(shell root-config --libs) -fopenmp -lEG -lRooFit -lRooFitCore -lgflags -lglog -lceres
LIBS = $(shell root-config --libs) -fopenmp -lEG -lRooFit -lRooFitCore -lgflags -lglog -lceres

ProjDir=${PWD}

CoreVersion=20Apr24
CoreDir=${AMSCore}/subj/nogaus/${CoreVersion}

all :
	#make ${ProjDir}/exe/mdst
	#make ${ProjDir}/exe/generator
	make ${ProjDir}/exe/processor
	make ${ProjDir}/exe/fitting
	#make ${ProjDir}/exe/ksMap
	make ${ProjDir}/exe/ksFit
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

${ProjDir}/obj/mdst.o : ${ProjDir}/inc/mdst.C ${ProjDir}/inc/analyzer.C ${ProjDir}/inc/analyzer.h
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/ana.d

${ProjDir}/exe/mdst : ${ProjDir}/obj/mdst.o
	$(LDP) -o $@ $+ $(LIBP) -L${ProjDir}/lib -lClassDef

${ProjDir}/obj/generator.o : ${ProjDir}/inc/generator.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/generator.d

${ProjDir}/exe/generator : ${ProjDir}/obj/generator.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/processor.o : ${ProjDir}/inc/processor.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/processor.d

${ProjDir}/exe/processor : ${ProjDir}/obj/processor.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fitting.o : ${ProjDir}/inc/fitting.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fitting.d

${ProjDir}/exe/fitting : ${ProjDir}/obj/fitting.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/ksMap.o : ${ProjDir}/inc/ksMap.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/ksMap.d

${ProjDir}/exe/ksMap : ${ProjDir}/obj/ksMap.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/ksFit.o : ${ProjDir}/inc/ksFit.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/ksFit.d

${ProjDir}/exe/ksFit : ${ProjDir}/obj/ksFit.o
	$(LDP) -o $@ $+ $(LIBP)

clean :
	@rm -fv ${ProjDir}/obj/*
	@rm -fv ${ProjDir}/exe/*
	@rm -fv ${ProjDir}/log/*
