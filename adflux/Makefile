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
    CXX = g++ -std=c++14 -O2 -Wall $(UABI) $(PARA) -Wno-sign-compare -Wno-unused-variable
    LDP = $(shell root-config --ld) -std=c++14 $(PARA)
    LDS = $(shell root-config --ld) -std=c++14 $(PARA)
endif

ifeq ($(CXX),)
    $(error Compiler not set correctly)
endif

LIBP = $(shell root-config --libs) -fopenmp -lEG -lRooFit -lRooFitCore -lgflags -lglog -lceres
LIBS = $(shell root-config --libs) -fopenmp -lEG -lRooFit -lRooFitCore -lgflags -lglog -lceres

ProjDir=${PWD}

CoreVersion=vdev
CoreVersion=20Jan15
CoreDir=${AMSCore}/subj/adflux/${CoreVersion}

all :
	#make ${ProjDir}/exe/lvtme
	make ${ProjDir}/exe/adflux_acc
	make ${ProjDir}/exe/adflux_flx
	@rm -rf ${CoreDir}
	@mkdir -p ${CoreDir}
	@mkdir -p ${CoreDir}/env
	@mkdir -p ${CoreDir}/lst
	@mkdir -p ${CoreDir}/lib
	@mkdir -p ${CoreDir}/inc
	@cp -fa ${ProjDir}/env/setup_amsenv6.sh ${CoreDir}/env/amsenv.sh
	@cp -fa ${ProjDir}/lst/* ${CoreDir}/lst
	@cp -fa ${ProjDir}/lib/* ${CoreDir}/lib
	@cp -fa ${ProjDir}/inc/* ${CoreDir}/inc
	@cp -fa ${ProjDir}/exe/* ${CoreDir}
	@echo -e "CORE: ${CoreDir}"

${ProjDir}/obj/lvtme.o : ${ProjDir}/inc/lvtme.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/lvtme.d

${ProjDir}/exe/lvtme : ${ProjDir}/obj/lvtme.o
	$(LDP) -o $@ $+ $(LIBP) -L${ProjDir}/lib -lClassDefRun

${ProjDir}/obj/adflux_acc.o : ${ProjDir}/inc/adflux_acc.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/adflux_acc.d

${ProjDir}/exe/adflux_acc : ${ProjDir}/obj/adflux_acc.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/adflux_flx.o : ${ProjDir}/inc/adflux_flx.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/adflux_flx.d

${ProjDir}/exe/adflux_flx : ${ProjDir}/obj/adflux_flx.o
	$(LDP) -o $@ $+ $(LIBP)

clean :
	@rm -fv ${ProjDir}/obj/*
	@rm -fv ${ProjDir}/exe/*
	@rm -fv ${ProjDir}/log/*
