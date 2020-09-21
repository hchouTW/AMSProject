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
CoreDir=${AMSCore}/subj/apflux/${CoreVersion}

all :
	make ${ProjDir}/exe/prflux_acc
	make ${ProjDir}/exe/prflux_flx
	make ${ProjDir}/exe/apflux_acc
	make ${ProjDir}/exe/apflux_app
	make ${ProjDir}/exe/apflux_tme
	make ${ProjDir}/exe/phys_tme
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

${ProjDir}/obj/prflux_acc.o : ${ProjDir}/inc/prflux_acc.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/prflux_acc.d

${ProjDir}/exe/prflux_acc : ${ProjDir}/obj/prflux_acc.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/prflux_flx.o : ${ProjDir}/inc/prflux_flx.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/prflux_flx.d

${ProjDir}/exe/prflux_flx : ${ProjDir}/obj/prflux_flx.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/apflux_acc.o : ${ProjDir}/inc/apflux_acc.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apflux_acc.d

${ProjDir}/exe/apflux_acc : ${ProjDir}/obj/apflux_acc.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/apflux_app.o : ${ProjDir}/inc/apflux_app.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apflux_app.d

${ProjDir}/exe/apflux_app : ${ProjDir}/obj/apflux_app.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/apflux_tme.o : ${ProjDir}/inc/apflux_tme.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apflux_tme.d

${ProjDir}/exe/apflux_tme : ${ProjDir}/obj/apflux_tme.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/phys_tme.o : ${ProjDir}/inc/phys_tme.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/phys_tme.d

${ProjDir}/exe/phys_tme : ${ProjDir}/obj/phys_tme.o
	$(LDP) -o $@ $+ $(LIBP)

clean :
	@rm -fv ${ProjDir}/obj/*
	@rm -fv ${ProjDir}/exe/*
	@rm -fv ${ProjDir}/log/*
