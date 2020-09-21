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
	make ${ProjDir}/exe/apflux_acc
	make ${ProjDir}/exe/apflux_flx
	make ${ProjDir}/exe/apflux_tme
	make ${ProjDir}/exe/apflux_rlt
	make ${ProjDir}/exe/apflux_phy
	make ${ProjDir}/exe/psphys_tme
	make ${ProjDir}/exe/apphys_tme
	#make ${ProjDir}/exe/phys_tme
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

${ProjDir}/obj/apflux_acc.o : ${ProjDir}/inc/apflux_acc.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apflux_acc.d

${ProjDir}/exe/apflux_acc : ${ProjDir}/obj/apflux_acc.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/apflux_flx.o : ${ProjDir}/inc/apflux_flx.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apflux_flx.d

${ProjDir}/exe/apflux_flx : ${ProjDir}/obj/apflux_flx.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/apflux_tme.o : ${ProjDir}/inc/apflux_tme.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apflux_tme.d

${ProjDir}/exe/apflux_tme : ${ProjDir}/obj/apflux_tme.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/apflux_rlt.o : ${ProjDir}/inc/apflux_rlt.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apflux_rlt.d

${ProjDir}/exe/apflux_rlt : ${ProjDir}/obj/apflux_rlt.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/apflux_phy.o : ${ProjDir}/inc/apflux_phy.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apflux_phy.d

${ProjDir}/exe/apflux_phy : ${ProjDir}/obj/apflux_phy.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/psphys_tme.o : ${ProjDir}/inc/psphys_tme.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/psphys_tme.d

${ProjDir}/exe/psphys_tme : ${ProjDir}/obj/psphys_tme.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/apphys_tme.o : ${ProjDir}/inc/apphys_tme.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/apphys_tme.d

${ProjDir}/exe/apphys_tme : ${ProjDir}/obj/apphys_tme.o
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
