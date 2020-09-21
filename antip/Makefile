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

CoreVersion=vdev
CoreVersion=19Nov03
CoreDir=${AMSCore}/subj/antip/${CoreVersion}

all :
	#make ${ProjDir}/exe/mdst
	#make ${ProjDir}/exe/issprof
	#make ${ProjDir}/exe/physana_l
	#make ${ProjDir}/exe/physana_m
	#make ${ProjDir}/exe/physana_i
	#make ${ProjDir}/exe/physana_hl1
	#make ${ProjDir}/exe/physana_hl9
	#make ${ProjDir}/exe/physana_hfs
	#make ${ProjDir}/exe/physana_all
	make ${ProjDir}/exe/physana_time
	#make ${ProjDir}/exe/physana_mes
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

${ProjDir}/obj/issprof.o : ${ProjDir}/inc/issprof.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/issprof.d

${ProjDir}/exe/issprof : ${ProjDir}/obj/issprof.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_l.o : ${ProjDir}/inc/physana_l.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_l.d

${ProjDir}/exe/physana_l : ${ProjDir}/obj/physana_l.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_m.o : ${ProjDir}/inc/physana_m.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_m.d

${ProjDir}/exe/physana_m : ${ProjDir}/obj/physana_m.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_i.o : ${ProjDir}/inc/physana_i.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_i.d

${ProjDir}/exe/physana_i : ${ProjDir}/obj/physana_i.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_hl1.o : ${ProjDir}/inc/physana_hl1.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_hl1.d

${ProjDir}/exe/physana_hl1 : ${ProjDir}/obj/physana_hl1.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_hl9.o : ${ProjDir}/inc/physana_hl9.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_hl9.d

${ProjDir}/exe/physana_hl9 : ${ProjDir}/obj/physana_hl9.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_hfs.o : ${ProjDir}/inc/physana_hfs.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_hfs.d

${ProjDir}/exe/physana_hfs : ${ProjDir}/obj/physana_hfs.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_all.o : ${ProjDir}/inc/physana_all.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_all.d

${ProjDir}/exe/physana_all : ${ProjDir}/obj/physana_all.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_time.o : ${ProjDir}/inc/physana_time.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_time.d

${ProjDir}/exe/physana_time : ${ProjDir}/obj/physana_time.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/physana_mes.o : ${ProjDir}/inc/physana_mes.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/physana_mes.d

${ProjDir}/exe/physana_mes : ${ProjDir}/obj/physana_mes.o
	$(LDP) -o $@ $+ $(LIBP)

clean :
	@rm -fv ${ProjDir}/obj/*
	@rm -fv ${ProjDir}/exe/*
	@rm -fv ${ProjDir}/log/*
