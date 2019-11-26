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

LIBP = $(shell root-config --libs) -fopenmp -lEG -lgflags -lglog -lceres
LIBS = $(shell root-config --libs) -fopenmp -lEG -lgflags -lglog -lceres

ProjDir=${PWD}

CoreVersion=vdev
#CoreVersion=19Jul08
CoreDir=${AMSCore}/subj/trsys/${CoreVersion}

all :
	make ${ProjDir}/exe/mdst
	make ${ProjDir}/exe/fit_mscat
	make ${ProjDir}/exe/fit_eloss
	make ${ProjDir}/exe/fit_geom
	make ${ProjDir}/exe/fit_vel_tf
	make ${ProjDir}/exe/fit_vel_rh
	make ${ProjDir}/exe/fit_phys_tf
	make ${ProjDir}/exe/fit_phys_rh
	make ${ProjDir}/exe/fit_mutr_tf
	make ${ProjDir}/exe/fit_mutr_rh
	make ${ProjDir}/exe/fit_merge
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

${ProjDir}/obj/fit_mscat.o : ${ProjDir}/inc/fit_mscat.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_mscat.d

${ProjDir}/exe/fit_mscat : ${ProjDir}/obj/fit_mscat.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_eloss.o : ${ProjDir}/inc/fit_eloss.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_eloss.d

${ProjDir}/exe/fit_eloss : ${ProjDir}/obj/fit_eloss.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_geom.o : ${ProjDir}/inc/fit_geom.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_geom.d

${ProjDir}/exe/fit_geom : ${ProjDir}/obj/fit_geom.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_vel_tf.o : ${ProjDir}/inc/fit_vel_tf.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_vel_tf.d

${ProjDir}/exe/fit_vel_tf : ${ProjDir}/obj/fit_vel_tf.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_vel_rh.o : ${ProjDir}/inc/fit_vel_rh.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_vel_rh.d

${ProjDir}/exe/fit_vel_rh : ${ProjDir}/obj/fit_vel_rh.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_phys_tf.o : ${ProjDir}/inc/fit_phys_tf.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_phys_tf.d

${ProjDir}/exe/fit_phys_tf : ${ProjDir}/obj/fit_phys_tf.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_phys_rh.o : ${ProjDir}/inc/fit_phys_rh.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_phys_rh.d

${ProjDir}/exe/fit_phys_rh : ${ProjDir}/obj/fit_phys_rh.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_mutr_tf.o : ${ProjDir}/inc/fit_mutr_tf.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_mutr_tf.d

${ProjDir}/exe/fit_mutr_tf : ${ProjDir}/obj/fit_mutr_tf.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_mutr_rh.o : ${ProjDir}/inc/fit_mutr_rh.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_mutr_rh.d

${ProjDir}/exe/fit_mutr_rh : ${ProjDir}/obj/fit_mutr_rh.o
	$(LDP) -o $@ $+ $(LIBP)

${ProjDir}/obj/fit_merge.o : ${ProjDir}/inc/fit_merge.C
	$(CXX) -c -o $@ -g $< -MMD
-include ${ProjDir}/obj/fit_merge.d

${ProjDir}/exe/fit_merge : ${ProjDir}/obj/fit_merge.o
	$(LDP) -o $@ $+ $(LIBP)

clean :
	@rm -fv ${ProjDir}/obj/*
	@rm -fv ${ProjDir}/exe/*
	@rm -fv ${ProjDir}/log/*
