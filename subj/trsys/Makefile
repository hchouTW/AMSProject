-include $(AMSMKfile)

TagetDir = $(AMSCoreSubj)/trsys/vdev

#Version=18Jan29
Version=18Feb05
ClassLib = $(AMSCoreProd)/${Version}/lib/libClassDef.so

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/env
	@mkdir -p $(TagetDir)/obj
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/lst
	@cp -fa $(PWD)/lst/* $(TagetDir)/lst 2> /dev/null
	@cp -fa $(AMSENV) $(TagetDir)/env/amsenv.sh 2> /dev/null
	@cp -fa $(PWD)/src/hit_ams02_fill.C $(TagetDir)/src/hit_ams02_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_ams02_fill.C $(TagetDir)/src/prop_ams02_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_ams02_fit.C $(TagetDir)/src/prop_ams02_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/fit_ams02_fill.C $(TagetDir)/src/fit_ams02_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/fit_ams02_fit.C $(TagetDir)/src/fit_ams02_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_smc_fill.C $(TagetDir)/src/prop_smc_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/ceres_test.C $(TagetDir)/src/ceres_test.C 2> /dev/null
	#make $(TagetDir)/hit_ams02_fill
	#make $(TagetDir)/prop_ams02_fill
	#make $(TagetDir)/prop_ams02_fit
	make $(TagetDir)/fit_ams02_fill
	#make $(TagetDir)/fit_ams02_fit
	#make $(TagetDir)/ceres_test
	#make $(TagetDir)/prop_smc_fill

$(TagetDir)/obj/hit_ams02_fill.o : $(TagetDir)/src/hit_ams02_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/hit_ams02_fill.d

$(TagetDir)/hit_ams02_fill : $(TagetDir)/obj/hit_ams02_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/prop_ams02_fill.o : $(TagetDir)/src/prop_ams02_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/prop_ams02_fill.d

$(TagetDir)/prop_ams02_fill : $(TagetDir)/obj/prop_ams02_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/prop_ams02_fit.o : $(TagetDir)/src/prop_ams02_fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/prop_ams02_fit.d

$(TagetDir)/prop_ams02_fit : $(TagetDir)/obj/prop_ams02_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/fit_ams02_fill.o : $(TagetDir)/src/fit_ams02_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/fit_ams02_fill.d

$(TagetDir)/fit_ams02_fill : $(TagetDir)/obj/fit_ams02_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/fit_ams02_fit.o : $(TagetDir)/src/fit_ams02_fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/fit_ams02_fit.d

$(TagetDir)/fit_ams02_fit : $(TagetDir)/obj/fit_ams02_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/prop_smc_fill.o : $(TagetDir)/src/prop_smc_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/prop_smc_fill.d

$(TagetDir)/prop_smc_fill : $(TagetDir)/obj/prop_smc_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/ceres_test.o : $(TagetDir)/src/ceres_test.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/ceres_test.d

$(TagetDir)/ceres_test : $(TagetDir)/obj/ceres_test.o
	$(LD) -o $@ $+ $(LIBP)


clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/prop
	@rm -fv $(TagetDir)/fit
	@rm -fv $(TagetDir)/hit_ams02_fill
	@rm -fv $(TagetDir)/prop_ams02_fill
	@rm -fv $(TagetDir)/prop_ams02_fit
	@rm -fv $(TagetDir)/fit_ams02_fill
	@rm -fv $(TagetDir)/fit_ams02_fit
	@rm -fv $(TagetDir)/prop_smc_fill
