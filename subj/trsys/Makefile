-include $(AMSMKfile)

# Test Ceres Solver
LIBP += -lgflags -lglog -lceres
CXX += -D__CeresSolver__
LD  += -D__CeresSolver__    

TagetDir = $(AMSCoreSubj)/trsys/vdev

Version=18Feb27
ClassLib = $(AMSCoreProd)/${Version}/lib/libClassDef.so
#ClassLib = -L$(AMSCoreProd)/${Version}/lib -lClassDef

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/env
	@mkdir -p $(TagetDir)/obj
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/lst
	@cp -fa $(PWD)/lst/* $(TagetDir)/lst 2> /dev/null
	@cp -fa $(AMSENV) $(TagetDir)/env/amsenv.sh 2> /dev/null
	@cp -fa $(PWD)/src/landau.C $(TagetDir)/src/landau.C 2> /dev/null
	@cp -fa $(PWD)/src/track_fill.C $(TagetDir)/src/track_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/track_fit.C $(TagetDir)/src/track_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/hit_fill.C $(TagetDir)/src/hit_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/hit_fit.C $(TagetDir)/src/hit_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_fill.C $(TagetDir)/src/prop_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_fit.C $(TagetDir)/src/prop_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/fit_fill.C $(TagetDir)/src/fit_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/fit_fit.C $(TagetDir)/src/fit_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_smc_fill.C $(TagetDir)/src/prop_smc_fill.C 2> /dev/null
	#make $(TagetDir)/landau
	make $(TagetDir)/track_fill
	make $(TagetDir)/track_fit
	#make $(TagetDir)/hit_fill
	#make $(TagetDir)/hit_fit
	#make $(TagetDir)/prop_fill
	#make $(TagetDir)/prop_fit
	#make $(TagetDir)/fit_fill
	#make $(TagetDir)/fit_fit
	#make $(TagetDir)/prop_smc_fill

$(TagetDir)/obj/landau.o : $(TagetDir)/src/landau.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/landau.d

$(TagetDir)/landau : $(TagetDir)/obj/landau.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/track_fill.o : $(TagetDir)/src/track_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/track_fill.d

$(TagetDir)/track_fill : $(TagetDir)/obj/track_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/track_fit.o : $(TagetDir)/src/track_fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/track_fit.d

$(TagetDir)/track_fit : $(TagetDir)/obj/track_fit.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/hit_fill.o : $(TagetDir)/src/hit_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/hit_fill.d

$(TagetDir)/hit_fill : $(TagetDir)/obj/hit_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/hit_fit.o : $(TagetDir)/src/hit_fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/hit_fit.d

$(TagetDir)/hit_fit : $(TagetDir)/obj/hit_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/prop_fill.o : $(TagetDir)/src/prop_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/prop_fill.d

$(TagetDir)/prop_fill : $(TagetDir)/obj/prop_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/prop_fit.o : $(TagetDir)/src/prop_fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/prop_fit.d

$(TagetDir)/prop_fit : $(TagetDir)/obj/prop_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/fit_fill.o : $(TagetDir)/src/fit_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/fit_fill.d

$(TagetDir)/fit_fill : $(TagetDir)/obj/fit_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/fit_fit.o : $(TagetDir)/src/fit_fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/fit_fit.d

$(TagetDir)/fit_fit : $(TagetDir)/obj/fit_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/prop_smc_fill.o : $(TagetDir)/src/prop_smc_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/prop_smc_fill.d

$(TagetDir)/prop_smc_fill : $(TagetDir)/obj/prop_smc_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/ceres_test.o : $(TagetDir)/src/ceres_test.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $< -MMD
-include $(TagetDir)/obj/ceres_test.d

$(TagetDir)/ceres_test : $(TagetDir)/obj/ceres_test.o
	$(LD) -o $@ $+ $(LIBP)


clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/landau
	@rm -fv $(TagetDir)/track_fill
	@rm -fv $(TagetDir)/track_fit
	@rm -fv $(TagetDir)/hit_fill
	@rm -fv $(TagetDir)/hit_fit
	@rm -fv $(TagetDir)/prop_fill
	@rm -fv $(TagetDir)/prop_fit
	@rm -fv $(TagetDir)/fit_fill
	@rm -fv $(TagetDir)/fit_fit
	@rm -fv $(TagetDir)/prop_smc_fill
