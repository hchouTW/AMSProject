-include $(AMSMKfile)

# Test Ceres Solver
LIBP += -lgflags -lglog -lceres

TagetDir = $(AMSCoreSubj)/trsys/vdev
TagetDir = $(AMSCoreSubj)/trsys/test

#Version=18Oct17
Version=19Jan09
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
	@cp -fa $(PWD)/src/time_fill.C $(TagetDir)/src/time_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/time_fit.C $(TagetDir)/src/time_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/hit_fill.C $(TagetDir)/src/hit_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/hitTDq_fit.C $(TagetDir)/src/hitTDq_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/hitt_fit.C $(TagetDir)/src/hitt_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/hitc_fit.C $(TagetDir)/src/hitc_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_fill.C $(TagetDir)/src/prop_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_fit.C $(TagetDir)/src/prop_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/bta_fill.C $(TagetDir)/src/bta_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/bta_fit.C $(TagetDir)/src/bta_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/track_fill.C $(TagetDir)/src/track_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/track_fit.C $(TagetDir)/src/track_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/mu_fill.C $(TagetDir)/src/mu_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/mu_fit.C $(TagetDir)/src/mu_fit.C 2> /dev/null
	#make $(TagetDir)/time_fill
	#make $(TagetDir)/time_fit
	#make $(TagetDir)/hit_fill
	make $(TagetDir)/hitTDq_fit
	#make $(TagetDir)/hitt_fit
	#make $(TagetDir)/hitc_fit
	#make $(TagetDir)/prop_fill
	#make $(TagetDir)/prop_fit
	#make $(TagetDir)/bta_fill
	#make $(TagetDir)/bta_fit
	#make $(TagetDir)/track_fill
	#make $(TagetDir)/track_fit
	#make $(TagetDir)/mu_fill
	#make $(TagetDir)/mu_fit

$(TagetDir)/obj/time_fill.o : $(TagetDir)/src/time_fill.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/time_fill.d

$(TagetDir)/time_fill : $(TagetDir)/obj/time_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/time_fit.o : $(TagetDir)/src/time_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/time_fit.d

$(TagetDir)/time_fit : $(TagetDir)/obj/time_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/hit_fill.o : $(TagetDir)/src/hit_fill.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/hit_fill.d

$(TagetDir)/hit_fill : $(TagetDir)/obj/hit_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/hitTDq_fit.o : $(TagetDir)/src/hitTDq_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/hitTDq_fit.d

$(TagetDir)/hitTDq_fit : $(TagetDir)/obj/hitTDq_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/hitt_fit.o : $(TagetDir)/src/hitt_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/hitt_fit.d

$(TagetDir)/hitt_fit : $(TagetDir)/obj/hitt_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/hitc_fit.o : $(TagetDir)/src/hitc_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/hitc_fit.d

$(TagetDir)/hitc_fit : $(TagetDir)/obj/hitc_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/prop_fill.o : $(TagetDir)/src/prop_fill.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/prop_fill.d

$(TagetDir)/prop_fill : $(TagetDir)/obj/prop_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/prop_fit.o : $(TagetDir)/src/prop_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/prop_fit.d

$(TagetDir)/prop_fit : $(TagetDir)/obj/prop_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/bta_fill.o : $(TagetDir)/src/bta_fill.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/bta_fill.d

$(TagetDir)/bta_fill : $(TagetDir)/obj/bta_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/bta_fit.o : $(TagetDir)/src/bta_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/bta_fit.d

$(TagetDir)/bta_fit : $(TagetDir)/obj/bta_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/track_fill.o : $(TagetDir)/src/track_fill.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/track_fill.d

$(TagetDir)/track_fill : $(TagetDir)/obj/track_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/track_fit.o : $(TagetDir)/src/track_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/track_fit.d

$(TagetDir)/track_fit : $(TagetDir)/obj/track_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/mu_fill.o : $(TagetDir)/src/mu_fill.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/mu_fill.d

$(TagetDir)/mu_fill : $(TagetDir)/obj/mu_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/mu_fit.o : $(TagetDir)/src/mu_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/mu_fit.d

$(TagetDir)/mu_fit : $(TagetDir)/obj/mu_fit.o
	$(LD) -o $@ $+ $(LIBP)

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/time_fill
	@rm -fv $(TagetDir)/time_fit
	@rm -fv $(TagetDir)/hit_fill
	@rm -fv $(TagetDir)/hitq_fit
	@rm -fv $(TagetDir)/hitt_fit
	@rm -fv $(TagetDir)/prop_fill
	@rm -fv $(TagetDir)/prop_fit
	@rm -fv $(TagetDir)/bta_fill
	@rm -fv $(TagetDir)/bta_fit
	@rm -fv $(TagetDir)/track_fill
	@rm -fv $(TagetDir)/track_fit
	@rm -fv $(TagetDir)/mu_fill
	@rm -fv $(TagetDir)/mu_fit
	@rm -fv $(TagetDir)/fit_fit
