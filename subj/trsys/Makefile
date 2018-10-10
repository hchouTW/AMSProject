-include $(AMSMKfile)

# Test Ceres Solver
LIBP += -lgflags -lglog -lceres

TagetDir = $(AMSCoreSubj)/trsys/vdev
TagetDir = $(AMSCoreSubj)/trsys/test

#Version=18Sep21
Version=18Oct03
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
	@cp -fa $(PWD)/src/hitq_fit.C $(TagetDir)/src/hitq_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/hitt_fit.C $(TagetDir)/src/hitt_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_fill.C $(TagetDir)/src/prop_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_fit.C $(TagetDir)/src/prop_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/track_fill.C $(TagetDir)/src/track_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/track_fit.C $(TagetDir)/src/track_fit.C 2> /dev/null
	#make $(TagetDir)/time_fill
	#make $(TagetDir)/time_fit
	#make $(TagetDir)/hit_fill
	#make $(TagetDir)/hitq_fit
	#make $(TagetDir)/hitt_fit
	#make $(TagetDir)/prop_fill
	#make $(TagetDir)/prop_fit
	make $(TagetDir)/track_fill
	make $(TagetDir)/track_fit

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

$(TagetDir)/obj/hitq_fit.o : $(TagetDir)/src/hitq_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/hitq_fit.d

$(TagetDir)/hitq_fit : $(TagetDir)/obj/hitq_fit.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/obj/hitt_fit.o : $(TagetDir)/src/hitt_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/hitt_fit.d

$(TagetDir)/hitt_fit : $(TagetDir)/obj/hitt_fit.o
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

$(TagetDir)/obj/track_fill.o : $(TagetDir)/src/track_fill.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/track_fill.d

$(TagetDir)/track_fill : $(TagetDir)/obj/track_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/track_fit.o : $(TagetDir)/src/track_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/track_fit.d

$(TagetDir)/track_fit : $(TagetDir)/obj/track_fit.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)


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
	@rm -fv $(TagetDir)/track_fill
	@rm -fv $(TagetDir)/track_fit
	@rm -fv $(TagetDir)/fit_fit
