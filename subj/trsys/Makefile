-include $(AMSMKfile)

# Test Ceres Solver
LIBP += -lgflags -lglog -lceres

TagetDir = $(AMSCoreSubj)/trsys/vdev
TagetDir = $(AMSCoreSubj)/trsys/test3

Version=18Mar23
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
	@cp -fa $(PWD)/src/hit_fit.C $(TagetDir)/src/hit_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_fill.C $(TagetDir)/src/prop_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_fit.C $(TagetDir)/src/prop_fit.C 2> /dev/null
	@cp -fa $(PWD)/src/track_fill.C $(TagetDir)/src/track_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/track_fit.C $(TagetDir)/src/track_fit.C 2> /dev/null
	#make $(TagetDir)/time_fill
	#make $(TagetDir)/time_fit
	#make $(TagetDir)/hit_fill
	#make $(TagetDir)/hit_fit
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

$(TagetDir)/obj/hit_fit.o : $(TagetDir)/src/hit_fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/hit_fit.d

$(TagetDir)/hit_fit : $(TagetDir)/obj/hit_fit.o
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
	@rm -fv $(TagetDir)/hit_fit
	@rm -fv $(TagetDir)/prop_fill
	@rm -fv $(TagetDir)/prop_fit
	@rm -fv $(TagetDir)/track_fill
	@rm -fv $(TagetDir)/track_fit
	@rm -fv $(TagetDir)/fit_fit
