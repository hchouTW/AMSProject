-include $(AMSMKfile)

# Test Ceres Solver
LIBP += -lgflags -lglog -lceres

TagetDir = $(AMSCoreSubj)/appStudy/vdev
TagetDir = $(AMSCoreSubj)/appStudy/test

Version=19Jan21
ClassLib = $(AMSCoreProd)/${Version}/lib/libClassDef.so

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/env
	@mkdir -p $(TagetDir)/obj
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/lst
	@cp -fa $(PWD)/lst/* $(TagetDir)/lst 2> /dev/null
	@cp -fa $(AMSENV) $(TagetDir)/env/amsenv.sh 2> /dev/null
	@cp -fa $(PWD)/src/fill.C $(TagetDir)/src/fill.C 2> /dev/null
	@cp -fa $(PWD)/src/fit.C $(TagetDir)/src/fit.C 2> /dev/null
	make $(TagetDir)/fill
	make $(TagetDir)/fit

$(TagetDir)/obj/fill.o : $(TagetDir)/src/fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/fill.d

$(TagetDir)/fill : $(TagetDir)/obj/fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/fit.o : $(TagetDir)/src/fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/fit.d

$(TagetDir)/fit : $(TagetDir)/obj/fit.o
	$(LD) -o $@ $+ $(LIBP)

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/fill
	@rm -fv $(TagetDir)/fit
