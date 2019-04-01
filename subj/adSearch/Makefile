-include $(AMSMKfile)

# Test Ceres Solver
CXX  += -DGLOG_NO_ABBREVIATED_SEVERITIES
LIBP += -lgflags -lglog -lceres
LIBS += -lgflags -lglog -lceres

TagetDir = $(AMSCoreSubj)/adSearch/vdev
TagetDir = $(AMSCoreSubj)/adSearch/test

Version=19Apr01
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
	@cp -fa $(PWD)/src/fill.C $(TagetDir)/src/fill.C 2> /dev/null
	@cp -fa $(PWD)/src/fit.C $(TagetDir)/src/fit.C 2> /dev/null
	make $(TagetDir)/fill
	make $(TagetDir)/fit

$(TagetDir)/obj/fill.o : $(TagetDir)/src/fill.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/fill.d

$(TagetDir)/fill : $(TagetDir)/obj/fill.o
	$(LDP) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/fit.o : $(TagetDir)/src/fit.C
	$(CXX) -c -o $@ $(FLGS) $< -MMD
-include $(TagetDir)/obj/fit.d

$(TagetDir)/fit : $(TagetDir)/obj/fit.o
	$(LDP) -o $@ $+ $(LIBP) $(ClassLib)

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/fill
	@rm -fv $(TagetDir)/fit
