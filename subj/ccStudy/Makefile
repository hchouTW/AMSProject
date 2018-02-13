-include $(AMSMKfile)

TagetDir = $(AMSCoreSubj)/ccStudy/vdev

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
	@cp -fa $(PWD)/src/cc_fill.C $(TagetDir)/src/cc_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/cc_fit.C $(TagetDir)/src/cc_fit.C 2> /dev/null
	make $(TagetDir)/cc_fill
	make $(TagetDir)/cc_fit

$(TagetDir)/obj/cc_fill.o : $(TagetDir)/src/cc_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/cc_fill.d

$(TagetDir)/cc_fill : $(TagetDir)/obj/cc_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

$(TagetDir)/obj/cc_fit.o : $(TagetDir)/src/cc_fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/cc_fit.d

$(TagetDir)/cc_fit : $(TagetDir)/obj/cc_fit.o
	$(LD) -o $@ $+ $(LIBP)

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/cc_fill
	@rm -fv $(TagetDir)/cc_fit
