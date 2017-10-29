-include $(AMSMKfile)

Version=17SEP27
TagetDir = $(AMSCoreSubj)/trsys/testing/vdev
#ClassLib = $(AMSCoreProd)/${Version}/lib/ClassDef.so
ClassLib = /data1/hchou/${Version}/lib/ClassDef.so

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/env
	@mkdir -p $(TagetDir)/obj
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/lst
	@cp -fa $(PWD)/lst/* $(TagetDir)/lst 2> /dev/null
	@cp -fa $(AMSENV) $(TagetDir)/env/amsenv.sh 2> /dev/null
	#@cp -fa $(PWD)/src/prop.C $(TagetDir)/src/prop.C 2> /dev/null
	#@cp -fa $(PWD)/src/fit.C $(TagetDir)/src/fit.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_ams02_fill.C $(TagetDir)/src/prop_ams02_fill.C 2> /dev/null
	@cp -fa $(PWD)/src/prop_ams02_fit.C $(TagetDir)/src/prop_ams02_fit.C 2> /dev/null
	#make $(TagetDir)/prop
	#make $(TagetDir)/fit
	make $(TagetDir)/prop_ams02_fill
	make $(TagetDir)/prop_ams02_fit

$(TagetDir)/obj/prop.o : $(TagetDir)/src/prop.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/prop.d

$(TagetDir)/prop : $(TagetDir)/obj/prop.o
	$(LD) -o $@ $+ $(LIBS)

$(TagetDir)/obj/fit.o : $(TagetDir)/src/fit.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/fit.d

$(TagetDir)/fit : $(TagetDir)/obj/fit.o
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
	$(LD) -o $@ $+ $(LIBS)

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/prop
	@rm -fv $(TagetDir)/fit
	@rm -fv $(TagetDir)/prop_ams02_fill
	@rm -fv $(TagetDir)/prop_ams02_fit
