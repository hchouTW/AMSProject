-include $(AMSMKfile)

TagetDir = $(AMSCoreSubj)/trsys/testing/vdev

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/env
	@mkdir -p $(TagetDir)/obj
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/lst
	@cp -fa $(PWD)/lst/* $(TagetDir)/lst 2> /dev/null
	@cp -fa $(AMSENV) $(TagetDir)/env/amsenv.sh 2> /dev/null
	@cp -fa $(PWD)/src/prop.C $(TagetDir)/src/prop.C 2> /dev/null
	make $(TagetDir)/prop

$(TagetDir)/obj/prop.o : $(TagetDir)/src/prop.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/prop.d

$(TagetDir)/prop : $(TagetDir)/obj/prop.o
	$(LD) -o $@ $+ $(LIBS)

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/prop
