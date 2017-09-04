-include $(AMSMKfile)

ClassLib = $(AMSCoreProd)/vdev/lib/ClassDef.so
TagetDir = $(AMSCoreSubj)/template

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/env
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/obj
	@mkdir -p $(TagetDir)/lst
	@cp -fa $(PWD)/lst/* $(TagetDir)/lst 2> /dev/null
	@cp -fa $(AMSENV) $(TagetDir)/env/amsenv.sh 2> /dev/null
	@cp -fa $(PWD)/src/YiSlcNtuple.h   $(TagetDir)/src/YiSlcNtuple.h 2> /dev/null
	@cp -fa $(PWD)/src/YiSlcNtuple.tcc $(TagetDir)/src/YiSlcNtuple.tcc 2> /dev/null
	@cp -fa $(PWD)/src/YiSlcNtuple.C   $(TagetDir)/src/YiSlcNtuple.C 2> /dev/null
	make $(TagetDir)/YiSlcNtuple

$(TagetDir)/obj/YiSlcNtuple.o : $(TagetDir)/src/YiSlcNtuple.C $(TagetDir)/src/YiSlcNtuple.tcc $(TagetDir)/src/YiSlcNtuple.h
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/YiSlcNtuple.d

$(TagetDir)/YiSlcNtuple : $(TagetDir)/obj/YiSlcNtuple.o $(ClassLib)
	$(LD) -o $@ $+ $(LIBS)

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/YiSlcNtuple
