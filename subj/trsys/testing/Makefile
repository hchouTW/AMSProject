-include $(AMSMKfile)

TagetDir = $(AMSCoreSubj)/trsys/testing

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/obj
	@cp -fa $(PWD)/src/prop.C   $(TagetDir)/src/prop.C 2> /dev/null
	make $(TagetDir)/prop

$(TagetDir)/obj/prop.o : $(TagetDir)/src/prop.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/prop.d

$(TagetDir)/prop : $(TagetDir)/obj/prop.o
	$(LD) -o $@ $+ $(LIBS)

clean :
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/prop
