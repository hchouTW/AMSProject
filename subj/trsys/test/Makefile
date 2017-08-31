-include $(AMSMKfile)

TagetDir = $(AMSCoreSubj)/trsys/test

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/obj
	#@cp -fa $(PWD)/src/testing.h   $(TagetDir)/src/testing.h 2> /dev/null
	@cp -fa $(PWD)/src/testing.C   $(TagetDir)/src/testing.C 2> /dev/null
	make $(TagetDir)/testing

#$(TagetDir)/obj/testing.o : $(TagetDir)/src/testing.C $(TagetDir)/src/testing.h
$(TagetDir)/obj/testing.o : $(TagetDir)/src/testing.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/testing.d

$(TagetDir)/testing : $(TagetDir)/obj/testing.o
	$(LD) -o $@ $+ $(LIBS)

clean :
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/testing
