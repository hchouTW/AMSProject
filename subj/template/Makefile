-include $(AMSMKfile)
#CXX = $(ICCCXX)
CXX = $(GCCCXX)

#ClassLib = /afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production/V.2017Jul27/lib/ClassDef.so
#TagetDir = /afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/analysis/template

ClassLib = $(AMSCoreProd)/vdev/lib/ClassDef.so
TagetDir = $(AMSCoreSubj)/template

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/obj
	@cp -fa $(PWD)/src/YiSlcNtuple.h   $(TagetDir)/src/YiSlcNtuple.h 2> /dev/null
	@cp -fa $(PWD)/src/YiSlcNtuple.tcc $(TagetDir)/src/YiSlcNtuple.tcc 2> /dev/null
	@cp -fa $(PWD)/src/YiSlcNtuple.C   $(TagetDir)/src/YiSlcNtuple.C 2> /dev/null
	make $(TagetDir)/YiSlcNtuple

$(TagetDir)/obj/YiSlcNtuple.o : $(TagetDir)/src/YiSlcNtuple.C $(TagetDir)/src/YiSlcNtuple.tcc $(TagetDir)/src/YiSlcNtuple.h
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/YiSlcNtuple.d

$(TagetDir)/YiSlcNtuple : $(TagetDir)/obj/YiSlcNtuple.o $(ClassLib)
	$(LD) -o $@ $+ $(LIBP)

clean :
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/YiSlcNtuple
