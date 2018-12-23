-include $(AMSMKfile)

CXX += -DGLOG_NO_ABBREVIATED_SEVERITIES
LD += -DGLOG_NO_ABBREVIATED_SEVERITIES
LIBP += -lgflags -lglog -lceres

Version = vdev
#Version = 18Oct17
TagetDir = $(AMSCoreProd)/$(Version)

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/env
	@mkdir -p $(TagetDir)/lib
	@mkdir -p $(TagetDir)/obj
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/lst
	@cp -fa $(PWD)/lst/* $(TagetDir)/lst 2> /dev/null
	@cp -fa $(AMSENV) $(TagetDir)/env/amsenv.sh 2> /dev/null
	@cp -fa $(PWD)/src/testing.C $(TagetDir)/src/testing.C 2> /dev/null
	@cp -fa $(PWD)/src/ClassDef.C $(TagetDir)/src/ClassDef.C 2> /dev/null
	@cp -fa $(PWD)/src/ClassDef.h $(TagetDir)/src/ClassDef.h 2> /dev/null
	@cp -fa $(PWD)/src/LinkDef.h $(TagetDir)/src/LinkDef.h 2> /dev/null
	@cp -fa $(PWD)/src/YiProdNtuple.h $(TagetDir)/src/YiProdNtuple.h 2> /dev/null
	@cp -fa $(PWD)/src/YiProdNtuple.tcc $(TagetDir)/src/YiProdNtuple.tcc 2> /dev/null
	@cp -fa $(PWD)/src/YiProdNtuple.C $(TagetDir)/src/YiProdNtuple.C 2> /dev/null
	#make $(TagetDir)/testing
	make $(TagetDir)/obj/dict.o
	make $(TagetDir)/lib/libClassDef.so
	make $(TagetDir)/YiProdNtuple

$(TagetDir)/obj/testing.o : $(TagetDir)/src/testing.C
	$(CXX) -c -o $@ $(FLGS) -g $< -MMD
-include $(TagetDir)/obj/testing.d

$(TagetDir)/testing : $(TagetDir)/obj/testing.o
	$(LD) -o $@ $+ $(LIBP)

$(TagetDir)/src/dict.C : $(TagetDir)/src/ClassDef.h $(TagetDir)/src/LinkDef.h
	rootcint -f $@ -c -p $(FLGS) $+
	#rootcling -f $@ -c -p $(FLGS) $+ # only for ROOT6

$(TagetDir)/obj/dict.o : $(TagetDir)/src/dict.C
	$(CXX) -c -o $@ $(FLGS) $<

$(TagetDir)/lib/libClassDef.so : $(TagetDir)/src/dict.C
	$(CXX) -fPIC -shared $^ -o $@ $(FLGS)
	#@cp -fa $(TagetDir)/src/dict_rdict.pcm $(TagetDir)/lib/dict_rdict.pcm # only for ROOT6

$(TagetDir)/obj/YiProdNtuple.o : $(TagetDir)/src/YiProdNtuple.C $(TagetDir)/src/YiProdNtuple.tcc $(TagetDir)/src/YiProdNtuple.h
	$(CXX) -c -o $@ $(FLGS) -g $< -MMD 
-include $(TagetDir)/obj/YiProdNtuple.d

$(TagetDir)/YiProdNtuple : $(TagetDir)/obj/YiProdNtuple.o
	$(LD) -o $@ $+ $(LIBP) -L$(TagetDir)/lib -lClassDef

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/lib/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/testing
	@rm -fv $(TagetDir)/YiProdNtuple
