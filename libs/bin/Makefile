-include $(AMSMKfile)

FLGS += -I$(AMSProjLibs)/CPPLibs
FLGS += -I$(AMSProjLibs)/ROOTLibs
FLGS += -I$(AMSProjLibs)/TRACKLibs

all :
	make TRACK.o
	make libTRACK.so

TRACK.o :  ../TRACKLibs/TRACKLibs.h 
	$(CXX) -c -o $@ $(FLGS) $<

libTRACK.so :  ../TRACKLibs/TRACKLibs.h ../ROOTLibs/ROOTLibs.h ../CPPLibs/CPPLibs.h 
	$(CXX) -fPIC -shared $^ -o $@ $(FLGS)

clean :
