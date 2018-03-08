-include $(AMSMKfile)


#ifeq "$(CXX)" "icc"
#CFLIB=-lifcore
#CF=ifort
#else
#CFLIB=-lgfortran
#CF=gfortran
#endif

CXXFLAGS= -g  $(CPPFLAGS) 

all: libTRDVertex.a 

SRCS= Cell.C Segment.C HitState.C TRD2DHit.C Line.C Vertex2D.C Vertex3D.C TRDVertex.C

MYOBJSS := $(SRCS:$(SRC)%.C=$(BIN)%.o)

SimpleHitState.o: SimpleHitState.h
	@echo Compiling  $< ...	
	$(CXX) -c $(CXXFLAGS)  $< -o $@

$(BIN)%.o: $(SRC)%.C
	@echo Compiling  $< ...	
	$(CXX) -c $(CXXFLAGS)  $< -o $@

libTRDVertex.a: SimpleHitState.o $(MYOBJSS)
	@echo Compilling $< ......
	ar rcs $@ $^

#rootdict.C: HitState.h SimpleHitState.h linkdef.h
#	$(ROOTSYS)/bin/rootcint -f rootdict.C -c -p $(CPPFLAGS) $^

clean:
	rm -f libTRDVertex.a *.o
