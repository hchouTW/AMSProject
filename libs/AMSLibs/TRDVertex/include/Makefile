-include $(AMSMKfile)

all: libTRDVertex.a 

SRCS = Cell.C Segment.C HitState.C TRD2DHit.C Line.C Vertex2D.C Vertex3D.C TRDVertex.C

MYOBJSS := $(SRCS:$(SRC)%.C=$(BIN)%.o)

SimpleHitState.o: SimpleHitState.h
	@echo Compiling  $< ...	
	$(CXX) -c $< -o $@

$(BIN)%.o: $(SRC)%.C
	@echo Compiling  $< ...	
	$(CXX) -c $< -o $@

libTRDVertex.a: SimpleHitState.o $(MYOBJSS)
	@echo Compilling $< ......
	ar rcs $@ $^

clean:
	rm -f libTRDVertex.a *.o
