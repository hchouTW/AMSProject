-include $(AMSMKfile)

all: lib64/libTRACKSys.a

lib64/%s.o: include/%.C
	@echo Compiling  $< ...	
	$(CXX) -c $< -o $@

lib64/TRACKSys.o: include/TRACKSys.h
	@echo Compiling  $< ...	
	$(CXX) -c $< -o $@

lib64/libTRACKSys.a: lib64/TRACKSys.o
	@echo Compilling $< ......
	ar rcs $@ $^

clean:
	rm -f lib64/libTRACKSys.a lib64/*.o
