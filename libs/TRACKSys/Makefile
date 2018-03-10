-include $(AMSMKfile)

all: lib64/libTRACKSys.a

lib64/TRACKSys.o: include/TRACKSysLibs.h
	@echo Compiling  $< ...	
	$(CXX) -c $< -o $@

lib64/libTRACKSys.a: lib64/TRACKSys.o
	@echo Compilling $< ......
	ar rcs $@ $^

clean:
	rm -f lib64/libTRACKSys.a lib64/TRACKSys.o
