-include $(AMSMKfile)


# Test Ceres Solver
#FLGS += -I/usr/include/eigen3 -I/usr/include/ceres/internal/miniglog
#FLGS += -I/usr/include/eigen3 -I/usr/local/include
#LIBP += -L/usr/lib64 -lceres -lgflags
#LIBP +=  -L/usr/local/lib -lglog -L/usr/lib64 -lceres -lgflags
#ICCCXX += -D__CeresSolver__
#GCCCXX += -D__CeresSolver__
#LD += -D__CeresSolver__
#LDS += -D__CeresSolver__

TagetDir = $(AMSCoreSubj)/appStudy/vdev

Version=18Feb13
ClassLib = $(AMSCoreProd)/${Version}/lib/libClassDef.so

all :
	@mkdir -p $(TagetDir)
	@mkdir -p $(TagetDir)/env
	@mkdir -p $(TagetDir)/obj
	@mkdir -p $(TagetDir)/src
	@mkdir -p $(TagetDir)/lst
	@cp -fa $(PWD)/lst/* $(TagetDir)/lst 2> /dev/null
	@cp -fa $(AMSENV) $(TagetDir)/env/amsenv.sh 2> /dev/null
	@cp -fa $(PWD)/src/app_fill.C $(TagetDir)/src/app_fill.C 2> /dev/null
	make $(TagetDir)/app_fill

$(TagetDir)/obj/app_fill.o : $(TagetDir)/src/app_fill.C
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
-include $(TagetDir)/obj/app_fill.d

$(TagetDir)/app_fill : $(TagetDir)/obj/app_fill.o
	$(LD) -o $@ $+ $(LIBP) $(ClassLib)

#$(TagetDir)/obj/app_fit.o : $(TagetDir)/src/app_fit.C
#	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $< -MMD
#-include $(TagetDir)/obj/app_fit.d
#
#$(TagetDir)/app_fit : $(TagetDir)/obj/app_fit.o
#	$(LD) -o $@ $+ $(LIBP)

clean :
	@rm -fv $(TagetDir)/env/*
	@rm -fv $(TagetDir)/obj/*
	@rm -fv $(TagetDir)/src/*
	@rm -fv $(TagetDir)/lst/*
	@rm -fv $(TagetDir)/app_fill
