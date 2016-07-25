CXX    		?= g++
CFLAGS		= -c -Wall $(INCLUDE) -std=gnu++11
DEBUG 		= -g3 -ggdb3
OPTIM 		= -O2 -funroll-loops -fvariable-expansion-in-unroller -fopenmp -march=native
USE_FLAGS   = -Wall $(OPTIM)
OBJS_UBV	= main_all_bands.o iso.o imf.o lib.o configFile_ubv.o clust_ubv.o imf_pflamm.o
OBJS_ACS	= main_all_bands.o iso.o imf.o lib.o configFile_acs.o clust_acs.o imf_pflamm.o
EXE_UBV		= gCMD_0.21.5_ubv
EXE_ACS		= gCMD_0.21.5_acs
INCLUDE     = 
RUN_DIR		= 
LIB     	= -larmadillo -lgsl -lgomp -lgslcblas


all: gcmd_ubv gcmd_acs

gen: generator move_gen

gcmd_ubv: $(OBJS_UBV)
	$(CXX)  $(USE_FLAGS) $(LFLAGS) -o $(EXE_UBV) $(OBJS_UBV) $(LIB)

gcmd_acs: $(OBJS_ACS)
	$(CXX)  $(USE_FLAGS) $(LFLAGS) -o $(EXE_ACS) $(OBJS_ACS) $(LIB)

main_all_bands.o: iso.h imf.h main_all_bands.cpp lib.h configFile.h clust.h
	$(CXX) $(CFLAGS) $(USE_FLAGS) main_all_bands.cpp

iso.o: iso.cpp iso.h lib.h
	$(CXX) $(CFLAGS) $(USE_FLAGS) iso.cpp

lib.o:	lib.cpp lib.h
	$(CXX) $(CFLAGS) $(USE_FLAGS) lib.cpp

imf.o:	imf.h lib.h imf.cpp
	$(CXX) $(CFLAGS) $(USE_FLAGS) imf.cpp

fileio.o: fileio.h lib.h fileio.cpp
	$(CXX) $(CFLAGS) $(USE_FLAGS) fileio.cpp

configFile_ubv.o: configFile.h lib.h configFile.cpp filters_UBV.h
	$(CXX) $(CFLAGS) $(USE_FLAGS) -DFILTERS_UBV configFile.cpp -o configFile_ubv.o

configFile_acs.o: configFile.h lib.h configFile.cpp filters_ACS.h
	$(CXX) $(CFLAGS) $(USE_FLAGS) -DFILTERS_ACS configFile.cpp -o configFile_acs.o

clust_ubv.o: clust.h imf.h configFile.h iso.h filters_UBV.h clust.cpp
	$(CXX) $(CFLAGS) $(USE_FLAGS) -DFILTERS_UBV clust.cpp -o clust_ubv.o

clust_acs.o: clust.h imf.h configFile.h iso.h filters_ACS.h clust.cpp
	$(CXX) $(CFLAGS) $(USE_FLAGS) -DFILTERS_ACS clust.cpp -o clust_acs.o

imf_pflamm.o: imf_pflamm.h imf_pflamm.c
	$(CXX) $(CFLAGS) $(USE_FLAGS) imf_pflamm.c
#########################################
# utilities
#########################################

clean:
	rm -f *.o *~ $(EXE)

tar:
	tar cfv $(EXE).tar *.h *.cpp *.c changelog makefile

zip:
	gzip -f $(EXE).tar

arch: tar zip

move:
	ln -sf $(CURDIR)/$(EXE) $(RUN_DIR)$(EXE)


