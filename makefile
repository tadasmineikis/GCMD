CC    		?= g++
CFLAGS		= -c -Wall $(INCLUDE) -std=gnu++11
DEBUG 		= -g3 -ggdb3
OPTIM 		= -O2 -funroll-loops -fvariable-expansion-in-unroller -fopenmp -march=native
USE_FLAGS       = -Wall $(OPTIM)
OBJS  		= main_all_bands.o iso.o imf.o lib.o configFile.o clust.o imf_pflamm.o
EXE   		= gCMD_0.21.5
INCLUDE         = 
RUN_DIR		= 
LIB     	= -larmadillo -lgsl -lgomp -lgslcblas


all: gcmd

gen: generator move_gen

gcmd: $(OBJS)
	$(CC)  $(USE_FLAGS) $(LFLAGS) -o $(EXE) $(OBJS) $(LIB)

main_all_bands.o: iso.h imf.h main_all_bands.cpp lib.h configFile.h clust.h
	$(CC) $(CFLAGS) $(USE_FLAGS) main_all_bands.cpp

iso.o: iso.cpp iso.h lib.h
	$(CC) $(CFLAGS) $(USE_FLAGS) iso.cpp

lib.o:	lib.cpp lib.h
	$(CC) $(CFLAGS) $(USE_FLAGS) lib.cpp

imf.o:	imf.h lib.h imf.cpp
	$(CC) $(CFLAGS) $(USE_FLAGS) imf.cpp

configFile.o: configFile.h lib.h configFile.cpp filters.h
	$(CC) $(CFLAGS) $(USE_FLAGS) configFile.cpp

fileio.o: fileio.h lib.h fileio.cpp
	$(CC) $(CFLAGS) $(USE_FLAGS) fileio.cpp

clust.o: clust.h imf.h configFile.h iso.h filters.h clust.cpp
	$(CC) $(CFLAGS) $(USE_FLAGS) clust.cpp

imf_pflamm.o: imf_pflamm.h imf_pflamm.c
	$(CC) $(CFLAGS) $(USE_FLAGS) imf_pflamm.c
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


