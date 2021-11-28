CC=gcc
FF=gfortran
CFLAGS=-Wall -Wextra -O2 
#libs
CFLAGS+=-lm -fopenmp
#TODO extra to reduce useless warnings
CFLAGS+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-switch
CFLAGS+=-Wno-unused-label #-Wno-unused-function
##TODO OMP CONFIG
#OMP_CANCELLATION=true
#export OMP_CANCELLATION

#SYSTEM CONFIGURATION
UNAME=$(shell uname -a | tr -c -d \[:alnum:\] | tr \[:lower:\] \[:upper:\] ) #upper uname-a
TMPDIR=/run/user/$(shell id -u)/

CONSTS = -DTMPDIR='"$(TMPDIR)"'     #TODO PUT IN /raccolta
MACROS = -DDEBUGPRINT="if(FALSE)" -DDEBUG="if(FALSE)" -DCONSISTENCY_CHECKS="if(FALSE)" -DVERBOSE="if(FALSE)"  -DDEBUGCHECKS="if(FALSE)"
MACROSDBG = -DCONSISTENCY_CHECKS="if(TRUE)"  -DDEBUGCHECKS="if(TRUE)" -DVERBOSE="if(TRUE)" -DDEBUG="if(TRUE)"
UNDEF := $(shell echo $(MACROSDBG) | tr " " "\n" | grep -oe '-D.*=' | tr -d "=" |  sed s/-D/-U/ )

CFLAGS+=$(RUNTIME)
#bind to source original project
srcDir=$(realpath ../src)
srcInclDir=$(realpath ../src/include)
coreSrcs:=$(srcDir)/SpGEMV_OMP.c $(srcDir)/commons/*.c  $(srcDir)/lib/*.c $(srcInclDir)/*.h 
testSrcs:=SpGEMV_OMP_test.c  $(coreSrcs)

##CBLAS BUILD NEEDED
CBLAS_INCLUDE=CBLAS_LAPACK/include 
#EXTRA LIBS TO RELINK AFTER FORTRAN COMPILE -- SUPPORT DEBIAN - FEDORA
#EXTRALIBS=/usr/lib64/flexiblas/libflexiblas_openblas-openmp.so
ifeq ($(findstring DEBIAN,$(UNAME)),DEBIAN)
	EXTRALIBS =/usr/lib/x86_64-linux-gnu/libm.so.6
	EXTRALIBS+=/usr/lib/x86_64-linux-gnu/libgomp.so.1
	EXTRALIBS+=/usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0
	EXTRALIBS+=/usr/lib/x86_64-linux-gnu/libpthread.so.0
	EXTRALIBS+=/usr/lib/x86_64-linux-gnu/libc.so.6
	#CFLAGS+= -fPIC #TODO NOT WORKING ON DELL FED BUILD.a
    STATIC_LIBS_DIR=lib_deb
else    #FEDORA
	EXTRALIBS =/usr/lib64/libm.so.6
	EXTRALIBS+=/usr/lib64/libgomp.so.1
	EXTRALIBS+=/usr/lib64/libgomp.so.1.0.0
	EXTRALIBS+=/usr/lib64/libpthread.so.0
	EXTRALIBS+=/usr/lib64/libc.so.6
    STATIC_LIBS_DIR=lib_fed
endif
CBLAS_SLIB=CBLAS_LAPACK/$(STATIC_LIBS_DIR)/libcblas.a
BLAS_SLIB=CBLAS_LAPACK/$(STATIC_LIBS_DIR)/librefblas.a

objs := $(shell  grep -Eo '.*\..*:\s' Makefile | grep -v -e '@' -e PHONY | awk -F: '{print $1}' | tr '\n:' ' ' )
all: $(objs)
#TESTS FOR ALL SPGEMV FUNCS
test_SpGEMV_OMP_Stats.o:	$(testSrcs) 
	$(CC)   $(CFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) -o $@

test_SpGEMV_OMP.o: 	$(testSrcs)
	$(CC)   $(CFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) \
      -DAVG_TIMES_ITERATION=1 $(UNDEF) -DVERBOSE="if(TRUE)" -o $@

#OMP_DISPLAY_ENV=VERBOSE
#export OMP_DISPLAY_ENV
test_SpGEMV_OMP_DBG.o: 	$(testSrcs)
	$(CC)   $(CFLAGS) -O0 -I$(srcDir)/include/  $(filter-out %.h,$^) -ggdb $(CONSTS) $(MACROS) \
	 $(UNDEF)  $(MACROSDBG) -DAVG_TIMES_ITERATION=1  -o $@

test_CBLAS_SpGEMV_OMP.elf: 	$(testSrcs)
	mkdir -p /tmp/oldBuilds
	mv *.o /tmp/oldBuilds
	#create objectfile not linked
	$(CC)  $(CFLAGS) -I$(srcDir)/include/ -I$(CBLAS_INCLUDE) $(filter-out %.h,$^) -c $(CONSTS) $(MACROS)  -DCBLAS_TESTS \
	  -ggdb $(UNDEF) $(MACROSDBG) -DAVG_TIMES_ITERATION=1
	#link it with static compiled cblas and source version library
	$(FF) -O2 -frecursive -o $@ *.o  $(CBLAS_SLIB) $(BLAS_SLIB) $(EXTRALIBS)
	rm -f *.o   #remove newly generated obj file for FF link
	mv /tmp/oldBuilds/* . 

ompGetICV.o:	$(srcDir)/commons/ompGetICV.c $(srcInclDir)/ompGetICV.h
	$(CC) $(CFLAGS) $^ -o $@ -ggdb -D OMP_GET_ICV_MAIN 

testAll:
	$(shell bash ../scripts/testAllMatrixGroups.sh | tee - > log/testAll$GRID_ROWSXGRID_COLS )
	#$(shell cat ../../data/matrixGroups.list | xargs -L 1 sh -c 'echo $0 ${@};./test_CBLAS_SpGEMV_OMP.elf $0 ${@}' ) #TODO FIX

sparseMatrixToImage.o:	sparseMatrixToImage.c sparseMatrixToImage.h $(coreSrcs)
	$(CC) $(CFLAGS) -I$(srcDir)/include -I. $(filter-out %.h,$^) -o $@ $(CONSTS) $(MACROS) \
	    -D MAIN_SPMAT_IMG 
sparseMatrixToImageTest.o:	sparseMatrixToImage.c sparseMatrixToImage.h $(coreSrcs)
	$(CC) $(CFLAGS) -I$(srcDir)/include -I. $(filter-out %.h,$^) -o $@ $(CONSTS) $(MACROS) \
	    -D MAIN_SPMAT_IMG 	-DTEST  -O0 -ggdb
clean:
	rm -i *.elf
	rm -f *.o

try:   
	echo $(TMPDIR) 
	echo $(objs) 
	echo $(UNDEF)
	@echo uname:        $(UNAME)
		
	

.PHONY: all clean testAll
