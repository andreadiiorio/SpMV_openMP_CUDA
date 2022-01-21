CC=gcc
FF=gfortran

CWALL=-Wall -Wextra -O2 
CLIBS=-lm -fopenmp
#TODO extra to reduce useless warnings
CWALL+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-switch
CWALL+=-Wno-unused-label #-Wno-unused-function
CFLAGS=$(CWALL) $(CLIBS) $(RUNTIME) -x c

#SYSTEM CONFIGURATION
UNAME=$(shell uname -a | tr -c -d \[:alnum:\] | tr \[:lower:\] \[:upper:\] ) #upper uname-a
TMPDIR=/run/user/$(shell id -u)/
#MACROS
CONSTS = -DTMPDIR='"$(TMPDIR)"'  -DAVG_TIMES_ITERATION=25   #TODO PUT IN /raccolta
MACROS = -DDEBUGPRINT="if(FALSE)" -DDEBUG="if(FALSE)" -DCONSISTENCY_CHECKS="if(FALSE)" -DVERBOSE="if(FALSE)"  -DDEBUGCHECKS="if(FALSE)" 
MACROSDBG = -DCONSISTENCY_CHECKS="if(TRUE)"  -DDEBUGCHECKS="if(TRUE)" -DVERBOSE="if(TRUE)" -DDEBUG="if(TRUE)"
MACROSDBG += -D DOUBLE_VECT_DIFF_EARLY_EXIT
UNDEF := $(shell echo $(MACROSDBG) | tr " " "\n" | grep -oe '-D.*=' | tr -d "=" |  sed s/-D/-U/ )
CFLAGSDBG = -ggdb #-Og #-O0	#TODO -O0 cause REDUCTION MISWRITE...FAULTY RESOULT (OR SERIAL ERROR)
#GATHER ALL DEBUG STUFF IN A VAR
DBG=$(CFLAGSDBG) $(UNDEF) -UAVG_TIMES_ITERATION -DAVG_TIMES_ITERATION=1 $(MACROSDBG)  

#CUDA NVCC
NVCC 			 = /usr/local/cuda/bin/nvcc	#CPATH
INCLUDE_CUDA	 = -I/usr/local/cuda/include/ 
INCLUDE_CUDA	+= -I/opt/cuda/11.2/samples/common/inc/ 
INCLUDE_CUDA	+= -I/usr/local/cuda/samples/common/inc
WALL_CUDA		:= --Wreorder --Wext-lambda-captures-this # --Wdefault-stream-launch
WALL_CUDA		+= --Wno-deprecated-declarations --Wno-deprecated-gpu-targets 
WALL_CUDA		+= -Xcompiler -Wno-literal-suffix	#TODO not for C files f*!* warning
#custom host compiler binDir set, TODO gcc11 not yet supported :(
GCC_SET			:= -ccbin=/home/andysnake/DATA/SW/gcc_s/gcc-9.4.0/build/bin/gcc	#/usr/bin/gcc
G++_SET			:= -ccbin=/home/andysnake/DATA/SW/gcc_s/gcc-9.4.0/build/bin/g++	#/usr/bin/gcc
CCAPABILITY		 = -gencode arch=compute_50,code=sm_50 		#NV960M 			DELL
#CCAPABILITY	 = -gencode arch=compute_75,code=sm_75 	#Quadro RTX 5000	SERVER
OPTS_CUDA		 =  -m64 $(GCC_SET) $(MACROS)	$(CCAPABILITY)
#SUMMARIZE CUDA OPTs
NVCCFLAGS		 = $(WALL_CUDA) $(OPTS_CUDA) -Xcompiler -fopenmp $(INCLUDE_CUDA)
NVCCFLAGSDBG	 = -g -G #export CUDA_DEBUGGER_SOFTWARE_PREEMPTION=1; set cuda api_failures stop
DBG_CUDA=$(NVCCFLAGSDBG) $(UNDEF) -UAVG_TIMES_ITERATION -DAVG_TIMES_ITERATION=1 $(MACROSDBG)  


#bind to source original project
srcDir		=$(realpath ../src)
srcInclDir	=$(realpath ../src/include)
coreSrcs	:=$(srcDir)/commons/*.c  $(srcDir)/lib/*.c $(srcInclDir)/*.h 
ompSrcs		:=$(coreSrcs) $(srcDir)/SpMV*OMP.c
cudaSrcs	:=$(coreSrcs) $(srcDir)/commons/*.cu $(srcDir)/Sp*.cu
testSrcsOmp	:=SpMV_test.cu  $(ompSrcs)
testSrcsCuda:=SpMV_test.cu  $(cudaSrcs)
testSrcsAll	:=SpMV_test.cu  $(ompSrcs) $(cudaSrcs)

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
 
#objs := $(shell  grep -Eo '.*\..*:\s' Makefile | grep -v -e '@' -e PHONY | awk -F: '{print $1}' | tr '\n:' ' ' ) #also .elf
objs := $(shell  grep -Eo '.*\.o:\s' Makefile | grep -v -e '@' -e PHONY | awk -F: '{print $1}' | tr '\n:' ' ' )
all: $(objs)
##FAST TESTS - PERFGATHER
test_SpMV_OMP_Stats.o:	$(testSrcsOmp) 
	$(CC) -o $@  $(CFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) 

test_SpMV_OMP_Stats_NoReductionSIMD.o:	$(testSrcsOmp) 
	$(CC)  -o $@ $(CFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS)  -DSIMD_ROWS_REDUCTION=FALSE

test_SpMV_OMP.o: 	$(testSrcsOmp)
	$(CC) -o $@  $(CFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) \
      -UAVG_TIMES_ITERATION  -DAVG_TIMES_ITERATION=1 $(UNDEF) -DVERBOSE="if(TRUE)" 
##CUDA VERSIONS
test_SpMV_CUDA_Stats.o:	$(testSrcsAll) 
	$(NVCC) -o $@  $(NVCCFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) 

test_SpMV_CUDA.o: 	$(testSrcsAll)  
	$(NVCC) -o $@  $(NVCCFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) \
      -UAVG_TIMES_ITERATION  -DAVG_TIMES_ITERATION=1 $(UNDEF) -DVERBOSE="if(TRUE)" 
test_SpMV_CUDA_DBG.o: 	$(testSrcsAll) 
	$(NVCC) -o $@  $(NVCCFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) \
      $(UNDEF) -DVERBOSE="if(TRUE)" $(DBG_CUDA)
##DEBUG BUILDS
test_SpMV_OMP_DBG.o: 	$(testSrcsOmp)
	$(CC) -o $@  $(CFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) $(DBG) 

test_SpMV_OMP_DBG_NoReductionSIMD.o: 	$(testSrcsOmp)
	$(CC) -o $@  $(CFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) \
	 $(DBG)  -DSIMD_ROWS_REDUCTION=FALSE

test_CBLAS_SpMV_OMP.elf: 	$(testSrcsOmp) SpMV_CBLAS.c
	mkdir -p /tmp/oldBuilds
	mv *.o /tmp/oldBuilds
	#create objectfile not linked
	$(CC) -c $(CFLAGS) -I$(srcDir)/include/ -I$(CBLAS_INCLUDE) $(filter-out %.h,$^) \
	  $(CONSTS) $(MACROS)  -DCBLAS_TESTS $(DBG)
	#link it with static compiled cblas and source version library
	$(FF) -O2 -frecursive -o $@ *.o  $(CBLAS_SLIB) $(BLAS_SLIB) $(EXTRALIBS)
	rm -f *.o   #remove newly generated obj file for FF link
	mv /tmp/oldBuilds/* . 

test_utils.o: 	$(srcDir)/commons/*.c $(srcDir)/lib/*.c $(srcInclDir)/*.h 
	$(CC) -o $@  $(CFLAGS) -I$(srcDir)/include/  $(filter-out %.h,$^) $(CONSTS) $(MACROS) $(DBG) \
	-DSPARSEUTILS_MAIN_TEST -O0

ompGetICV.o:	$(srcDir)/commons/ompGetICV.c $(srcInclDir)/ompGetICV.h
	$(CC) -o $@ $(CFLAGS) $(filter-out %.h,$^) -ggdb -D OMP_GET_ICV_MAIN 

testAll:
	$(shell bash ../scripts/testAllMatrixGroups.sh | tee - > log/testAll$GRID_ROWSXGRID_COLS )
	#$(shell cat ../../data/matrixGroups.list | xargs -L 1 sh -c 'echo $0 ${@};./test_CBLAS_SpMV_OMP.elf $0 ${@}' ) #TODO FIX

sparseMatrixToImage.o:	sparseMatrixToImage.c sparseMatrixToImage.h $(coreSrcs)
	$(CC) -o $@ $(CFLAGS) -I$(srcDir)/include -I. $(filter-out %.h,$^) $(CONSTS) $(MACROS) \
	    -D MAIN_SPMAT_IMG 
sparseMatrixToImageTest.o:	sparseMatrixToImage.c sparseMatrixToImage.h $(coreSrcs)
	$(CC) -o $@ $(CFLAGS) -I$(srcDir)/include -I. $(filter-out %.h,$^) $(CONSTS) $(MACROS) \
	    -D MAIN_SPMAT_IMG 	-DTEST  $(CFLAGSDBG)
clean:
	rm -f *.elf
	rm -f *.o

try:   
	echo $(TMPDIR) 
	echo $(objs) 
	echo $(UNDEF)
	@echo uname:        $(UNAME)
		
	

.PHONY: all clean testAll
