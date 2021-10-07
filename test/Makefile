CC=gcc
FF=gfortran
CFLAGS=-Wall -Wextra -ggdb -O2
#libs
CFLAGS+=-lm
#TODO extra to reduce useless warnings
CFLAGS+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable

#bind to source original project
srcDir=../src
##CBLAS BUILD NEEDED
CBLAS_INCLUDE=CBLAS_LAPACK/include 
CBLAS_SLIB=CBLAS_LAPACK/lib/libcblas.a
BLAS_SLIB=CBLAS_LAPACK/lib/librefblas.a
#EXTRA LIBS
#EXTRALIBS=/usr/lib64/flexiblas/libflexiblas_openblas-openmp.so
EXTRALIBS =/usr/lib64/libm.so.6
EXTRALIBS+=/usr/lib64/libgomp.so.1
EXTRALIBS+=/usr/lib64/libgomp.so.1.0.0
EXTRALIBS+=/usr/lib64/libpthread.so.0
EXTRALIBS+=/usr/lib64/libc.so.6

test_SpGEMV_OMP: SpGEMV_OMP_test.c  $(srcDir)/commons/utils.c  $(srcDir)/lib/parser.c  $(srcDir)/lib/mmio.c
	#create objectfile not linked
	$(CC) -fopenmp $(CFLAGS) -I$(srcDir)/include/ -I$(CBLAS_INCLUDE) $^ -c 
	#link it with static compiled cblas and source version library
	$(FF) -O2 -frecursive -o $@.elf *.o  $(CBLAS_SLIB) $(BLAS_SLIB) $(EXTRALIBS)

clean:
	rm -i *.o 

objs = test_SpGEMV_OMP.o
all: $(objs)
.PHONY: all
