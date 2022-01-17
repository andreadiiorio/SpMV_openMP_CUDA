CC=gcc
CFLAGS=-Wall -Wextra -ggdb -O2 #-std=c99 #TODO reallocarray
#libs
CFLAGS+=-lm -fopenmp
INCLUDE=-Iinclude
#TODO extra to reduce useless warnings
CFLAGS+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable
CFLAGS+=-Wno-unused-label -Wno-switch #-Wno-unused-function
##TODO OMP CONFIG
#OMP_CANCELLATION=true
#export OMP_CANCELLATION

#SYSTEM CONFIGURATION
UNAME=$(shell uname -a | tr -c -d \[:alnum:\] | tr \[:lower:\] \[:upper:\] ) #upper uname-a
TMPDIR=/run/user/$(shell id -u)/

CONSTS = -DTMPDIR='"$(TMPDIR)"'     #TODO PUT IN /raccolta
MACROS = -DDEBUGPRINT="if(FALSE)" -DCONSISTENCY_CHECKS="if(TRUE)"
MACROS+= -DVERBOSE="if(FALSE)"

UNDEF := $(shell echo $(MACROSDBG) | tr " " "\n" | grep -oe '-D.*=' | tr -d "=" |  sed s/-D/-U/ )


#CUDA
NVCC 			 = /usr/local/cuda/bin/nvcc	#CPATH
INCLUDE_CUDA	 = -I/usr/local/cuda/include/ #TODO SERVERONLY -I/opt/cuda/11.2/samples/common/inc/ 
INCLUDE_CUDA	 = -I/usr/local/cuda/samples/common/inc

WALL_CUDA		:= --Wreorder --Wext-lambda-captures-this # --Wdefault-stream-launch
WALL_CUDA		+= --Wno-deprecated-declarations --Wno-deprecated-gpu-targets 
WALL_CUDA		+= -Xcompiler -Wno-literal-suffix	#TODO not for C warning...
GCC_SET			:= -ccbin=/home/andysnake/DATA/SW/gcc_s/gcc-9.4.0/build/bin/gcc	#/usr/bin/gcc
G++_SET			:= -ccbin=/home/andysnake/DATA/SW/gcc_s/gcc-9.4.0/build/bin/g++	#/usr/bin/gcc
CCAPABILITY		 = -gencode arch=compute_50,code=sm_50 -gencode arch=compute_75,code=sm_75 	#computeCapabilities
OPTS_CUDA		 =  -m64 $(GCC_SET) $(MACROS)	$(CCAPABILITY)	-g -G #export CUDA_DEBUGGER_SOFTWARE_PREEMPTION=1; set cuda api_failures stop
#SUMMARIZE CUDA OPTs
NVCCFLAGS		 = $(WALL_CUDA) $(OPTS_CUDA) -Xcompiler -fopenmp

objs := $(shell  grep -Eo '.*\..*:\s' Makefile | grep -v -e '@' -e PHONY | awk -F: '{print $1}' | tr '\n:' ' ' )
all: $(objs)

SpGEMV_OMP.o: $(filter-out %cuda%, main.cu SpGEMV_*_OMP.c  commons/*.c lib/*.c include/*.h)
	$(CC) -o $@ $(CFLAGS) $(CONSTS) $(MACROS) $(INCLUDE) -x c $(filter-out %.h,$^) $(RUNTIME)


INCLUDE 	   += $(INCLUDE_CUDA)
SpGEMV_CUDA.o: main.cu Sp*.c*  commons/*  include/*.h lib/*.c
	$(NVCC) -o $@ $(NVCCFLAGS) $(CONSTS) $(MACROS) $(INCLUDE) $(filter-out %.h ,$^) $(RUNTIME)

##CUize:
##	find -type l -name "*.cu" | xargs -n 1 rm -f || echo "nothing to deCUize"
##	$(shell find -name "*.c" | xargs -n 1 -I% ln -s -f % %u)
clean:
	rm -i *.o 


#MAKEFILE TESTS FEATURES
headersAll = $(find -name "*.h")
try: $(headersAll) main_OMP.c include/*
	@echo ciao $(filter-out %.c, $^)
	python3 -c 'from os import environ as env;print(env.get("OMP_CANCELLATION"))'

.PHONY: all
