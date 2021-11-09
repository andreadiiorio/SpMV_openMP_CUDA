CC=gcc
CFLAGS=-Wall -Wextra -ggdb -O2 #-std=c99 #TODO reallocarray
#libs
CFLAGS+=-lm -fopenmp

#TODO extra to reduce useless warnings
CFLAGS+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable
CFLAGS+=-Wno-unused-label

MACROS = -DDEBUGPRINT="if(FALSE)" -DCONSISTENCY_CHECKS="if(TRUE)"
MACROS+= -DVERBOSE="if(FALSE)"

objs = SpGEMM_OMP.o SpGEMM_CUDA.o

all: $(objs)
.PHONY: all

SpGEMV_OMP.o: main_OMP.c SpGEMV_OMP.c  commons/*.c lib/*.c include/*.h
	$(CC) $(CFLAGS) $(MACROS) -I include/ $(filter-out %.h,$^) -o $@ # -O0
clean:
	rm -i *.o 


#MAKEFILE TESTS FEATURES
headersAll = $(find -name "*.h")
try: $(headersAll) main_OMP.c include/*
	@echo ciao $(filter-out %.c, $^)
	python3 -c 'from os import environ as env;print(env.get("OMP_CANCELLATION"))'
