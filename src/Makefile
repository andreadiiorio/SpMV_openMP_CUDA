CC=gcc
CFLAGS=-Wall -Wextra -ggdb -O2 #-std=c99 #TODO reallocarray
#libs
CFLAGS+=-lm -fopenmp

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

CFLAGS+=$(RUNTIME)
objs = SpGEMV_OMP.o #SpGEMM_CUDA.o

all: $(objs)
.PHONY: all

SpGEMV_OMP.o: main_OMP.c SpGEMV_OMP.c  commons/*.c lib/*.c include/*.h
	$(CC) $(CFLAGS) $(CONSTS) $(MACROS) -I include/ $(filter-out %.h,$^) -o $@
clean:
	rm -i *.o 


#MAKEFILE TESTS FEATURES
headersAll = $(find -name "*.h")
try: $(headersAll) main_OMP.c include/*
	@echo ciao $(filter-out %.c, $^)
	python3 -c 'from os import environ as env;print(env.get("OMP_CANCELLATION"))'
