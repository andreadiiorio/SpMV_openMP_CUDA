CC=gcc
CFLAGS=-Wall -Wextra -ggdb -O2
#libs
CFLAGS+=-lm

#TODO extra to reduce useless warnings
CFLAGS+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable

objs = SpGEMV_OMP.o SpGEMV_CUDA.o test.o

all: $(objs)
.PHONY: all

SpGEMV_OMP.o: SpGEMV_OMP.c  commons/utils.c lib/parser.c lib/mmio.c
	$(CC) -fopenmp $(CFLAGS) -I include/ $^ -o $@
clean:
	rm -i *.o 
