objs  = SpGEMV_OMP.o SpGEMV_CUDA.o test.o
srcs  = $(find -name "*.c" )
srcs += $(find -name "*.h" )

all: $(objs)
.PHONY: all

SpGEMV_OMP.o: $(srcs)
	$(MAKE) -C src $@
	ln -fs src/$@ .
clean:
	$(MAKE) -C src $@
	rm -i *.o 
