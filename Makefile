objs  = SpGEMV_OMP.o SpGEMV_CUDA.o test.o
srcs  = $(find -name "*.c" )
srcs += $(find -name "*.h" )

all: $(objs)
.PHONY: all

SpGEMV_OMP.o: $(srcs)
	$(MAKE) -C src $@
	ln -fs src/$@ .
test_SpGEMV_OMP.o:	$(srcs)
	$(MAKE) -C test $@
	ln -fs test/$@ .
clean:
	$(MAKE) -C src $@
	rm -i *.o 
