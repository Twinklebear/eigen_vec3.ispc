C = gcc
ISPC? = ispc
CFLAGS = -std=c99 -Wall -Wextra -pedantic -O3 -march=native -m64
LDFLAGS = -lm

eigen_vec3.o: eigen_vec3.ispc
	$(ISPC) --arch=x86-64 $^ -o $@ -h eigen_vec3.h

.PHONY: demo
demo: demo.out

demo.out: eigen_vec3.o demo.o
	$(C) $^ -o $@ $(LDFLAGS)

# The error bench program requires Joachim Kopp's hybrid implementation to
# benchmark against. Grab the source www.mpi-hd.mpg.de/personalhomes/globes/3x3/
# and drop it into this directory then build error_bench
.PHONY: error_bench
error_bench: error_bench.out


error_bench.out: eigen_vec3.o error_bench.o dsyevc3.o dsyevh3.o dsyevq3.o dsytrd3.o
	$(C) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm *.o eigen_vec3.h *.out

