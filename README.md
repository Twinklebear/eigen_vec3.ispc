eigen\_vec3.ispc
---
This is a port of Joachim Kopp's excellent hybrid 3x3 real-symmetric eigenvector and eigenvalue solver
over to ispc. This version uses [ispc](https://ispc.github.io/) to solve systems in parallel on the CPU's
SIMD lanes and to easily extend to multiple threads.

All credit for the implementation should go to the original author, Kopp, who introduced the hybrid method in:

**Efficient numerical diagonalization of hermitian 3x3 matrices**
Int. J. Mod. Phys. C 19 (2008) 523-548
arXiv.org: [physics/0610206](http://arxiv.org/abs/physics/0610206)

The source for the various methods discussed in the paper is provided [online](www.mpi-hd.mpg.de/personalhomes/globes/3x3/)
and is LGPGL v2 licensed.

If you require a fast double precision implementation Kopp's C or FORTRAN implementation is likely a better
choice, at the time of writing double-precision doesn't perform very well in ispc. However if you
don't need perfect accuracy on some (admittedly non-rigourous) benchmarks I found that the single
precision ispc version can compute eigenvectors and eigenvalues within 0.00005 of Kopp's double precision C version
about 99.5% of the time (see [error\_bench.c](error_bench.c)). To run some performance demos
and see an example of usage from C see [demo.c](demo.c).

