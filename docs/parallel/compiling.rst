Compiling the code
==================


Compiling for OpenMP
--------------------

The default `make.arch` examples for DFTB+ should enable OpenMP parallelism, but
in order to gain from this you also require an efficient thread parallelised
BLAS library. These include `<MKL https://software.intel.com/en-us/mkl>`_ or
`<OpenBLAS https://www.openblas.net/>`_

The `<available dftb+ https://www.dftbplus.org/download/dftb-181/>`_ binaries
are compiled with OpenMP support and use the OpenBLAS libraries.

You will also usually need to specify the number of parallel threads that the
code can use. This can be done by setting the OMP_NUM_THREADS environment
variable::
  
  * For the bash shell::
    
      export OMP_NUM_THREADS=<number of threads to use>
     
  * For the csh or tcsh shell::
    
      setenv OMP_NUM_THREADS <number of threads to use>

Before running DFTB+

Compiling for MPI
-----------------

In order to use Message Passing Interface (MPI) enabled DFTB+, you will require

1. A working MPI 2 (or later) installation (this will include a wrapper for Fortran
   compilation)
2. The ScaLAPACK library
3. either serial or thread aware LAPACK and BLAS libraries

Depending on your system, optimized MPI libraries for your network may be
available.

Building DFTB+
--------------

You will need to either edit `make.config` to enable::

  WITH_MPI := 1

or compile with ::

  make WITH_MPI=1
