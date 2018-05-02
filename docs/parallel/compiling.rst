Compiling the code
==================


Compiling for OpenMP
--------------------

The default `make.arch` examples for DFTB+ should enable OpenMP parallelism, but
in order to gain from this you will also require an efficient thread
parallelised BLAS library. These include `MKL
<https://software.intel.com/en-us/mkl>`_ or `OpenBLAS
<https://www.openblas.net/>`_.

The available `dftb+ binaries <https://www.dftbplus.org/download/dftb-181/>`_
are compiled with OpenMP support and use the OpenBLAS libraries.

You will also usually need to specify the number of parallel threads that the
code can use. This can be done by setting the OMP_NUM_THREADS environment
variable
  
* For the bash shell::
    
    export OMP_NUM_THREADS=<number of threads to use>
     
* For the csh or tcsh shells::
    
    setenv OMP_NUM_THREADS <number of threads to use>

before running DFTB+.


Compiling with MPI
------------------

In order to use Message Passing Interface (MPI) enabled DFTB+, you will require

#. A working MPI 2 (or later) installation (this should include a wrapper for
   Fortran compilation, which is named some variant on mpifort or mpif90 for
   most MPI implementations)

#. The ScaLAPACK library (or a library such as MKL that supplies this
   functionality)

#. Either serial or thread aware LAPACK and BLAS libraries (or equivalents)

Depending on your system, optimised MPI libraries for your network may be
available (contact your system manager for details).


Obtaining the source
^^^^^^^^^^^^^^^^^^^^

If you have downloaded the official 18.1 (or later) release from the `DFTB+
website <http://www.dftb-plus.info/>`_, all essential components for compiling
the code in parallel with MPI support are included (excluding the prerequisites
listed above).

If instead you obtain the code from the `github <https://github.com/dftbplus>`_
repository, the `mpifx <https://github.com/dftbplus/mpifx>`_ and `scalapackfx
<https://github.com/dftbplus/scalapackfx>`_ libraries will be required. These
are included as submodules from the main code, and can be fetched with::

  git submodule update --init --recursive

There is more information on the structure and use of `submodules online
<https://github.com/blog/2104-working-with-submodules>`_.


Building DFTB+ with MPI support enabled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You will need to either edit the file `make.config` to enable::

  WITH_MPI := 1

or compile with ::

  make WITH_MPI=1

to produce an MPI enabled binary. In this case, we suggest serial LAPACK and
BLAS libraries are used (since the main parallelism comes from ScaLAPACK).
