Introduction
============


Why run in parallel?
--------------------

There are two main reasons for using parallel computers:

#. Faster throughput of results
#. Ability to simulate larger systems

These are related to two concepts in parallel computing, `strong` and `weak`
parallel scaling of software. Strong scaling is defined as the how the time to
solve a problem varies with the number of computing processors for a fixed total
problem size. Weak scaling is how the time varies when the problem size
increases with the number of processors.

There can also be restrictions on the amount of available memory, and one
solution is to use a parallel machine where the amount of available memory
increases with the number of computing cores (usually a distributed memory
machine, but shared memory systems often also have a substantial amount of
memory per computing core).


Types of parallel hardware
--------------------------

The types of parallel computer that DFTB+ currently can make use of
are

#. **shared memory** machines, where all processes have access to the same
   data. These are typically small to medium size systems. Most modern CPUs have
   multiple cores, so fall into this category even for desktop machines.

#. **distributed memory** which consists of network connected machines. These are
   typically larger scale systems with higher numbers of processors and a
   dedicated high speed network between them.

#. **GPU accelerated** where either shared memory or distributed
   systems have access to graphical processor units, which provide
   extremely parallel calculations but currently only for ground-state
   calculations.
 
The different system types require distinct program models to make use of the
hardware (however code designed for a distributed memory system can often be
also used for shared memory architectures).


Shared memory parallel
----------------------

The default compilation options for DFTB+ produce an OpenMP enabled parallel
program for shared memory systems. The `make.arch` file for compiling the code
should

#. Include the necessary compiler and linking options for OpenMP. The
   supported `make.arch` examples already do this.

#. A thread parallel LAPACK and BLAS is required and should be specified in
   `make.arch`, along with any extra thread communication libraries. Most
   modern implementations of LAPACK (MKL, openBLAS, ATLAS BLAS, etc.) support
   shared memory parallelism.


Distributed memory parallel
---------------------------

DFTB+ can also be compiled to use MPI parallelism, typically for distributed
memory machines. This is usually required for larger parallel computers and
problem sizes.

This requires additional computational and communication libraries to be
available on your system. The system administrator(s) of your machine may be
able to help locate or configure these for you. The required packages are

* MPI : openMPI and MPICH are common options, but there may be a vendor
  supplied library for your network that has better performance
    
* ScaLAPACK
      
* LAPACK and BLAS : optimised serial implementations 

Sections of the code are currently unable to operate with MPI parallelism
(particularly the excited state calculations), but the majority of the
functionality is the same as the shared memory version.


Hybrid parallelism
------------------

DFTB+ can be compiled with both MPI and openMP parallelism combined. However
using this can require system specific settings for thread affinity to provide
efficiency and this is beyond the scope of this tutorial.

GPU accelerated
---------------

Either the MAGMA (single node) or ELPA (distributed memory) solvers
and libraries are required.
