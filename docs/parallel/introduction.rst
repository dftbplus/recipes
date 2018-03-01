============
Introduction
============

Why run in parallel?
--------------------

There are two main reasons for using parallel computers::

  1. Faster throughput of results
  2. Ability to simulate larger systems

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

Shared memory parallel
----------------------

The default compilation options for DFTB+ produce an OpenMP enabled parallel
program. The `make.arch` file for compiling the code should ::

  1. Include the neccessary compiler and linking options for openMP. The
     supported `make.arch` examples already do this.
  2. A thread parallel LAPACK and BLAS is required and should be specified in
     `make.arch`, along with any extra thread communication libraries. Most
     modern implementations of LAPACK (MKL, openBLAS, ATLAS BLAS, etc.) support
     shared memory parallelism.

Distributed memory parallel
---------------------------

DFTB+ can also be compiled to use MPI parallelism. This is usually required for
larger parallel computers and problem sizes.

This requires additional computational and communication libraries to be
available on your system. The system administrator(s) of your machine may be
able to help locate or configure these for you. The required packages are::

  * MPI : openMPI and MPICH are common options, but there may be a vendor
    supplied library for your network that has better performance
    
  * ScaLAPACK
      
  * LAPACK and BLAS : optimised serial implementations 

Sections of the code are currently unable to operate with MPI parallelism
(particularlty the excited state calculations).

Hybrid parallelism
------------------

DFTB+ can be compiled with both MPI and openMP parallelism combined. However
using this can require thread affinity for efficency and this is beyond the
scope of this tutorial.
