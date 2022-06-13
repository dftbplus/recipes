.. _sec-mpi:

MPI parallel calculations
=========================

Architecture
------------

MPI parallel calculations can be run on shared memory machines, but
for larger calculations, distributed memory systems are used. These
consist of separate computing nodes (often with multiple processors)
networked together, typically with a shared file and queuing system.

Processor groups
----------------

* Use of the `Groups` keyword in DFTB+ to improve scaling and efficiency for
  calculations with spin polarisation and/or k-points.


Solvers
-------

ScaLAPACK solvers
^^^^^^^^^^^^^^^^^

The `QR`, `DivideAndConquer` and `RelativelyRobust` dense ScaLAPACK
eigensolvers are available for the serial and MPI parallel
DFTB+. However, if possible, we suggest using the ELSI solvers for
their better performance and capabilities.

ELSI solvers
^^^^^^^^^^^^

Some of the [ELSI]_ solvers which can be used are

  +--------+-------------------+--------------------+---------------+
  | Solver | Use case          | Scaling with atoms | Eigenvalues   |
  +========+===================+====================+===============+
  | ELPA   | General (dense)   | O(N :sup:`3` )     | available     |
  +--------+-------------------+--------------------+---------------+
  | PEXSI  | General (sparse)  | between            | not available |
  |        |                   | O(N:sup:`2`) and   |               |
  |        |                   | O(N)               |               |
  +--------+-------------------+--------------------+---------------+
  | NTPoly | Better for gapped | O(N)               | not available |
  |        | systems (sparse)  |                    |               |
  +--------+-------------------+--------------------+---------------+
  

The PEXSI solver is not included in the condaforge binaries, but can
be built separately.


Examples
--------



Note about transport calculations
---------------------------------

For :ref:`electronic transport <sec-transport>`, MPI parallelism is also supported.


References
==========

.. [ELSI] ELSI -- An open infrastructure for electronic structure
           solvers, Yu *et al.* (2020) DOI: `j.cpc.2020.107459
           <https://doi.org/10.1016/j.cpc.2020.107459>`_
