.. _sec-interfaces:

###########################
Interfaces with other codes
###########################

In this section some of the methods for communication between DFTB+ and external
software are discussed. This can be used to improve ease of use and allows you
to invoke DFTB+ via software that you may be more more familiar with, such as
the Atomic Simulation Environment (`ASE <https://wiki.fysik.dtu.dk/ase/>`_) or
to access DFTB+ functionality via your own code in languages like Python.

Use of interface communication also offers the possibility of expanding the
applications and functionality of DFTB+. For example, DFTB+ can serve as an
energy/force engine for an external driver which then could perform calculations
like geometry optimisation, coupled quantum / molecular mechanics calculations
or more advanced molecular dynamics methods at the DFTB level.

.. toctree::
   :maxdepth: 1

   ase/index.rst
   pyapi/pyapi.rst
