************
Introduction
************

This document is a collection of examples demonstrating the usage of the quantum
mechanical atomistic software package `DFTB+ <http://www.dftbplus.org>`_.

The examples assume that you have the `latest stable version
<http://www.dftbplus.org/download/dftb-stable/>`_ of DFTB+ installed (although
many of the examples may run also with older versions).  Additionally the
examples require some parameterisation files (Slater-Koster files), which can be
downloaded from `dftb.org <http://www.dftb.org>`_ (see below).

The description of the recipes often only shows the relevant part of the
input. In order to obtain the full input files and to run the examples, you
should download the archive file containing all the inputs of the individual
recipes.

.. only :: builder_html or readthedocs

   See :download:`archive with all inputs <_downloads/archives/recipes.tar.bz2>`.

Once unpacked, the full self-contained input files can be found at locations
given in square brackets after the section title of each example (e.g. [Input:
`recipes/basics/firstcalc/`]) in the root folder of the archive.


Getting Slater-Koster data
~~~~~~~~~~~~~~~~~~~~~~~~~~

After unpacking the archive you can also automatically download all the
Slater-Koster-files needed by the examples by issuing the command::

  ./scripts/get_slakos

in the root folder of the archive. You should then be able to run each example
by calling ``dftb+`` from the corresponding input folder (assuming that dftb+ is
installed in your executable path).

For some examples, there is also a supplied script file in the directory to run
multistage calculations, this is called ``run.sh``.
