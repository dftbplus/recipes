************
Introduction
************

This document is a collection of examples demonstrating the usage of the quantum
mechanical atomistic software package `DFTB+ <http://www.dftbplus.org>`_.

The examples assume that you have the `latest stable version
<http://www.dftbplus.org/download/dftb-stable/>`_ of DFTB+ installed (although
many of them may run also with older versions).  Additionally they require some
parameterisation files (Slater-Koster files), which can be downloaded
from `dftb.org <http://www.dftb.org>`_.

The recipes often only show the relevant part of the input. In order to obtain
the full input and in order to run the examples, you should download the archive
containing all the inputs of the individual recipes.

.. only :: builder_html or readthedocs

   See :download:`archive with all inputs <_downloads/archives/recipes.tar.bz2>`.

In the various recipes, the directory, where the full self-containing input can
be found within this archive, is indicated in square brackets right after the
section title (e.g. [Input: `recipes/basics/firstcalc/`]).

After unpacking the archive, you can also automatically download all
Slater-Koster-files needed by the examples by issuing ::

  ./scripts/get_slakos

in the root folder of the archive. Then you should be able to run each example
by calling ``dftb+`` from the corresponding input folder.
