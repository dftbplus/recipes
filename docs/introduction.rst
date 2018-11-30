.. _sec-introduction:

************
Introduction
************

This document is a collection of examples demonstrating the usage of the quantum
mechanical atomistic software package `DFTB+ <http://www.dftbplus.org>`_.

Before you start
================

The examples assume that you have the `latest stable version
<http://www.dftbplus.org/download/dftb-stable/>`_ of DFTB+ installed (although
many of them may run also with older versions).  Additionally they require some
parameterisation files (Slater-Koster files), which can be downloaded
from `dftb.org <http://www.dftb.org>`_.

The recipes often only show the relevant parts of the input. In order to obtain
the full input and in order to run the examples yourself, download the archive
containing all the inputs of the individual recipes.

.. only :: builder_html or readthedocs

   Download :download:`archive with all inputs
   <_downloads/archives/recipes.tar.bz2>`.

.. only :: not (builder_html or readthedocs)

   The archive can be downloaded from the `online version of the DFTB+ recipes
   <https://dftbplus-recipes.readthedocs.io/>`_ at
   https://dftbplus-recipes.readthedocs.io/.
   
In each recipe we indicate in square brackets after the section titles the
corresponding directories in the archive, where the self-contained input can be
found (e.g. [Input: `recipes/basics/firstcalc/`]).

After unpacking the archive, you should download all Slater-Koster-files needed
by the examples by issuing ::

  ./scripts/get_slakos

in the root folder of the archive. After this you should be able to run each
example. The `run.sh` scripts in the example folders contain the individual
commands needed to run the full example.


Where to start
==============

The individual chapters are more or less independent from each other, so you may
directly go to the one closest to your interests. However, if you are new to
DFTB+, please make sure to work through the relevant chapters in
:ref:`sec-basics` first.

The recipes serve the puprose of enabling you to start with a given
functionality of DFTB+ and are, therefore, rather short and focused. Please
always consult also the corresponding sections of the `DFTB+ manual
<http://www.dftbplus.org/fileadmin/DFTBPLUS/public/dftbplus/latest/manual.pdf>`_
for further details and possibilities.

Please note that the example outputs in the recipes may have been created with
older versions of DFTB+ and may, therefore, differ slightly in their format from
the ones you obtain. The corresponding inputs in the archive should work with
the last stable release of DFTB+ without any changes, though.
