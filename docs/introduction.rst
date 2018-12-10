.. _sec-introduction:

************
Introduction
************

This document is a collection of examples demonstrating the usage of the
atomistic quantum mechanical software `DFTB+ <http://www.dftbplus.org>`_.

Before you start
================

The examples assume that you have the `latest stable version
<http://www.dftbplus.org/download/dftb-stable/>`_ of DFTB+ installed (although
many of the recipes may also work with older versions of the code).
Additionally the examples require some parameterisation files (Slater-Koster
files), which can be downloaded from `dftb.org <http://www.dftb.org>`_.

The recipes in this document often only show the relevant parts of the input. In
order to obtain the full input files and in order to run the examples yourself,
please download the archive containing all the inputs of the individual recipes.

.. only :: builder_html or readthedocs

   Download :download:`archive with all inputs
   <_downloads/archives/recipes.tar.bz2>`.

.. only :: not (builder_html or readthedocs)

   This can be downloaded from the `online version of the DFTB+ recipes
   <https://dftbplus-recipes.readthedocs.io/>`_ at
   https://dftbplus-recipes.readthedocs.io/.
   
In each recipe we indicate where to find the corresponding directories in the
archive with square brackets after the section title (e.g. [Input:
`recipes/basics/firstcalc/`]).

Getting Slater-Koster data
~~~~~~~~~~~~~~~~~~~~~~~~~~

After unpacking the archive, you must also download all Slater-Koster-files
needed by the examples, by using the supplied script as ::

  ./scripts/get_slakos

in the root folder of the archive. You should then be able to run each example
by calling ``dftb+`` from the corresponding input folder (assuming that dftb+ is
installed in your executable path).  For some examples, there is also a supplied
script file in the directory to run examples of multistage calculations, called
``run.sh`` which contains the individual commands needed to run the full
example.



Where to start
==============

The individual chapters are more or less independent from each other, so you may
go directly to the one relevant to your interests. However, if you are new to
DFTB+, please make sure to work through the relevant introductory examples in
the :ref:`sec-basics` chapters first.

The recipes are to introduce you to specific functionalities of DFTB+ and so
are, therefore rather short and focused. Please also always consult the
corresponding sections of the `DFTB+ manual
<http://www.dftbplus.org/fileadmin/DFTBPLUS/public/dftbplus/latest/manual.pdf>`_
for further details and possibilities.

Please note that the example outputs in the recipes may have been created with
older versions of DFTB+ and therefore could differ slightly in format from
output of the most recent code. The corresponding inputs in the archive should
work, without any changes, with the last stable release of DFTB+.
