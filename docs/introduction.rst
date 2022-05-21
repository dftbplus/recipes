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

After unpacking the archive, you must also download all
Slater-Koster-files needed by the examples in this tutorial, by using
the supplied script as ::

  ./scripts/get_slakos

in the root folder of the archive. You should then be able to run each example
by calling ``dftb+`` from the corresponding input folder (assuming that dftb+ is
installed in your executable path).  For some examples, there is also a supplied
script file in the directory to run examples of multistage calculations, called
``run.sh`` which contains the individual commands needed to run the full
example.



Installing DFTB+ from conda-forge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DFTB+ is available with the cross platform package manager `conda
<https://en.wikipedia.org/wiki/Conda_(package_manager)>`_ on the
`conda-forge <https://conda-forge.org>`_ channel.

If you have no conda installation yet, we recommend you bootstrap an
installation with the conda-forge distribution `miniforge
<https://github.com/conda-forge/miniforge/releases/latest>`_ or the
anaconda distribution `miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_. If you are
unfamiliar with conda, their `getting started guide
<https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_
introduces the environment.

Installing ``dftb+`` from the conda-forge channel can then be achieved
by adding conda-forge to your channels with:

.. code-block:: none

   conda config --add channels conda-forge

We suggest using the `mamba installer
<https://mamba.readthedocs.io/>`_ with conda-forge, as we have
experienced dependency resolution problems with the original conda
installer::

  conda install -n base mamba

Once the conda-forge channel has been enabled, ``dftb+`` can be
installed with:

.. code-block:: none

   mamba install 'dftbplus=*=nompi_*'

If you prefer to install an MPI parallel version you have to
explicitly request it with

.. code-block:: none

   mamba install 'dftbplus=*=mpi_mpich_*'

for MPICH, or

.. code-block:: none

   mamba install 'dftbplus=*=mpi_openmpi_*'

with Open MPI.

There are some differences in functionality between the serial/OpenMP
and MPI versions of the code. The non-MPI version supports more
excited state methods, while the MPI version has better parallelism
for many tasks.

You may want to set up separate environments for the MPI and non-MPI
dftb+ versions::

  conda create --name dftbplus
  conda activate dftbplus
  mamba install 'dftbplus=*=nompi_*'
  conda deactivate

and likewise for your choice of MPI::

  conda create -n dftbplusMPI
  mamba install -n dftbplusMPI 'dftbplus=*=mpi_openmpi_*'

Then list available environments::

  conda info --envs

Additional components like the dptools and the Python API, are
available as separate packages on the same channel. You can install
them with

.. code-block:: none

   mamba install dftbplus-tools dftbplus-python

It is possible to list all of the versions of ``dftb+`` and its
additional components that are available on your platform with:

.. code-block:: none

   mamba search 'dftbplus*' --channel conda-forge


Where to start with the tutorials
=================================

The individual chapters of this document are more or less independent
from each other, so you may go directly to the relevant one for your
interests. However, if you are new to DFTB+, please make sure to work
through the relevant introductory examples in the :ref:`sec-basics`
chapters first.

The recipes are to introduce you to specific functionalities of DFTB+
and so are therefore rather short and focused. Please also always
consult the corresponding sections of the `DFTB+ manual
<http://www.dftbplus.org/fileadmin/DFTBPLUS/public/dftbplus/latest/manual.pdf>`_
for further details and possibilities.

Please note that the example outputs in the recipes may have been
created with older versions of DFTB+ and therefore could differ
slightly in format from output of the most recent code. The
corresponding inputs in the archive should work, without any changes,
with the last stable release of DFTB+.
