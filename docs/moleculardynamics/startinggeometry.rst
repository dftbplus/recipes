.. highlight:: none

Preparing for an MD calculation
===============================

.. only :: builder_html or readthedocs

   [Input: `recipes/moleculardynamics/initialstructre/`]

The initial structure for starting a molecular dynamics simulation should
usually be fully structurally relaxed. This is done to remove any excess
potential energy, which would otherwise rapidly convert into kinetic energy (so
heating up the system).

This example relaxes the geometry of a molecule (PTCDA) to a local miniumum,
prior to performing molecular dynamics.

After running the relaxation, the output file ``geo_end.gen`` is now suitable to
use as a starting structure for MD (or for other property calculations).

.. _sec-md-vib-modes:

Vibrational modes
=================

.. only :: builder_html or readthedocs

   [Input: `recipes/moleculardynamics/vibrations/`]

Once at a structural minimum, the quasi-harmonic vibrational modes of the atoms
in the system can be calculated. These can then be compared with the power
spectrum of the system, as determined by molecular dynamics.

The vibrational modes can be obtained from the mass-weighted Hessian matrix of
the system. This can be calculated by ``DFTB+`` from finite difference second
derivatives of the energy with respect to atom positions (technically, the first
derivatives of the forces are actually used).

This derivative calculation can be enabled as an alternative choice of the
Driver::
  
  Driver = SecondDerivatives {
      Delta = 1E-4
  }

The provided example also sets the size of the finite difference step used to
differentiate the energy as 0.0001 atomic units.

The modes are known as quasi-harmonic, since depending on the numerical step
size these derivatives will include a contribution from higher even derivatives
of the function. Hence, for accurate evaluation of harmonic frequencies, the
smallest stable choice of the step size should be used (however, a careful
choice of this value can sometimes be used to include anharmonicity in the
calculated vibrational energies, see for example `Jones and Briddon 1998
<https://doi.org/10.1016/S0080-8784(08)63058-6>`_).

Running the code will produce a file called ``hessian.out`` which contains the
derivative information.

Calculating the modes
~~~~~~~~~~~~~~~~~~~~~

The ``modes`` code can read the Hessian matrix and calculate the vibrational
modes. An example input file is included in this directory using the
``hessian.out`` file that is produced by the DFTB+ calculation using the
supplied ``dftb_in.hsd``.

The general structure of the input for modes is similar to DFTB+, also being in
the hsd format. The code expects a ``modes_in.hsd`` file, which contains blocks
for the geometry, Hessian and the Slater-Koster files (used to get the mass of
the atoms). A typical input looks like::


  # Needs the equilibrium geometry, at which the Hessian had been calculated
  Geometry = GenFormat { 
     <<< geom.gen
  }
  
  DisplayModes = {
    PlotModes = -20:-1 # Take the top 10 modes
    Animate = Yes      # make xyz files showing the atoms moving
  }

  # You need to specify the SK-files, as the mass of the elements is needed
  SlaterKosterFiles = Type2FileNames {
    Prefix = "../../slakos/mio-ext/"
    Separator = "-"
    Suffix = ".skf"
  }
  
  # Include the Hessian, which was calculated by DFTB+
  Hessian = {
    <<< "hessian.out"
  }

  # This file uses the 3rd input format of the modes code
  InputVersion = 3

Here, the code produces animations of several modes, which can be viewed for
example using the `jmol <http://jmol.sourceforge.net/>`_ cross platform viewer
::

   jmol mode_114.xyz &

for the highest frequency stretch mode of this molecule.

The vibrational modes can also be obtained by processing the results of molecular dynamics (see :ref:`sec-md-analysis`)
