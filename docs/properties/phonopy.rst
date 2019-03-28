Phonon calculations with phonopy
================================

The `phonopy <http://atztogo.github.io/phonopy/>`_ code can calculate a range of
harmonic and quasi-harmonic vibrational properties and from version 2.0 onwards
supports DFTB+. Information about how to install phonopy is `available
<http://atztogo.github.io/phonopy/install.html>`_.

The below examples were tested with phonopy v2.1.2.

Phonon band structure
~~~~~~~~~~~~~~~~~~~~~

.. only :: builder_html or readthedocs

   [Input: `recipes/properties/phononbs/`]

The diamond lattice has very high symmetry, hence a phonon band structure can be
obtained with a single calculation. A conventional unit cell with relaxed
lattice constant is provided in `geo.gen` (the phonopy input assumes this is the
name of the suplied starting geometry).

The DFTB+ input needs to calculate the atomic forces and also to write a results
tag file, hence the `dftb_in.hsd` input contains::
  
  Analysis = {
    # required option for phonopy
    CalculateForces = Yes
  }
  Options = {
    # Required options for storing data needed by phonopy
    WriteResultsTag = Yes
  }

The (single) distorted geometry and the required `phonopy_disp.yaml` file is then
generated with the command ::

   phonopy -d --dim="4 4 4" --dftb+

This constructs a 4x4x4 supercell of the primitive cell and saves the
undistorted and distorted supercells as `geo.genS` and `geo.genS-001`
respectively. Note that you should test the required number of repeats in
supercell required to reach convergence of the band structure (or other
properties of interest).

For the single (in this case) relevant `geo.genS-*` file, calculate the DFTB
forces on the atoms, retaining the resulting `results.tag` file. This example
uses the mio Slater-Koster paramters for carbon.

Then create the required `FORCE_SETS` file ::

  phonopy -f results.tag --dftb+  ...

This assumes that the `results.tag` and `phonopy_disp.yaml` files are in the
same directory (for more complex cases the calcuation should be run in the same
directory as the phonopy files and the path to directories containing the DFTB+
output files).

Then specify the path in the Brillouin zone you are interested in (see the
`phonopy documentation
<https://atztogo.github.io/phonopy/setting-tags.html#band-structure-related-tags>`_). Then
post-process the phonopy data, providing the dimensions of the the supercell
repeat either on the command line or in the settings file (a `DIM` file)::

   phonopy -p band.conf --dim="4 4 4" --dftb+

Finally create a band structure in gnuplot format ::

  phonopy-bandplot --gnuplot band.yaml > band.dat

The resulting band structure for the mio carbon model is shown in
:numref:`_fig_phonopy_diamond`.

  .. _fig_phonopy_diamond:
  .. figure:: ../_figures/properties/phonopy/band.png
     :height: 40ex
     :align: center
     :alt: Resulting phonon band structure for diamond.
