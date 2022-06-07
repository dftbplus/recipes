.. _sec-surfaces:

########################
Semi-infinite structures
########################

.. only :: builder_html or readthedocs

[Input: `recipes/boundaryconditions/surfaces/Si100`]

The same methods used for :ref:`sec-transport` calculations can also
be applied to simulate structures like the surface of a semi-infinite
crystal, or the ends of polymer chains, nano-wires or nanotubes.

Periodic boundary slab
======================

A supercell input (`Si100_slab.hsd`) is provided for the Si (100)
surface with buckled dimers on top of a (thin) bulk-like region. In
this case there are (obviously?) two surfaces on the top and bottom of
the slab, and most of the atoms in between are kept fixed during
optimisation.

There are several potential problems with this geometry

  * There are two surfaces, which are potentially close enough to
    influence each other, either through the bulk of the cell or
    across the (artificial) vacuum layer separating to the next cell.

  * Charge transfer from the bulk to the surface might be happening in
    the real system on top of a thick crystal, but this model
    structure is *forced* to be **charge neutral**

  * The dipole moment of the slab can cause problems with
    self-consistency (this particular example has inversional
    symmetry, so prevents a dipole occurring).

Green's function embedding
==========================

The alternative is to represent a semi-infinite bulk structure by
replacing the slab with a `contact` made of bulk material, and
simulating atoms on top of this. There are more details in the
:ref:`sec-transport` section, but this calculation proceeds in two
steps ::

  1. Perform a special periodic calculation to evaluate the bulk
     contact Fermi energy and atomic potentials (see the provided
     'bulk.hsd' input)

  2. Use the resulting `shiftcont_bulk.bin` file to calculate the
     surface on top of this semi-infinite bulk (see the provided
     'surface.hsd' input)

There are several unusual features in both calculations, compared to a
supercell/periodic input

  * The geometry (see `Si_2x1.gen`) contains the surface atoms first,
    then layers of the bulk material (in a specific order).

  * Both calculations have an extra `Transport {` block, specifying
    which atoms are bulk and which are surface (`contact` and `device`
    respectively).

Other parts of the bulk calculation input looks like a usual periodic
calculation (k-points), but the surface calculation uses some other
features::

  Solver = GreensFunction{}
  Electrostatics = Poisson {
    MinimalGrid [Angstrom] = 0.3 0.3 0.3
    PoissonThickness [AA] = 50
  }

This is now using the libNEGF Green's function solver for an open
system, instead of the usual method of finding
eigenvalues. Additionally a Poisson solver is used, instead of the
usual Ewald method, and in this example it solves the electrostatics
in a region extending 50 Angstroms above the top of the surface.

The k-points are also somewhat different, the structure is periodic in
the *xy* plane, so these are like a usual supercell. But in the *z*
direction there should not be any k-point integration, as this is open
(a single point with a shift of 0.0 is used)::

  KPointsAndWeights = SupercellFolding {
    10 0 0
    0 10 0
    0 0 1
    0.5 0.5 0.0
  }

(superficially this looks similar to the slab geometry input, but in
that case the vacuum 'slab' between the silicon surfaces removes the
need for k-point integration along `z`).

The geometry optimisation is also limited to the top layers of the
structure (atoms 1 to 16). Unlike the slab calculation, the reason for
this is that the analytical forces (evaluated by the Helmann-Feynman
theorem) become incorrect close to the contact, but are accurate
beyond a few layers of bulk-like material in the device region.

Hence, it's important to include a buffer of geometrically fixed
bulk-like material in the `device region` below the surface atoms.
