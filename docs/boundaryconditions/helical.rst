.. highlight:: none

******************
Helical geometries
******************

.. only :: builder_html or readthedocs

   [Input: `recipes/scripts/repeatGen.pl`]

This section describes calculations for geometries with structures
which can be described using helical (`objective`) boundary conditions.

Currently DFTB+ only supports non-SCC calculations with a combination
of helical and C\ :sub:`n` rotational symmetries.

A gen format for helical geometries looks something like::

    2  H
    C
    1 1    0.2756230044E+01    0.2849950460E+01    0.1794011798E+01
    2 1    0.2756230044E+01    0.2849950460E+01    0.3569110265E+00
    0 0 0
    0.2140932670E+01 18.0 10

This example is for a (10,0) carbon nanotube with a repeat
distance along the tube of 2.141 Angstroms and the next layer of atoms
is rotated by 18 degrees, with an additional C\ :sub:`10` rotational
symmetry. In this case the structure rotational repeats every 360 / 10
= 36 degrees, which is commensurate with twice the 18 degree helical
twist.

  .. figure:: ../_figures/boundaryconditions/helical/tube.svg
     :height: 40ex
     :align: center
     :alt: Repeated nanotube a) along and b) across the tube axis

     Repeated nanotube a) along and b) across the tube axis with the 2
     atom objective cell marked.

The supplied script `repeatGen.pl` can extend helical structures to
make larger sections of the geometry. Eventually this functionality
will be included in the dp_tools package, but for the moment these
repeats need to be added to the end of the structure for this older
script ::

  2  H
  C
  1 1    0.2756230044E+01    0.2849950460E+01    0.1794011798E+01
  2 1    0.2756230044E+01    0.2849950460E+01    0.3569110265E+00
  0 0 0
  0.2140932670E+01 18.0 10
  3 10

Then `repeatGen.pl struct.gen > tmp.gen` to make a repeated geometry
with the full ring geometry (the usual gen2xyz script *can* process
this: `gen2xyz tmp.gen`).

This example is can also be represented as a compact supercell, but
what happens if the twist angle is not 18 degrees (i.e. a simple
rational ratio with 36 degrees)?


Ribbon structures
-----------------

.. only :: builder_html or readthedocs

    [Input: `recipes/boundaryconditions/helical/`]

Symmetric ribbon geometries can also be represented with helical
boundary conditions. Both a conventional supercell and a helical
geometry for the same structure are supplied. The helical case has a
180 degree rotation symmetry to make the other side of the ribbon. The
supercell has simple translational symmetry along the ribbon with a
cell containing twice as many atoms.

  .. figure:: ../_figures/boundaryconditions/helical/ribbon.svg
     :height: 40ex
     :align: center
     :alt: Ribbon geometry as a supercell and one possible helical
           cell consisting of half of the supercell.

The k-points for a helical geometry represent sampling along the screw
axis (and a shift)::

  KPointsAndWeights = HelicalUniform {20 0.5}

In this case generating 20 points along the reciprocal space (and
automatically x2 points around the C\ :sub:`2` rotation). In this
case, the translational operation is directly equivalent to the
supercell, so both calculations give equivalent band structures.

  .. figure:: ../_figures/boundaryconditions/helical/helicalbs.svg
     :height: 40ex
     :align: center
     :alt: Band structure of helical cell, showing band that are
           symmetric (g) and anti-symmetric (u) under a C2 rotation. In
           this case, the band structure for a supercell is identical
           if the labels are ignored.

     Band structure of helical cell, showing band that are
     symmetric (g) and anti-symmetric (u) under a C\ :sub:`2`
     rotation. In this case, the band structure for a supercell is
     identical if the labels are ignored.

However, in the case of the helical structure, the k-points for the
second C\ :sub:`2` symmetry operation allow us to separate bands that
are symmetric (`g`) and anti-symmetric (`u`) with respect to that
operation.

We can of course optimise the atomic coordinates, but also we can
twist the structure, by changing the helical angle, to produce
geometries that cannot be easily represented as a supercell (a 1
degree twist requires 360 repeats along the axis to be equivalent
under translation).

Try twisting the ribbon and relaxing it's geometry, does the energy go
up or down compared to a flat ribbon?
