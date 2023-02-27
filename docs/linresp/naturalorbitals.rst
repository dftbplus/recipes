.. highlight:: none

***********************************************
Visualization of natural orbitals with waveplot
***********************************************

[Input: `recipes/linresp/natural-orbitals`]

The majority of excited states are not well described by the
transition of an electron from a single Kohn-Sham orbital to another
single Kohn-Sham orbital. In other words, the excited state is in
general a superposition of many determinants. In such a situation, the
excited state weight in the file EXC.DAT will be smaller than
one. Natural transition orbitals (NTO) allow for a compact
representation of the excitation in such a scenario :cite:`luzanov1976application`. They are obtained
from a singular value decomposition of the transition density matrix
and deliver `hole` and `electron` orbitals for a given excited
state. We discuss their visualization in the context of the
Benzene-Tetracyanoethylene dimer we looked at :ref:`earlier <benz_dimer>`. Strictly
speaking, natural orbitals are not really required for this system,
since for example the lowest excited state is well described by a HOMO
to LUMO transition, without contributions from other single-particle
transitions. We still discuss the NTO feature for this example, because
it allows us to compare the results with the earlier calculations.

We first run DFTB+ with the following modified input:

.. literalinclude:: ../_archives/recipes/linresp/natural-orbitals/dftb_in.hsd
   :lines: 21-38

Here ``CalculateForces`` requests to compute excited forces for state
``StateOfInterest``. During this calculation also the NTO for that state
are created and written to a file (always called `excitedOrbs.bin`), if
``WriteEigenvectors`` in the ``Casida`` block is set to ``Yes``. In the
``Options`` block, the keyword ``WriteDetailedXml`` is required for the
plot. In the ``Analysis`` block, ``WriteEigenvectors`` advises the code to
print out the ground state molecular orbitals also.

Let us look at an excerpt of the generated file `detailed.xml`::

 <excitedoccupations>
  <spin1>
   <k1>
    -1.00000341372109 -5.201183337079720E-002 -5.170515789354190E-002 .....
    .....................
    5.170515780471528E-002 5.201183337096685E-002 1.00638172624593

   </k1>
  </spin1>
 </excitedoccupations>

We see the state occupations :math:`n_{i}` of the NTO. Negative values
indicate hole orbitals and positive values electron orbitals. Most
excited states (unless you have degeneracy) have only one important
(:math:`n_{i}\approx -1.0`) hole and one important
(:math:`n_{i}\approx 1.0`) electron NTO, even though the CI expansion
might include a large number of single-particle transitions.

As seen in the section  :ref:`here <sec-basics-waveplot>`, we can now visualize these orbitals using the waveplot code. The input file `waveplot_in.hsd` has this form:

.. literalinclude:: ../_archives/recipes/linresp/natural-orbitals/waveplot_in.hsd
   :emphasize-lines: 4,13

Important is here the line ``EigenvecBin`` which reads the created file
with natural transition orbitals instead of the ground state MO. The
keyword ``PlottedLevels`` in this example requests the first and last
orbital to be plotted (note that the curly bracket contains the
indices of the NO, not the occupation). These are exactly the natural
transition orbitals with highest occupation. If you do the plot, you
will realize that the hole orbital is essentially Kohn-Sham state 37
(the HOMO), while the electron orbital is Kohn-Sham state 38 (the
LUMO).
