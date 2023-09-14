Solvation effects
=================

Chemistry mainly happens in the condensed phase in solvation.
To model this environment effect we can employ an implicit solvation model.


Partition coefficients
----------------------

We will again use the `tributyl phosphate <https://pubchem.ncbi.nlm.nih.gov/compound/31357>`__ system for this chapter starting from the lowest conformer computed at DFT level of theory.\ :cite:`salthammer2022`

.. tab-set::

   .. tab-item:: coordinates
      :sync: xyz

      .. literalinclude:: data/tbnp/crenso-woctanol.xyz
         :caption: dft-lowest.xyz

   .. tab-item:: image
      :sync: png

      .. image:: data/tbnp/crenso-woctanol.png
         :scale: 50%


For a partition coefficient we need to evaluate the solvation free energy for the two phases, we will evaluate the environmental relevant octanol-water partition coefficient (log\ *K*\ :sub:`ow`).
The formula for the partition coefficients can needs the free energies for the molecule in water and octanol

.. math::

   \log K_\text{ow} = (G_\text{water} - G_\text{octanol}) / (k_B T) \log e

.. tip::

   The value :math:`\log e / (k_B T)` at 298 K is 460.200 E\ :sub:`h`\ :sup:`–1`.


We will setup a calculation using the ALPB solvation model,\ :cite:`ehlert2021-1` using the GFN2-xTB and DFTB3-D4 using the inputs below for water.

.. tab-set::

   .. tab-item:: GFN2-xTB
      :sync: xtb

      .. literalinclude:: data/gfn2-xtb-water.hsd
         :caption: dftb_in.hsd
         :language: shell
         :emphasize-lines: 7

   .. tab-item:: DFTB3-D4
      :sync: dftb

      .. literalinclude:: data/dftb-3ob-water.hsd
         :caption: dftb_in.hsd
         :language: shell
         :emphasize-lines: 33


.. admonition:: Exercise
   :class: info

   Perform the calculation for the free energy of water and octanol.
   You might need to adjust the parameter file name for the octanol solvent parameters.

   How well does the calculated partition coefficient compare to the experimental value of 4.0?
   Is the value well reproduced, which effects are not accounted for?


Summary
-------

.. admonition:: You learned...
   :class: important

   - to setup an implicit solvation model with DFTB and xTB
   - perform an optimization accounting for solvation effects
   - calculate a partition coefficient from free energies
