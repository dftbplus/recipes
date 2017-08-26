.. highlight:: none

Preparing for an MD calculation
===============================

The initial structure for starting a molecular dynamics simulation should
usually be structurally relaxed. This is to remove excess potential energy which
would rapidly exchange into kinetic energy. 

The example below relaxes the geometry of a molecule (PTCDA) to a local
miniumum, prior to performing molecular dynamics.

.. only:: builder_html
   
   See :download:`the full input
   <../inputs/moleculardynamics/relaxation/dftb_in.hsd>`


Vibrational modes
=================

Once at a structural minimum the quasi-harmonic vibrational modes can be
calculated, for future comparision with the power spectrum of the system.

.. only:: builder_html
   
   See :download:`the full input
   <../inputs/moleculardynamics/vibration/dftb_in.hsd>` and the
   :download:`geometry <../inputs/moleculardynamics/vibration/geom.gen>`

Calculating the modes
~~~~~~~~~~~~~~~~~~~~~

.. only:: builder_html
   
   See :download:`the modes input
   <../inputs/moleculardynamics/vibration/modes_in.hsd>`
