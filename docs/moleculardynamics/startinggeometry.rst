.. highlight:: none

Preparing for an MD calculation
===============================

The initial structure for starting a molecular dynamics simulation should
usually be structurally relaxed. This is to remove excess potential energy which
would rapidly exchange into kinetic energy. 

The example below relaxes the geometry of a molecule (PTCDA) to a local
miniumum, prior to performing molecular dynamics.

.. only:: builder_html
   
   See :download:`the full input <../_downloads/md1.dftb_in.hsd>` and the
   initial :download:`starting geometry <../_downloads/md1.geom.gen>`


Vibrational modes
=================

Once at a structural minimum the quasi-harmonic vibrational modes can be
calculated, for future comparision with the power spectrum of the system.

.. only:: builder_html
   
   See :download:`the full input <../_downloads/md2.dftb_in.hsd>` and the
   geometry :download:`geometry <../_downloads/md2.geom.gen>`

Calculating the modes
~~~~~~~~~~~~~~~~~~~~~

.. only:: builder_html
   
   See :download:`the modes input <../_downloads/md2.modes_in.hsd>`
