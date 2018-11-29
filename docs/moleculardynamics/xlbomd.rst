.. highlight:: none

******************
Speeding up SCC MD
******************

Instead of standard Born-Oppenheimer (BO) dynamics, there are two methods in use
to propagate the atomic positions and electronic structure, either the
Car-Parrinello scheme or extended Lagrangian (XL) Born-Oppenheimer
dynamics. DFTB+ supports the second of these methods, which (when stable)
produces equivalent dynamics to conventional Born-Oppenheimer molecular dynamics
but with similar (lower) costs to non-SCC DFTB.

This example uses the fast XL-BOMD method with a single SCC step for each
geometry.

.. only:: builder_html
   
   See :download:`the full input <../_downloads/md5.dftb_in.hsd>`

