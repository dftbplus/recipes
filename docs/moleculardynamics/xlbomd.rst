.. highlight:: none

******************
Speeding up SCC MD
******************

Instead of standard Born-Oppenheimer dynamics there are two methods in use to
propogate the atomic positions and electronic structure, either the
Car-Parrinello scheme or extended Lagrangian Born-Oppenheimer dynamics. DFTB+
supports the second of these methods, which (when stable) produces equivalent
dynamics to conventional Born-Oppenheimer molecular dynamics with similar costs
to non-SCC DFTB.

.. only:: builder_html
   
   See :download:`the full input <../_downloads/md5.dftb_in.hsd>`

