.. _sec-interfaces-ase-neb:

#########################
Nudged Elastic Band - NEB
#########################

The nudged elastic band (NEB) method is used to obtain saddle points as well as
minimum energy paths between known initial and final states. Intermediate images
are introduced (by interpolation), for which local energy minima are calculated
under the constraint of equal spacing between neighboring images (except for the
`climbing image` modification).
For a detailed description of all the features provided by ASE's NEB class
please consult the corresponding documentation (`NEB class
<https://wiki.fysik.dtu.dk/ase/ase/neb.html>`_).

The umbrella inversion of a simple |NH3| (ammonia) molecule serves as a
prominent example for NEB energy barrier calculations.

.. |NH3| replace:: NH\ :sub:`3`\

.. toctree::
   :maxdepth: 1

   via File-IO <fileio/fileio.rst>
   via Socket-Communication <sockets/sockets.rst>
