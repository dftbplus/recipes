
.. _specifying-geometry:

Specifying the geometry
=======================


Partitioning the system
-----------------------

In the simplest case, transport calculations can be performed on an atomistic structure
comprising two semi-infinite contacts. The contacts are
usually 1D-periodic systems being terminated on one end by the device.

In order to carry out a transport calculation with DFTB+, the system must be carefully
partitioned by the user and the structure must contain, 

  A. The extended molecule (central region),
  #. Two *Principal Layers* (see below) of the first contact,  
  #. Two *Principal Layers* of the second contact.

  .. _fig_transport_geometry_partitioning:
  .. figure:: ../_figures/transport/geometry/device.png
     :height: 20ex
     :align: center
     :alt: Partioning of a device

     Partitioning of a simple system into extended molecule and
     perfect contacts.

The extended molecule contains the atomistic device itself *plus* those parts of
the attached contacts, which are directly influenced by the presence of the
device (we can call them *surfaces*).  Each contact, on the other hand, contains
those parts of the actual contacts, which are far enough from the device and are
not affected by it. Examples below will better clarify this point.

  .. _fig_transport_geometry_principal_layers:
  .. figure:: ../_figures/transport/geometry/layers.png
     :height: 20ex
     :align: center
     :alt: The subdivision of the structure into layers

     Subdivision of a general structure into principal layers (PL), two for each
     contact and an arbitrary number within the extended molecule.


Partitioning the system and contacts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most important concept to bear in mind is that of *principal layers* (PLs).
These are defined as contiguous groups of atoms that have interaction only with
atoms of adjacent PLs.  In practice these layers must contain a sufficient
number of atoms in order to ensure that Hamiltonian and overlap interactions
with second neighbouring PLs vanish or can be considered negligible.  The
subdivision into PLs is essential for the definition of the two contacts and
becomes useful in the computation of all Green's functions that exploit a
recursive algorithm. (See [PPD2008]_.)

As shown in the :numref:`fig_transport_geometry_principal_layers`, the extended
molecule can contain an arbitrary number of PLs (>=1), but the layers must
follow a sequential ordering.  The ordering of the PLs follow directly from the
ordering of the atoms in the DFTB+ structure.

Typically it is convenient to create the structure and then sort the atom along
the transport coordinate, before partitioning the system. In other cases, when
nanowires are constructed it is convenient to repeat a PL unit for the desired
length of the extended molecule.

The (perfect) contacts must be defined by two principal layers. 
Differently from the layers of the extended molecule, those defining the contacts must 
follow additional rules:

  A. The two principal layers of a given contact must be identical.
  #. They must be rigidly shifted images of each other.
  #. The first PL must be the closest to the device region. 


Numbering the atoms in the system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to ensure the above prescriptions, the numbering of the atoms must
follow a precise ordering.

The atoms of the central region must be specified first, before the atoms of the
contacts (see example below). The atoms of the contacts follow the central
region.  The atoms of each contact must be specified continously. (You specify
either all atoms of the left contact first, and then those of the right one, or
vica versa.)

The numbering inside the first principal layer defining each contact is
arbitrary, but the *same numbering must be applied to the second PL*. Thus, the
indices of the corresponding atoms in the two contact PLs must differ by the
same value, and the coordinates must differ by the same translation vector.
These prescriptions are checked by the code and an error message is issued if
the two PLs do not conform. It is important to remember to place the first PL
closer to the device, i.e. the principal layer closer to the device must contain
the atoms with lower indices.

  .. _fig_transport_geometry_numbering:
  .. figure:: ../_figures/transport/geometry/device_numbering.png
     :height: 20ex
     :align: center
     :alt: atomic_numbring

     Example for numbering the atoms. Atoms of the device have the
     lowest indices, followed by the atoms of each contact,
     respectively. The two principal layers of both contacts
     are shifted images of each other.

Supercell structures
^^^^^^^^^^^^^^^^^^^^

The code can compute transport on structures that have a periodicity in the transverse directions 
(with respect to transport). In this case the structure must be defined as a supercell and the 
rules listed above apply to each cell. 
The real-space Poisson solver of DFTB+ limits the supercell lattices to orthorombic types 
(all angles between supercell vectors must be 90 degrees).
In fact the supercell is always defined 3-dimensional, and the user should ensure that the 
dummy lattice vector along the transport direction is long enough to avoid superpositions between
images. 

The current version of DFTB+ only supports one supercell definition for the
entire system, including central region and contacts. It may be redundant to
observe that in this case the two contacts must be of the same periodicity.
 
Transport Block
^^^^^^^^^^^^^^^

The geometry must be defined in the transport block as specified in the following example::
 
    Transport {
      Device {
        AtomRange = 1 24 
      }
      Contact {
        Id = "source"
        AtomRange = 25 44 
      }
      Contact {
        Id = "drain"
        AtomRange = 45 58
      }
    }
