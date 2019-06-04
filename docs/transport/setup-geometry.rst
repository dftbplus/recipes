.. _setup-geometry:


Setup Geometry utility
======================

This feature of dftb+ helps preparing the structure as prescribed in the section before.
The tool can be invoked as a separate program `setupgeom` and requires an input file
with name `setup_in.hsd`. 
Here we give an example that works most easily with the help of the tool jmol 
(jmol link `http://jmol.sourceforge.net/`_). 
The example discussed here corresponds to the molecular junction that will be presented 
in a section below.

Input file
^^^^^^^^^^

The tool works by setting up an hsd input file like the following::

  Geometry = GenFormat{
    <<< "input.gen"
  }
  Transport{
    Contact{
      Id="s"
      SpecifiedPLs = 2
      Atoms = {46:61 70:85 94:109}
      ContactVector [Angstrom] = 4.94 0.0 0.0  
    }
    Contact{
      Id="d"
      SpecifiedPLs = 2
      Atoms = {118:133 142:157 166:181}
      ContactVector [Angstrom] = 4.94 0.0 0.0 
    }
    Task=SetupGeometry{}
  }

``Atoms``. It is necessary to specify all atoms belonging to the two contacts. 

``ContactVector``. It is also necessary to specify the contact supercell vectors
Note that contacts must extend along either `x`, `y` or `z`. 

``SetupGeometry``. This is the ``Task`` that must be invoked in order to 
perform the geometry preparation for transport.

``SpecifiedPLs``. It is used to specify how many contact PLs are provided. 
Possible values are 1 or 2. See below for details.

Selecting atoms using jmol  
^^^^^^^^^^^^^^^^^^^^^^^^^^

Any tool can be used for this purpose, eventually also a manual selection. However jmol 
serves very well for this task. Unfortunately, jmol does not recognize gen files, so you 
have to convert the starting geometry into a format readable by jmol. Usually xyz is best 
to this purpose. Use `gen2xyz` (or any other conversion tool such as `babel` or `ASE`). 

Open the structure in jmol and select the atoms of contact 1. 
If the contact already contains 2 PLs select both. At this stage it is likely that the 
two PLs are not following the correct numbering, namely they are **not** two shifted 
exact copies of each other. No worries! ``SetupGeometry`` will reorder the atom indices
of the two PLs in the correct way. 
In some cases you might have only 1 PL per contact. The tool can duplicate this PL 
as required by the input geometry. If this happens add the keyword
``NumPLsDefined = 1`` in the relevant ``Contact`` blocks.

  .. _fig_transport_setup-geometry_sel:
  .. figure:: ../_figures/transport/setup-geometry/atom-selected.png
     :width:  100%
     :align: center
     :alt: Screenshot showing the selected atoms of contact 1 

     Jmol screenshot showing the selected contact atoms 
 
In :numref:`fig_transport_setup-geometry_sel` it is possible to see the geometry to be 
processed for reordering with selected atoms belonging to contact 1.

Different strategies can be used to select the contact atoms in jmol. The easiest is 
probably using the `select` tool and using mouse (see :numref:`fig_transport_setup-geometry_tool`). 
Orient the molecule anduse the tool by holding SHIFT + LEFT Mouse Button,  
dragging the mouse to include all contact atoms. 

  .. _fig_transport_setup-geometry_tool:
  .. figure:: ../_figures/transport/setup-geometry/tool.png
     :width:  80%
     :align: center
     :alt: Screenshot of Jmol showing selection tool 

     Jmol screenshot showing selection tool

NOTE: Initially, when you click on the selection tool, all atoms will be selected 
and will appear highlighted. 
Choose the menu ``Display -> Select -> none`` to unselect all atoms.
Alternatively, open the ``Jmol Script Console`` and type::
  
  $ select none

Now you can select the contact atoms and then type in the console::
  
  $ print {selected}
  ({45:60 69:84 93:108})

The selected atoms are shown in a compact syntax that can be directly 
copy-pasted into ``setup_in.hsd``. 
**NOTE that this jmol command show atom indices starting from 0 and not from 1**. 
In this case use the following syntax in the ``setup_in.hsd`` input file::
  
  Atoms [zeroBased] = {45:60 69:84 93:108}

where the modifier `zeroBased` tells the code that atom indices start from 0. 
Repeat a similar sequence of commands for the other contact.

``ContactVector`` is needed so the code can understand the direction of the contact
and the supercell periodicity. Use the `measurements` tool of jmol in order to 
get the vector length (See :numref:`fig_transport_setup-geometry_sel`).

The user should provide the Slater-Koster files so the code can elaborate the
correct cutoff distances. These are specified in the same way as for the dftb+ code:: 

    Task = SetupGeometry{
      SlaterKosterFiles = Type2FileNames{
         prefix =  "PATH/"
         separator  = "-"
         suffix  = ".skf"
       }
    }    

The following behaviour are relevant.
``SpecifiedPLs = 2``: In this case `setupgeom` reorders the second PL and check that 
the distance between second-neighbour PLs is larger than the cutoff. An error is 
shown if this is not the case.

``SpecifiedPLs = 1``: In this case `setupgeom` build as many additional PLs as 
needed to fulfill the contact requirements.

In both cases the device region is further layered into PLs for the efficient 
iterative Green's function algorithm.
In most cases the SK tables have a rather large cutoff, extending as long as all 
Hamiltonian matrix elements are below 1e-5 a.u. (about 1 meV). 
In order to make transport calculations a little faster it is possible to shorten
slightly the SK cutoffs. A small decrease easily results in PLs with half of the
original value (the contact must be periodic) hence faster calculations, with very
small effect on the final results (e.g., transmission, ldos, currents).
The SK cutoff can be set with the block `TruncateSKRange` (see also dftb+ manual)::

  Transport{
    Task = SetupGeometry{
        TruncateSKRange = {
           SKMaxDistance [AA] = 5.0
           HardCutOff = Yes
        }
    } 
  }

Clearly in doing this, accuracy is traded for speed. In the case of C-C interactions,
the cutoff distance is about 5.17 Angstrom that is quite comparable with a 
smalled cutoff of 5.0 Angstrom. In any case the user should check and validate 
the results when selecting this option.

Once the input is ready convert the structure to `whatever.gen` and run setupgeom
As output you will find the structure ``processed.gen`` prepared for transport 
calculations and a file ``transport.hsd`` containing the ``Transport`` block
needed for the following contact calculations::

  Transport{
    Device{
      FirstLayerAtoms={  1 25 40 50 60 76 }
      AtomRange= 1 95
    }
    Contact{
      AtomRange= 96 143
    }
    Contact{
      AtomRange= 144 191
    }
  }

For consistency, the user should specify exactly the same `SKMaxDistance` in the input
file of dftb+.






















