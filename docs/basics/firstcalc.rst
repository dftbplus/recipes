.. highlight:: none

****************************
First calculation with DFTB+
****************************

This chapter should serve as a tutorial guiding you through your first
calculation with DFTB+. As an example, the equilibrium geometry of a water
molecule is calculated. In order to execute the calculation described here, you
will need the Slater-Koster set `mio`, which you can download from the DFTB.ORG
(http://www.dftb.org).

The description here is based on DFTB+ version 17.1, the input/output files in
other versions may slightly differ from those shown here.

Providing the input
===================

First you have to create the input for the code. DFTB+ accepts either
Human-readable Structured Data (HSD) or eXtended Markup Language (XML) input. In
this tutorial HSD will be used, and in this case the input file must be called
`dftb_in.hsd`.  The input file used in this howto looks as follows::

  Geometry = GenFormat { 
  3 C 
    O H 
  
    1 1  0.00000000000E+00 -0.10000000000E+01  0.00000000000E+00
    2 2  0.00000000000E+00  0.00000000000E+00  0.78306400000E+00
    3 2  0.00000000000E+00  0.00000000000E+00 -0.78306400000E+00 
  }
  
  Driver = ConjugateGradient {
    MovedAtoms = 1:-1
    MaxForceComponent = 1E-4
    MaxSteps = 100
    OutputPrefix = "geom.out"
  }

  Hamiltonian = DFTB {
    SCC = Yes
    SlaterKosterFiles {
      O-O = "O-O.skf"
      O-H = "O-H.skf"
      H-O = "H-O.skf"
      H-H = "H-H.skf"
    }
    MaxAngularMomentum {
      O = "p"
      H = "s"
    }
    Filling = Fermi {
      Temperature [Kelvin] = 0.0
    }
  }
  
  Options {}
  
  Analysis = {
    CalculateForces = Yes
  }
  
  ParserOptions {
    ParserVersion = 5
  } 

The order of the specified blocks in the HSD input is arbitrary. You are
free to capitalise the keywords and physical units as you like, since they are
case-insensitive. This is not valid, however, for string values, especially if
they are specifying file names.


Geometry
--------

The ``Geometry`` block contains the types and coordinates of the atoms in
your system.  The geometry of the system in the sample input file is provided in
the so called "gen" format, which was the traditional geometry input format of
the DFTB method. The formal description of this format can be found in the
DFTB+ manual.  The current example::

  Geometry = GenFormat {
    3  C                   # 3 atoms, non-periodic cluster
     O H                   # Two elements, 1 - O, 2 - H
    #  Index Type  Coordinates
         1    1    0.00000000000E+00  -0.10000000000E+01   0.00000000000E+00
         2    2    0.00000000000E+00   0.00000000000E+00   0.78306400000E+00
         3    2    0.00000000000E+00   0.00000000000E+00  -0.78306400000E+00

This specifies a cluster of 3 atoms of the type O and H (cluster
geometries having open boundary conditions, as distinct from periodic
supercell geometries). The Cartesian coordinates of the atoms in the
"gen" format are given in Angstroms.  The first column of integers
contains sequential numbering of the atoms in the system (the actual
values are ignored by the parser).  The second column contains the
type of each atom, given as the position of the appropriate element in
the element list of the second line of the "gen" data.  The
``GenFormat{}`` is not the only way to specify the geometry, you
should check the manual for other methods.

As demonstrated above, it is possible to put arbitrary comments in the
HSD input after a hash-mark (``#``) character. Everything between this
character and the end of the current line is ignored by the parser.

Very often, the geometry is stored in an external file. To save you
the copying and pasting from that file into the `dftb_in.hsd` file,
you can use the file inclusion feature of the HSD format::

  Geometry = GenFormat {
    <<< "input_geometry.gen"
  }

The ``<<<`` operator includes the specified file as raw data. (The file is not
checked for any HSD constructs.) In the example above, the file
`input_geometry.gen` *must* be in gen format.

Driver
------

After having specified the geometry of your system, you should decide
what DFTB+ will do with that geometry. The ``Driver`` environment
determines how the geometry should be changed (if at all) during the
calculation. If you only would like to make a static calculation, you
must either set it to an empty value like::

  Driver = {}   # Empty value for the driver

or omit the ``Driver`` block completely from `dftb_in.hsd`.

In the current example::

  # Do conjugate gradient optimisation
  Driver = ConjugateGradient {
    MovedAtoms = 1:-1               # Move all atoms in the system
    MaxForceComponent = 1.0e-4      # Stop if maximal force below 1.0e-4
    MaxSteps = 100                  # Stop after maximal 100 steps
    OutputPrefix = "geom.out"       # Final geometry in geom.out.{xyz,gen}
  } 

the molecule is relaxed using the conjugate gradient method. The
entire range of atoms from the first (atom 1) until and including the
last (-1) is allowed to move. Instead of ``1:-1`` you could also have
written::

  MovedAtoms = 1:3               # Atoms from the 1st until the 3rd

or::

  MovedAtoms = O H               # Select O and H atoms.

or::

  MovedAtoms = 1 2 3              # Explicitely listing all atom numbers.


In our case the geometry optimisation continues as long as the maximum
component of the force acting on the moving atoms is bigger than 1e-4
atomic units (Hartree per Bohr radius).  Numeric values are by default
interpreted to be in atomic units. However the HSD format offers the
possibility of using alternative units by specifying a unit modifier
before the equals sign. This is given in square brackets. For example
instead of the original atomic units, you could have used::

  MaxForceComponent [eV/AA] = 5.14e-3    # Force in Electronvolts/Angstrom

or::

  MaxForceComponent [Electronvolt/Angstrom] = 5.14e-3

see the manual for the list of accepted modifiers.

The ``MaxSteps`` keyword specifies the maximum number of geometry
optimisation steps that the program can take before stopping, even if
the specified tolerance for the maximal force component have not been
achieved by that stage of the calculation.

Finally, the ``OutputPrefix`` keyword specifies the name of the file
to be written that will contain the present geometry during the
optimisation (and then the final geometry at the end of the
calculation). The geometry is written in gen and xyz formats to the
files obtained by appending ".gen" and ".xyz" suffixes to the
specified name (`geom.out.gen` and `geom.out.xyz` in our case.)
The `dptools` package on the `DFTB+ website
<http://www.dftb-plus.info>`_ contains scripts to convert between the
gen and the xyz formats (and various other formats).


Hamiltonian
-----------

You have to decide upon the model used to describe your system in
order to calculate its properties. At the moment DFTB+ eases this
decision quite a lot, since it currently only supports types of
Density Functional based Tight Binding Hamiltonians (with some
extensions). In our example, the chosen self-consistent DFTB
Hamiltonian has the following properties::

  Hamiltonian = DFTB {                 # DFTB Hamiltonian
    SCC = Yes                          # Use self consistent charges
    SlaterKosterFiles {                # Specifying Slater-Koster files
      O-O = "O-O.skf"
      O-H = "O-H.skf"
      H-O = "O-H.skf"
      H-H = "H-H.skf"
    }
    MaxAngularMomentum {               # Maximal l-value of the various species
      O = "p"
      H = "s"
    }
    Filling = Fermi {                  # No electronic temperature
      Temperature [Kelvin] = 0.0
    }
  } 

In this example the charge self-consistent DFTB (SCC-DFTB) method is
used for the electronic structure (and calculating the total energy,
forces, etc.). This method includes the effect of charge transfer
between atoms of the system. In order to find the final ground state
of the system it has to iteratively solve the system, until the atomic
charges are self-consistently converged. Convergence is reached if the
difference between the charges used to build the Hamiltonian and the
charges obtained after the diagonalisation of the Hamiltonian is below
a certain tolerance (the default is 1e-4 electrons, but can be tuned
with the ``SCCTolerance`` option). If this level of convergence is not
reached within a certain number of iterations, the code calculates the
total energy using the charges obtained so far and stops with an
appropriate warning message. The maximal number of scc-iterations is
by default 100, but can be changed via the ``MaxSCCIterations``
option.


The tabulated integrals (together with other atomic and diatomic
parameters) necessary for building the DFTB Hamiltonian are stored in
the so called Slater-Koster files. Those files always describe the
interaction between atom pairs. Therefore, you have to specify, for
each pairwise combination of chemical elements in your system, the
corresponding Slater-Koster file::

  SlaterKosterFiles = {               # Specifying Slater-Koster files
    O-O = "O-O.skf"
    O-H = "O-H.skf"
    H-O = "O-H.skf"
    H-H = "H-H.skf"
  }

If you use a consistent file naming convention, you can avoid typing
all the file names by specifying only the generating pattern. The
input::

  SlaterKosterFiles = Type2FileNames {   # File names with two atom type names
    Prefix = ""             # No prefix before first type name
    Separator = "-"         # Dash between type names
    Suffix = ".skf"         # Suffix after second type name
  }

would generate exactly the same file names as in the example above. If
the Slater-Koster files are in a different directory from the
`dftb_in.hsd` input file, you can specify the path as a prefix::

  SlaterKosterFiles = Type2FileNames {    # File names from two atom type names
    Prefix = "/home/aradi/slako/mio-0-1/"  # Path as prefix
    Separator = "-"         # Dash between type names
    Suffix = ".skf"         # Suffix after second type name
  }

Historically the Slater-Koster file format did not contain any
information about which valence orbitals were considered when
generating the interaction tables, this can lead to data for
physically inappropriate orbitals being included in the files.
Therefore, you must provide the value of the highest orbital angular
momentum each element, specified as ``s``, ``p``, ``d`` or ``f``. This
information can be obtained from the documentation of the
Slater-Koster files. In the distributed standardised sets (available
at http://www.dftb.org) this information is contained in the
documentation appended to the end of each SK-file.

The default behaviour of the code is to assume that your system is
neutral (net electrical charge of 0). If you would like to calculate
charged systems, you have to use the ``Charge`` option. Similarly, the
system is assumed to be spin-unpolarised. You can however use the
option ``SpinPolarisation`` to change this standard behaviour.

The ``Filling`` option describes the method to use for filling up the
one electron levels with electrons. Here Fermi-Dirac statistics are
used. The filling functions usually requires further parameters (e.g
the temperature).

Analysis
--------

The ``Analysis`` block contains options to calculate (or display if otherwise
only calculated internally) a number of properties. In this example, while
forces are needed to optimise the geometry, these are not usually printed in
full, only the maximum value. The ``CalculateForces`` option enables printing of
the forces.

Options
-------

The ``Options`` block contains a few global settings for the code. In the
current example, no options are specified. You could even leave out the::

  Options {}

line in the input, since the default value for the ``Options`` block
is an empty block.


ParserOptions
-------------

This block contains options which are interpreted by the parser itself
and are not passed to the main program. The most important of those
options is the ``ParserVersion`` option, which tells the parser, for
which version of the parser the current input file was created for. If
this is not the current parser but an older one, the parser internally
automatically converts the old input to the new format.

The version number of the parser in the current DFTB+ code is always
printed out at the program start. It is a good habit to set this value
in your input files explicitly, like in our case::

  ParserVersion = 5

This allows you to use your input file with future versions of DFTB+
without adapting it by hand, if the input format has changed in the
more recent version.



Running DFTB+
=============


After creating the main input file, you should make sure that all the
other required files (Slater-Koster files, any files included in the
HSD input via ``<<<`` constructs, etc.) are at the right place. In our
example, only the Slater-Koster files need to be present. Since they
are specified without a path, they must be in the same directory as
the `dftb_in.hsd` file itself. This howto uses Slater-Koster files
from the `mio-0-1` SK-set.

In order to run the calculation, you should invoke DFTB+ without
any arguments in the directory containing the file `dftb_in.hsd`::

  dftb+

Assuming the binary `dftb+` lies in your search path, you should
obtain an output starting with::

  |===============================================================================
  |
  |  DFTB+ (Release 17.1)
  |
  |  Copyright (C) 2017  DFTB+ developers group
  |
  |===============================================================================
  |
  |  When publishing results obtained with DFTB+, please cite the following
  |  reference:
  |
  |  * B. Aradi, B. Hourahine and T. Frauenheim,
  |    DFTB+, a Sparse Matrix-Based Implementation of the DFTB Method,
  |    J. Phys. Chem. A, 111 5678 (2007).  [doi: 10.1021/jp070186p]
  |
  |  You should also cite additional publications crediting the parametrization
  |  data you use. Please consult the documentation of the SK-files for the
  |  references.
  |
  |===============================================================================
  
  
  ***  Parsing and initializing
  
  Parser version: 5
  
  Interpreting input file 'dftb_in.hsd'
  --------------------------------------------------------------------------------
  Reading SK-files:
    O-O.skf
    O-H.skf
    O-H.skf
    H-H.skf
  Done.
  
  
  Processed input in HSD format written to 'dftb_pin.hsd'
  
  Starting initialization...
  --------------------------------------------------------------------------------
  Mode:                        Conjugate gradient relaxation
  Self consistent charges:     Yes
  SCC-tolerance:                 0.100000E-04
  Max. scc iterations:                    100
  Ewald alpha parameter:         0.000000E+00
  Spin polarisation:           No
  Nr. of up electrons:             4.000000
  Nr. of down electrons:           4.000000
  Periodic boundaries:         No
  Diagonalizer:                Relatively robust (version 1)
  Mixer:                       Broyden mixer
  Mixing parameter:                  0.200000
  Maximal SCC-cycles:                     100
  Nr. of chrg. vec. in memory:              0
  Nr. of moved atoms:                       3
  Max. nr. of geometry steps:             100
  Force tolerance:               0.100000E-03
  Force evaluation method:     Traditional                                                                                                                                                                                             
  Electronic temperature:        0.100000E-07
  Initial charges:             Set automatically (system chrg:   0.000E+00)
  Included shells:             O:  s, p
                               H:  s
  Extra options:
                               Mulliken analysis
  Force type                   original
  
  
  --------------------------------------------------------------------------------
  
  ***  Geometry step: 0
  
      iSCC Total electronic   Diff electronic      SCC error    
      1   -0.39511797E+01    0.00000000E+00    0.88081627E+00
      2   -0.39705438E+01   -0.19364070E-01    0.55742893E+00
      3   -0.39841371E+01   -0.13593374E-01    0.32497352E-01
      4   -0.39841854E+01   -0.48242063E-04    0.19288772E-02
      5   -0.39841856E+01   -0.17020682E-06    0.87062163E-05
  
   Total Energy:                      -3.9798793068 H         -108.2980 eV
   Total Mermin free energy:          -3.9798793068 H         -108.2980 eV
   Maximal force component:            0.187090E+00
  >> Charges saved for restart in charges.bin
  
  --------------------------------------------------------------------------------
  
  ***  Geometry step: 1
  
    iSCC Total electronic   Diff electronic      SCC error    
      1   -0.40495559E+01    0.00000000E+00    0.92334735E-01
  .
  .
  . 

If this is the case, you have managed to run DFTB+ for the first
time. Congratulations!


Examining the output
====================

DFTB+ communicates through two channels with you: by printing information to
standard output (which you should probably redirect into a file to keep for
later evaluation) and by writing information into various files. In the
following, the most important of these files will be introduced and analysed


Standard output
---------------

The first thing appearing in standard output after the start of DFTB+ is the
program header::

  |===============================================================================
  |
  |  DFTB+ (Release 17.1)
  |
  |  Copyright (C) 2017  DFTB+ developers group
  |
  |===============================================================================
  |===============================================================================
  |
  |  When publishing results obtained with DFTB+, please cite the following
  |  reference:
  |
  |  * B. Aradi, B. Hourahine and T. Frauenheim,
  |    DFTB+, a Sparse Matrix-Based Implementation of the DFTB Method,
  |    J. Phys. Chem. A, 111 5678 (2007).  [doi: 10.1021/jp070186p]
  |
  |  You should also cite additional publications crediting the parametrization
  |  data you use. Please consult the documentation of the SK-files for the
  |  references.
  |
  |===============================================================================
  
  
  ***  Parsing and initializing
  
  Parser version: 5

This tells you which program you are using (DFTB+), which release (17.1) and the
paper(s) associated with the code. Then the version of the parser used in this
DFTB+ release is listed.

As already discussed above, it can be a good habit to set this version
number explicitly in your input inside the ``ParserOptions`` block,
so that::

  ParserOptions { 
    ParserVersion = 5
  }

Next, the parser starts to interpret your input, then reads in the
necessary SK-files and writes the full input settings to
`dftb_pin.hsd`::
  
  Interpreting input file 'dftb_in.hsd'
  --------------------------------------------------------------------------------
  Reading SK-files:
    O-O.skf
    O-H.skf
    O-H.skf
    H-H.skf
  Done.


  Processed input in HSD format written to 'dftb_pin.hsd'

You do not have to explicitly set all the possible options for DFTB+
in the input, as for most of them there are default values set by the
parser if not set in the input. If you want to know which default
values have been set for those missing specifications, you should look
at the processed input file `dftb_pin.hsd`, which contains the value
for all the possible input settings (see next the subsection).

At this point that the DFTB+ code is then initialised, and the most
important parameters of the calculation are then printed out::

  Mode:                        Conjugate gradient relaxation
  Self consistent charges:     Yes
  SCC-tolerance:                 0.100000E-04
  Max. scc iterations:                    100
  Ewald alpha parameter:         0.000000E+00
  Spin polarisation:           No
  Nr. of up electrons:             4.000000
  Nr. of down electrons:           4.000000
  Periodic boundaries:         No
  Diagonalizer:                Relatively robust (version 1)
  Mixer:                       Broyden mixer
  Mixing parameter:                  0.200000
  Maximal SCC-cycles:                     100
  Nr. of chrg. vec. in memory:              0
  Nr. of moved atoms:                       3
  Max. nr. of geometry steps:             100
  Force tolerance:               0.100000E-03
  Force evaluation method:     Traditional                                                                                                                                                                                             
  Electronic temperature:        0.100000E-07
  Initial charges:             Set automatically (system chrg:   0.000E+00)
  Included shells:             O:  s, p
                               H:  s
  Extra options:
                               Mulliken analysis
  Force type                   original
  

As you can see, all quantities (e.g. force tolerance, electronic
temperature) are converted to the internal units of DFTB+, namely
atomic units (with Hartree as the base energy unit).

Then the program starts::

  ***  Geometry step: 0
  
      iSCC Total electronic   Diff electronic      SCC error    
      1   -0.39511797E+01    0.00000000E+00    0.88081627E+00
      2   -0.39705438E+01   -0.19364070E-01    0.55742893E+00
      3   -0.39841371E+01   -0.13593374E-01    0.32497352E-01
      4   -0.39841854E+01   -0.48242063E-04    0.19288772E-02
      5   -0.39841856E+01   -0.17020682E-06    0.87062163E-05
  
   Total Energy:                      -3.9798793068 H         -108.2980 eV
   Total Mermin free energy:          -3.9798793068 H         -108.2980 eV
   Maximal force component:            0.187090E+00
  >> Charges saved for restart in charges.bin
  :  

Since this is an SCC calculation, DFTB+ has to iterate the charges
until the specified convergence criteria is fulfilled. In every
cycle, you get information about the values of the electronic energy,
its difference to the value in the previous SCC cycle, and the
discrepancy (error) between the charges used to build the Hamiltonian
and the charges obtained after its solution. This final value is
relevant to the tolerance specified in the input (``SCCTolerance``).

If the SCC cycle has converged, the total energy (including SCC and
repulsive contributions) is calculated, and similarly the total Mermin
free energy (this is the Helmholtz free energy, but where only the
electronic entropy is included). Additionally the biggest force
component in the system is indicated.

Then the driver changes the geometry of the system, and the
self-consistent cycle is repeated as before but for the new
geometry. This process continues as long as the geometry does not
converge::

  ***  Geometry step: 12
  
    iSCC Total electronic   Diff electronic      SCC error    
      1   -0.41505816E+01    0.00000000E+00    0.20115717E-02
      2   -0.41505816E+01   -0.21681791E-07    0.14908557E-02
      3   -0.41505816E+01   -0.26422777E-07    0.27122328E-07
  
   Total Energy:                      -4.0779379339 H         -110.9663 eV
   Total Mermin free energy:          -4.0779379339 H         -110.9663 eV
   Maximal force component:            0.280551E-05
  >> Charges saved for restart in charges.bin
  
   Geometry converged

If the geometry does not converge before the maximum number of
geometry steps is reached, the code will stop and you will get an
appropriate warning message.  Assuming the ``MaxSteps`` option had
been set to ``6`` in the input, you would obtain::

  ***  Geometry step: 6
  
    iSCC Total electronic   Diff electronic      SCC error    
      1   -0.41414806E+01    0.00000000E+00    0.12690850E-01
      2   -0.41414816E+01   -0.96478820E-06    0.93483401E-02
      3   -0.41414827E+01   -0.11442335E-05    0.17373439E-05
  
   Total Energy:                      -4.0774103506 H         -110.9520 eV
   Total Mermin free energy:          -4.0774103506 H         -110.9520 eV
   Maximal force component:            0.207962E-01
  >> Charges saved for restart in charges.bin
  WARNING!
  -> !!! Geometry did NOT converge!

dftb_pin.hsd
------------

As already mentioned, the processed input file `dftb_pin.hsd` is an input file
generated from your `dftb_in.hsd` by including the default values for all
unspecified options and converting some of the input quantities to atomic
units. For example, in our case in the ``ConjugateGradient`` block several
unspecified options would appear, for which sensible default values have been
set::

  Driver = ConjugateGradient {
    MovedAtoms = 1:-1
    MaxForceComponent = 1E-4
    MaxSteps = 100
    OutputPrefix = "geom.out"
    LatticeOpt = No
    MaxAtomStep = 0.20000000000000001
    AppendGeometries = No
    ConvergentForcesOnly = Yes
    Constraints = {}
  }

Similarly, in the ``DFTB{}`` block the switch for the orbital resolved
SCC, for example, had been set to the default value of ``No``::

  OrbitalResolvedSCC = No

Options which have been explicitly set in the input are unchanged. The file
`dftb_pin.hsd` is itself a valid HSD input file, and you can use it as input
(after renaming it to `dftb_in.hsd`) to re-run the calculation. It is always in
the format suitable for the current parser, even if the input in `dftb_in.hsd`
was for an older format (indicated by the appropriate ``ParserVersion``
option). Therefore, the ``ParserVersion`` option in the processed input file
`dftb_pin.hsd` is always set to the current version of the parser which
generated the file.


detailed.out
------------

This file contains detailed information about the properties of your
system. It is updated continuously during the run, by the end of the
calculation will contain values calculated during the last SCC
cycle. All the numerical values given in this file are in atomic
units, unless explicitly specified otherwise.

`detailed.out` contains (among other data) the number of the last
geometry step, a summary of the last SCC cycle and coordinates of any
moved atoms::

  Geometry optimization step: 12
   
  
  ********************************************************************************
    iSCC Total electronic   Diff electronic      SCC error    
      3   -0.41505816E+01   -0.26422777E-07    0.27122328E-07
  ********************************************************************************
   
   Coordinates of moved atoms (au):
      1      0.00000000     -1.35303527     -0.00000000
      2     -0.00000000     -0.26834536      1.47115110
      3      0.00000000     -0.26834536     -1.47115110

Then the net atomic charges for each atom follow (in case of |H2O|
showing a strong electron transfer from the each hydrogen to the
oxygen)::

   Net atomic charges (e)
    Atom       Net charge
       1      -0.59261515
       2       0.29630757
       3       0.29630757

.. |H2O| replace:: H\ :sub:`2`\ O
       
Then the energies of the individual electronic levels (orbitals) in
both Hartrees and electronvolts, followed by the occupation of the
individual single particle levels for all of the possible spin
channels. For spin unpolarised calculations (like this one) you will
get only one set of values, since the levels are spin restricted and
are twofold degenerate::

   Eigenvalues /H
     -0.84898606
     -0.41433754
     -0.31375444
     -0.25917545
      0.39926500
      0.55838451
   
   Eigenvalues /eV
    -23.10208606
    -11.27469810
     -8.53769263
     -7.05252282
     10.86455343
     15.19441557
   
   Fillings
       2.00000
       2.00000
       2.00000
       2.00000
       0.00000
       0.00000

In a collinear spin polarised calculation you would obtain separate
values for the spin up and spin down levels.

Then you obtain a count of the total number electrons in the system,
and the number of electrons on each atom, each atomic shell of the
atoms (s, p, d, etc.)  and each atomic orbital (labelled by their m\
:sub:`z` value) as calculated by Mulliken-analysis::

   Nr. of electrons (up):      8.00000000
   Atom populations (up)
    Atom       Population
       1       6.59261515
       2       0.70369243
       3       0.70369243
   
   l-shell populations (up)
    Atom Sh.   l       Population
       1   1   0       1.73421713
       1   2   1       4.85839802
       2   1   0       0.70369243
       3   1   0       0.70369243
   
   Orbital populations (up)
    Atom Sh.   l   m       Population
       1   1   0   0       1.73421713
       1   2   1  -1       1.68107958
       1   2   1   0       1.17731844
       1   2   1   1       2.00000000
       2   1   0   0       0.70369243
       3   1   0   0       0.70369243

In our case, due to the electronegativity difference, the hydrogen
atoms are positively charged (having only 0.704 electrons), while the
oxygen atom is negatively charged (6.59 electrons, instead of the
neutral state of 6 valence electrons).

The file then contains the Fermi energy, the different energy
contributions to the total energy and the total energy in Hartrees and
electron-volts. If you are calculating at a finite electronic
temperature, you should consider using the Mermin free energy instead
of the total energy::

   Fermi level:                        0.0700447751 H            1.9060 eV
   Band energy:                       -3.6725069692 H          -99.9340 eV
   TS:                                 0.0000000000 H            0.0000 eV
   Band free energy (E-TS):           -3.6725069692 H          -99.9340 eV
   Extrapolated E(0K):                -3.6725069692 H          -99.9340 eV
   Input/Output electrons (q):      8.00000000      8.00000000
   
   Energy H0:                         -4.1689433198 H         -113.4427 eV
   Energy SCC:                         0.0183617102 H            0.4996 eV
   Total Electronic energy:           -4.1505816095 H         -112.9431 eV
   Repulsive energy:                   0.0726436756 H            1.9767 eV
   Total energy:                      -4.0779379339 H         -110.9663 eV
   Total Mermin free energy:          -4.0779379339 H         -110.9663 eV

Between the two blocks of energy data, the input and output charges at
the last Hamiltonian diagonalisation are shown, so that you can check
that no charges get lost during the calculation.

This is then followed by a confirmation that the SCC convergence has
been reached in the last geometry step::

  SCC converged

You should always make sure that this is true, so that the properties
of your system have been calculated by using convergent
charges. Values obtained by using non convergent charges are usually
meaningless.

Finally you get the forces on the atoms in your system.  You get also
the maximal force component occurring in your system and the maximal
force occurring among the moved atoms. After this, the dipole moment of
the system (in atomic units and Debye) is printed where possible. The
end of the file will then show whether the geometry optimisation has
reached convergence, i.e., all force components on the moved atoms are
below the specified tolerance::

   Full geometry written in geom.out.{xyz|gen}
   
   Total Forces
    -1.0881602793401035E-026   6.8304105649286129E-008   4.3629613810658441E-012
    -1.9606916877279574E-016  -3.4153820160920390E-008  -2.8055131119641974E-006
     1.9606916878367734E-016  -3.4150285529999103E-008   2.8055087490097552E-006
   
   Maximal derivative component:       0.280551E-05 au
   Max force for moved atoms::         0.280551E-05 au
   
   Dipole moment  :   -0.00000000    0.64280367    0.00000000 au
   Dipole moment  :   -0.00000000    1.63384410    0.00000000 Debye
   
   Geometry converged

As indicated above, in the current case, the final relaxed geometries
can be found stored as xyz and gen format in the output files
`geom.out.xyz` and `geom.out.gen` (The package `dptools`, which
can be downloaded from the `DFTB+ website
<http://www.dftb-plus.info>`_ contains some scripts to convert between
xyz, gen and other geometry formats).


band.out
--------

For large systems, and especially for periodic systems with many
k-points, it can become quite difficult to get a good overview of the
one electron levels and their occupations in
`detailed.out`. Therefore, an extra file `band.out` is also
created, which contains this information in a more human readable
format::

  KPT            1  SPIN            1  KWEIGHT    1.0000000000000000
     -23.10209     2.00000
     -11.27470     2.00000
      -8.53769     2.00000
      -7.05252     2.00000
      10.86455     0.00000
      15.19442     0.00000

The eigenenergies are in units of electron volts. You can use the
scripts `dp_bands` in the `dptools` package to convert the data in
`band.out` to NXY-format, which can be visualised with common 2D
plotting tools.

Despite its name, the file `band.out` is also created for
non-periodic systems, containing the eigenenergies and occupation
numbers for molecular systems (You should ignore the k-point index
and the k-point weight in the first line in this case).


results.tag
-----------

If you want to process the results of DFTB+ with another program, you
should not extract the information from the standard output or the
human readable output files (`detailed.out`, `band.out`, etc.),
since their format could significantly change between subsequent
releases of DFTB+. By setting the ``WriteResultsTag`` to ``Yes`` in
the ``Options {}`` block::

  Options { 
    WriteResultsTag = Yes 
  }

you obtain the file `results.tag` at the end of your calculation,
which contains some of the most important data in a format easily
parsed by a script or a program. This file contains entries like::

  forces              :real:2:3,3
   -0.108816027934010E-025  0.683041056492861E-007  0.436296138106584E-011
   -0.196069168772796E-015 -0.341538201609204E-007 -0.280551311196420E-005
    0.196069168783677E-015 -0.341502855299991E-007  0.280550874900976E-005

In the first line the name of the quantity is given, followed by its
type (``real``, ``integer``, ``logical``). Then the rank of the
quantity is given (``0``: scalar, ``1``: vector, ``2``: rank 2 matrix,
etc.), followed by the size of each dimension. Following this, the
data for the given quantity is dumped as free format.


Other output files
------------------

There are also other output files not discussed in detail here. They are only
created, if appropriate choices in the ``Options`` or ``ExcitedState`` blocks
are set. Please consult the manual for further details.
